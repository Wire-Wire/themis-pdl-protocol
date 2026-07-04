"""Path 2 (EXACT, re-fetch): full archive dual-geometry recompute.

Mirrors build_themis_cache_v2.process_enc exactly (same windowing, cadence=MOM native,
fill masks, occupancy thresholds) but computes BOTH s_1D and s_radial per sample, so
Dn_1d reproduces the cached catalogue (validation) and Dn_rad is the corrected value on
identical data. Positions from cached raw_state; MOM+FGM+OMNI re-fetched (cdasws).
Per-encounter results cached -> resumable. Usage: python recompute_archive_radial.py [limit]
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import json, os, sys, time, datetime as dt_module, csv
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
sys.path.insert(0, P(r"H:\0mssl\0mssl519\src"))
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from cdasws import CdasWs
from pdl_pilot.mapping.s_mapper import compute_s_with_uncertainty, compute_bin_occupancy, compute_s_radial
from pdl_pilot.boundaries.shue1998 import shue1998_alpha
from pdl_pilot.config.schema import BinConfig, UncertaintyConfig
from radial_models import compute_s_radial_jbs

BINS = BinConfig(); UNC = UncertaintyConfig()
X_MIN = 5.0; SZA_MAX = 30.0; R_MIN = 8.0; R_MAX = 25.0; ENC_HOURS = 6; MAX_RETRIES = 3
WORKERS = 16
BASE = P(r"H:\0mssl\0mssl519\data_cache\themis_archive")
STATE_DIR = os.path.join(BASE, "raw_state")
CACHE_REF = os.path.join(BASE, "encounters")          # original 1D results (validation ref)
OUT_DIR = os.path.join(BASE, "encounters_radial2")
os.makedirs(OUT_DIR, exist_ok=True)
NEAR_VARIANTS = [(0.0, 0.05), (0.0, 0.1), (0.1, 0.2), (0.2, 0.4)]
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 0   # 0 = all

_tls = threading.local()
def get_cdas():
    if not hasattr(_tls, 'cdas'): _tls.cdas = CdasWs()
    return _tls.cdas
_lock = threading.Lock()
def tprint(m):
    with _lock: print(m, flush=True)

def _e2u(a):
    if hasattr(a, 'dtype') and np.issubdtype(a.dtype, np.datetime64):
        return (a.astype('datetime64[ns]') - np.datetime64('1970-01-01T00:00:00', 'ns')).astype(np.float64) / 1e9
    return np.asarray(a, dtype=np.float64)

def _interp(st, sd, tt, mg=60.0):
    if len(st) == 0 or len(sd) == 0: return np.full_like(tt, np.nan)
    r = np.interp(tt, st, sd.astype(np.float64))
    if mg > 0 and len(st) > 1:
        i = np.clip(np.searchsorted(st, tt), 0, len(st) - 1); ip = np.clip(i - 1, 0, len(st) - 1)
        d = np.minimum(np.abs(tt - st[i]), np.abs(tt - st[ip])); r[d > mg] = np.nan
    return r

def _fetch(ds, vs, t0, t1):
    for a in range(MAX_RETRIES + 1):
        try:
            _, d = get_cdas().get_data(ds, vs, t0, t1); return d
        except Exception:
            if a < MAX_RETRIES: time.sleep(1 + a * 2)
    return None

def _med(a, m):
    v = a[m]; v = v[np.isfinite(v)]
    return float(np.nanmedian(v)) if len(v) > 0 else None

def ratios(s, den, bm, near):
    nm = (s >= near[0]) & (s < near[1]); bgm = (s >= 0.6) & (s <= 1.0)
    nn = _med(den, nm); nb = _med(den, bgm); bn = _med(bm, nm); bb = _med(bm, bgm)
    Dn = (nn / nb) if (nn and nb and nb > 0) else None
    EB = (bn / bb) if (bn and bb and bb > 0) else None
    return Dn, EB, float(nm.mean())

# ---- screen raw_state -> qualifying days (mirror build Phase 2) ----
def screen_one(fname):
    parts = fname.replace('.npz', '').split('_')
    if len(parts) != 3: return []
    probe, ys, ms = parts
    try:
        d = np.load(os.path.join(STATE_DIR, fname), allow_pickle=True)
        t = d["times"]; x = d["x_re"]; y = d["y_re"]; z = d["z_re"]
    except Exception:
        return []
    if len(t) == 0: return []
    r = np.sqrt(x**2 + y**2 + z**2); sza = np.degrees(np.arctan2(np.sqrt(y**2 + z**2), x))
    good = (x > X_MIN) & (sza < SZA_MAX) & (r > R_MIN) & (r < R_MAX) & np.isfinite(x)
    if not np.any(good): return []
    wins = {}
    for i in range(len(t)):
        if not good[i]: continue
        dk = dt_module.datetime.fromtimestamp(t[i], tz=dt_module.timezone.utc).date()
        if dk not in wins or sza[i] < wins[dk][1]: wins[dk] = (float(t[i]), float(sza[i]))
    return [(str(dk), probe, ct, sz, int(ys), int(ms)) for dk, (ct, sz) in wins.items()]

def process(q):
    ds, probe, ct, _, year, month = q; pl = probe[-1]; eid = f"{ds}_{probe}"
    op = os.path.join(OUT_DIR, f"{eid}.json")
    if os.path.exists(op):
        try:
            return json.load(open(op))
        except Exception:
            os.remove(op)
    half = ENC_HOURS * 3600 / 2
    ts = dt_module.datetime.fromtimestamp(ct - half, tz=dt_module.timezone.utc)
    te = dt_module.datetime.fromtimestamp(ct + half, tz=dt_module.timezone.utc)
    t0s = ts.strftime("%Y-%m-%dT%H:%M:%SZ"); t1s = te.strftime("%Y-%m-%dT%H:%M:%SZ")
    base = {"encounter_id": eid, "date": ds, "probe": probe, "year": year, "month": month,
            "evaluable_1d": False, "evaluable_rad": False, "exclude_reason": None}
    try:
        sp = os.path.join(STATE_DIR, f"{probe}_{year}_{month:02d}.npz")
        if not os.path.exists(sp): base["exclude_reason"] = "no STATE"; return base
        sd = np.load(sp, allow_pickle=True); st = sd["times"]; t0 = ct - half - 600; t1 = ct + half + 600
        m = (st >= t0) & (st <= t1); stt = st[m]; stx = sd["x_re"][m]; sty = sd["y_re"][m]; stz = sd["z_re"][m]
        mom = _fetch(f"TH{pl.upper()}_L2_MOM", [f"th{pl}_peim_density"], t0s, t1s)
        if mom is None or f"th{pl}_peim_density" not in mom.data_vars: base["exclude_reason"] = "no MOM"; return base
        mc = [c for c in mom.coords if "epoch" in c.lower() or "time" in c.lower()]
        if not mc: base["exclude_reason"] = "no MOM time"; return base
        mt = _e2u(mom.coords[mc[0]].values); den = mom[f"th{pl}_peim_density"].values.astype(np.float64)
        if len(mt) < 10: base["exclude_reason"] = "insufficient MOM"; return base
        tt = mt
        fgm = _fetch(f"TH{pl.upper()}_L2_FGM", [f"th{pl}_fgs_gsm"], t0s, t1s)
        if fgm is None or f"th{pl}_fgs_gsm" not in fgm.data_vars: base["exclude_reason"] = "no FGM"; return base
        bv = fgm[f"th{pl}_fgs_gsm"].values.astype(np.float64)
        if bv.ndim < 2: base["exclude_reason"] = "FGM not 3c"; return base
        bmr = np.sqrt(np.sum(bv**2, axis=-1))
        fc = [c for c in fgm.coords if "epoch" in c.lower() or "time" in c.lower()]
        ft = _e2u(fgm.coords[fc[0]].values) if fc else None
        if ft is None or len(ft) != len(bmr):
            for d in fgm[f"th{pl}_fgs_gsm"].dims:
                if d in fgm.coords and len(fgm.coords[d].values) == len(bmr): ft = _e2u(fgm.coords[d].values); break
        if ft is None or len(ft) != len(bmr):
            ft = np.linspace(tt[0], tt[-1], len(bmr)) if len(bmr) > 0 else None
        if ft is None: base["exclude_reason"] = "FGM time fail"; return base
        bm = _interp(ft, bmr, tt, 30.0)
        xr = _interp(stt, stx, tt, 600.0); yr = _interp(stt, sty, tt, 600.0); zr = _interp(stt, stz, tt, 600.0)
        xm = float(np.nanmean(xr)); ym = float(np.nanmean(yr)); zm = float(np.nanmean(zr))
        sza = float(np.degrees(np.arctan2(np.sqrt(ym**2 + zm**2), xm))) if np.isfinite(xm) else 999.
        if sza > SZA_MAX: base["sza_deg"] = round(sza, 1); base["exclude_reason"] = f"SZA={sza:.0f}"; return base
        os0 = (ts - dt_module.timedelta(minutes=30)).strftime("%Y-%m-%dT%H:%M:%SZ")
        oe0 = (te + dt_module.timedelta(minutes=10)).strftime("%Y-%m-%dT%H:%M:%SZ")
        om = _fetch("OMNI_HRO_1MIN", ["BX_GSE", "BY_GSM", "BZ_GSM", "F", "Pressure", "Mach_num"], os0, oe0)
        dp = bz = bt = bx = by = ma = None
        if om is not None:
            def _o(k):
                if k not in om.data_vars: return None
                a = om[k].values.astype(np.float64); v = a[np.isfinite(a)]; v = v[(v > -1e30) & (v < 1e30) & (np.abs(v) < 9990)]
                return float(np.nanmedian(v)) if len(v) > 0 else None
            dp = _o("Pressure"); bz = _o("BZ_GSM"); bt = _o("F"); ma = _o("Mach_num"); bx = _o("BX_GSE"); by = _o("BY_GSM")
        if dp is None or dp <= 0: base["sza_deg"] = round(sza, 1); base["exclude_reason"] = "no OMNI Dp"; return base
        cd = float(np.degrees(np.arccos(min(abs(bx) / bt, 1.0)))) if bx is not None and bt and bt > 0 else None
        ck = (float(np.degrees(np.arctan2(by, bz))) % 360.0) if by is not None and bz is not None else None
        den[den <= 0] = np.nan; den[den > 1000] = np.nan; bm[bm <= 0] = np.nan; bm[bm > 500] = np.nan
        # ---- 1D geometry (reproduce) ----
        s1, _, _, mp, bs = compute_s_with_uncertainty(xr, dp_nPa=dp, bz_nT=bz or 0., mach_alfven=ma or 8., unc=UNC)
        occ1 = compute_bin_occupancy(s1, BINS)
        Dn1, EB1, _ = ratios(s1, den, bm, (0.2, 0.4))
        # ---- radial geometry (same-alpha BS placeholder) ----
        alpha = shue1998_alpha(dp, bz or 0.0)
        srad = compute_s_radial(xr, yr, zr, mp, bs, alpha)
        occr = compute_bin_occupancy(srad, BINS)
        Dnr, EBr, _ = ratios(srad, den, bm, (0.2, 0.4))
        # ---- radial geometry with proper Jelinek-2012 bow shock (Shue MP + Jelinek BS) ----
        sjbs = compute_s_radial_jbs(xr, yr, zr, mp, alpha, dp)
        occj = compute_bin_occupancy(sjbs, BINS)
        Dnj, EBj, _ = ratios(sjbs, den, bm, (0.2, 0.4))
        rm = float(np.sqrt(xm**2 + ym**2 + zm**2))
        cb = ("quasi-radial" if cd < 30 else "low-cone" if cd <= 45 else "intermediate" if cd <= 60 else "perpendicular") if cd is not None else "unknown"
        ev1 = occ1["near"] >= .05 and occ1["background"] >= .01 and Dn1 is not None and EB1 is not None
        evr = occr["near"] >= .05 and occr["background"] >= .01 and Dnr is not None and EBr is not None
        evj = occj["near"] >= .05 and occj["background"] >= .01 and Dnj is not None and EBj is not None
        e = {"encounter_id": eid, "date": ds, "probe": probe, "year": year, "month": month,
             "sza_deg": round(sza, 1), "cone_deg": round(cd, 1) if cd is not None else None,
             "clock_deg": round(ck, 1) if ck is not None else None, "cone_bin": cb,
             "dp_nPa": round(dp, 2), "bz_nT": round(bz, 2) if bz else None, "ma": round(ma, 1) if ma else None,
             "mp_standoff_re": round(mp, 2), "bs_standoff_re": round(bs, 2), "alpha": round(alpha, 4),
             "x_gsm_re": round(xm, 2), "r_re": round(rm, 2), "n_points": len(tt),
             "near_occ_1d": round(occ1["near"], 4), "bg_occ_1d": round(occ1["background"], 4),
             "Dn_1d": round(Dn1, 4) if Dn1 else None, "EB_1d": round(EB1, 4) if EB1 else None, "evaluable_1d": ev1,
             "near_occ_rad": round(occr["near"], 4), "bg_occ_rad": round(occr["background"], 4),
             "Dn_rad": round(Dnr, 4) if Dnr else None, "EB_rad": round(EBr, 4) if EBr else None, "evaluable_rad": evr,
             "near_occ_jbs": round(occj["near"], 4), "bg_occ_jbs": round(occj["background"], 4),
             "Dn_jbs": round(Dnj, 4) if Dnj else None, "EB_jbs": round(EBj, 4) if EBj else None, "evaluable_jbs": evj,
             "s_shift_mean": round(float(np.nanmean(np.abs(srad - s1))), 4)}
        for ne in NEAR_VARIANTS:
            dv, ev_, ov = ratios(srad, den, bm, ne)
            e[f"Dn_rad_{ne[0]}_{ne[1]}"] = round(dv, 4) if dv else None
            e[f"occ_rad_{ne[0]}_{ne[1]}"] = round(ov, 4)
        with open(op, 'w') as f: json.dump(e, f, default=str)
        return e
    except Exception as ex:
        base["exclude_reason"] = f"err:{ex}"; return base

if __name__ == "__main__":
    t_start = time.time()
    state_files = sorted([f for f in os.listdir(STATE_DIR) if f.endswith('.npz')])
    print(f"raw_state files: {len(state_files)}")
    qual = []
    for f in state_files: qual.extend(screen_one(f))
    seen = {}
    for ds, probe, ct, sz, y, m in qual:
        k = (ds, probe)
        if k not in seen or sz < seen[k][3]: seen[k] = (ds, probe, ct, sz, y, m)
    qual = sorted(seen.values(), key=lambda v: (v[0], v[1]))
    print(f"qualifying days: {len(qual)}")
    if LIMIT > 0:
        # validation sample: spread across the list
        step = max(1, len(qual) // LIMIT); qual = qual[::step][:LIMIT]
        print(f"LIMIT={LIMIT} -> processing {len(qual)} sampled encounters")
    done = [0]; retr = [0]
    cat = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process, q): q for q in qual}
        for fut in as_completed(futs):
            done[0] += 1
            try:
                e = fut.result(); cat.append(e)
                if e.get("evaluable_rad") or e.get("evaluable_1d"): retr[0] += 1
            except Exception:
                pass
            if done[0] % 200 == 0:
                tprint(f"  [{done[0]}/{len(qual)}] processed, {retr[0]} evaluable(1d|rad), {time.time()-t_start:.0f}s")
    # write combined catalogue
    OUTCSV = P(r"H:\0mssl\review\repair\option3\derived\archive_radial_catalogue_v2.csv")
    os.makedirs(os.path.dirname(OUTCSV), exist_ok=True)
    cols = ["encounter_id", "date", "probe", "year", "month", "sza_deg", "cone_deg", "clock_deg", "cone_bin",
            "dp_nPa", "bz_nT", "ma", "mp_standoff_re", "bs_standoff_re", "alpha", "x_gsm_re", "r_re", "n_points",
            "near_occ_1d", "bg_occ_1d", "Dn_1d", "EB_1d", "evaluable_1d",
            "near_occ_rad", "bg_occ_rad", "Dn_rad", "EB_rad", "evaluable_rad",
            "near_occ_jbs", "bg_occ_jbs", "Dn_jbs", "EB_jbs", "evaluable_jbs", "s_shift_mean"]
    for ne in NEAR_VARIANTS:
        cols += [f"Dn_rad_{ne[0]}_{ne[1]}", f"occ_rad_{ne[0]}_{ne[1]}"]
    full = [e for e in cat if "Dn_1d" in e or "evaluable_1d" in e]
    with open(OUTCSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
        for e in cat: w.writerow(e)
    n1 = sum(1 for e in cat if e.get("evaluable_1d")); nr = sum(1 for e in cat if e.get("evaluable_rad"))
    print(f"\nDONE {time.time()-t_start:.0f}s. processed={len(cat)} evaluable_1d={n1} evaluable_rad={nr}")
    print("WROTE", OUTCSV)
