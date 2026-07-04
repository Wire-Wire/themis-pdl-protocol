"""FROZEN ANALYSIS SUBSTRATE BUILDER (32-worker, I/O-parallel).

One-time re-fetch that CACHES the per-sample time-series for every archive encounter,
so all campaign runs (profile cube, bin-multiverse, favourable-conditions, HQ rebuild,
confounder screen) become network-free LOCAL re-aggregations of a single frozen source.

Per encounter it stores (themis_archive/substrate/<eid>.npz):
  arrays on the MOM timeline: t, x_re, y_re, z_re, n (cm^-3), bmag (nT), beta, p_th (nPa),
    p_b (nPa), p_tot (nPa), vmag (km/s)
  scalars: dp, bz, ma, bx, cone_deg, clock_deg, sza_deg, mp0, bs0, alpha, r_re, x_mean, probe, date, year, month
Any geometry's s and any bin's Dn/EB/beta-ratio/pressure-ratio are then computed locally.
Resumable (skips existing npz), atomic writes. Disjoint output dir from the running v2 job.
Usage: python build_substrate.py [limit]
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
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
WORKERS = 32                      # 32 logical cores; I/O-bound so threads>=cores is fine
MU0 = 4 * np.pi * 1e-7
BASE = P(r"H:\0mssl\0mssl519\data_cache\themis_archive")
STATE_DIR = os.path.join(BASE, "raw_state")
SUB_DIR = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")   # relocated into review for self-containment
os.makedirs(SUB_DIR, exist_ok=True)
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 0

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
    op = os.path.join(SUB_DIR, f"{eid}.npz")
    if os.path.exists(op):
        return ("cached", eid)
    half = ENC_HOURS * 3600 / 2
    ts = dt_module.datetime.fromtimestamp(ct - half, tz=dt_module.timezone.utc)
    te = dt_module.datetime.fromtimestamp(ct + half, tz=dt_module.timezone.utc)
    t0s = ts.strftime("%Y-%m-%dT%H:%M:%SZ"); t1s = te.strftime("%Y-%m-%dT%H:%M:%SZ")
    try:
        sp = os.path.join(STATE_DIR, f"{probe}_{year}_{month:02d}.npz")
        if not os.path.exists(sp): return ("skip", eid + ":no_state")
        sd = np.load(sp, allow_pickle=True); st = sd["times"]; t0 = ct - half - 600; t1 = ct + half + 600
        m = (st >= t0) & (st <= t1); stt = st[m]; stx = sd["x_re"][m]; sty = sd["y_re"][m]; stz = sd["z_re"][m]
        # MOM: density + total pressure + velocity (for beta/pressure/jet diagnostics)
        mom = _fetch(f"TH{pl.upper()}_L2_MOM",
                     [f"th{pl}_peim_density", f"th{pl}_peim_ptot", f"th{pl}_peim_velocity_gsm"], t0s, t1s)
        if mom is None or f"th{pl}_peim_density" not in mom.data_vars: return ("skip", eid + ":no_mom")
        mc = [c for c in mom.coords if "epoch" in c.lower() or "time" in c.lower()]
        if not mc: return ("skip", eid + ":no_mom_time")
        tt = _e2u(mom.coords[mc[0]].values)
        n = mom[f"th{pl}_peim_density"].values.astype(np.float64)
        if len(tt) < 10: return ("skip", eid + ":short_mom")
        ptot_ev = mom[f"th{pl}_peim_ptot"].values.astype(np.float64) if f"th{pl}_peim_ptot" in mom.data_vars else np.full(len(tt), np.nan)
        if f"th{pl}_peim_velocity_gsm" in mom.data_vars and mom[f"th{pl}_peim_velocity_gsm"].values.ndim >= 2:
            vmag = np.sqrt(np.sum(mom[f"th{pl}_peim_velocity_gsm"].values.astype(np.float64) ** 2, axis=-1))
        else:
            vmag = np.full(len(tt), np.nan)
        # FGM |B| -> interp to MOM timeline
        fgm = _fetch(f"TH{pl.upper()}_L2_FGM", [f"th{pl}_fgs_gsm"], t0s, t1s)
        if fgm is None or f"th{pl}_fgs_gsm" not in fgm.data_vars: return ("skip", eid + ":no_fgm")
        bv = fgm[f"th{pl}_fgs_gsm"].values.astype(np.float64)
        if bv.ndim < 2: return ("skip", eid + ":fgm_2d")
        bmr = np.sqrt(np.sum(bv ** 2, axis=-1))
        fc = [c for c in fgm.coords if "epoch" in c.lower() or "time" in c.lower()]
        ft = _e2u(fgm.coords[fc[0]].values) if fc else None
        if ft is None or len(ft) != len(bmr):
            for d in fgm[f"th{pl}_fgs_gsm"].dims:
                if d in fgm.coords and len(fgm.coords[d].values) == len(bmr): ft = _e2u(fgm.coords[d].values); break
        if ft is None or len(ft) != len(bmr):
            ft = np.linspace(tt[0], tt[-1], len(bmr)) if len(bmr) > 0 else None
        if ft is None: return ("skip", eid + ":fgm_time")
        bmag = _interp(ft, bmr, tt, 30.0)
        xr = _interp(stt, stx, tt, 600.0); yr = _interp(stt, sty, tt, 600.0); zr = _interp(stt, stz, tt, 600.0)
        xm = float(np.nanmean(xr)); ym = float(np.nanmean(yr)); zm = float(np.nanmean(zr))
        if not np.isfinite(xm): return ("skip", eid + ":no_pos")
        sza = float(np.degrees(np.arctan2(np.sqrt(ym ** 2 + zm ** 2), xm)))
        if sza > SZA_MAX: return ("skip", eid + f":sza{sza:.0f}")
        # OMNI
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
        if dp is None or dp <= 0: return ("skip", eid + ":no_dp")
        cone = float(np.degrees(np.arccos(min(abs(bx) / bt, 1.0)))) if bx is not None and bt and bt > 0 else np.nan
        clock = (float(np.degrees(np.arctan2(by, bz))) % 360.0) if by is not None and bz is not None else np.nan
        # clean + derive pressures/beta
        n[(n <= 0) | (n > 1000)] = np.nan
        bmag[(bmag <= 0) | (bmag > 500)] = np.nan
        p_b = (bmag * 1e-9) ** 2 / (2 * MU0) * 1e9                  # nPa
        p_th = ptot_ev * 1.6e-4                                      # eV/cm^3 -> nPa
        p_tot = p_th + p_b
        beta = np.where(p_b > 0, p_th / p_b, np.nan)
        mp0 = float(__import__('pdl_pilot.boundaries.shue1998', fromlist=['shue1998_standoff']).shue1998_standoff(dp, bz or 0.0))
        bs0 = float(__import__('pdl_pilot.boundaries.merka2005', fromlist=['merka2005_standoff']).merka2005_standoff(mp0, ma or 8.0))
        alpha = shue1998_alpha(dp, bz or 0.0)
        rm = float(np.sqrt(xm ** 2 + ym ** 2 + zm ** 2))
        tmp = op + f".tmp{threading.get_ident()}.npz"   # MUST end in .npz: np.savez_compressed appends .npz otherwise, breaking os.replace
        np.savez_compressed(tmp,
            t=tt.astype(np.float64), x_re=xr.astype(np.float32), y_re=yr.astype(np.float32), z_re=zr.astype(np.float32),
            n=n.astype(np.float32), bmag=bmag.astype(np.float32), beta=beta.astype(np.float32),
            p_th=p_th.astype(np.float32), p_b=p_b.astype(np.float32), p_tot=p_tot.astype(np.float32),
            vmag=vmag.astype(np.float32),
            dp=dp, bz=(bz if bz is not None else np.nan), ma=(ma if ma is not None else np.nan),
            bx=(bx if bx is not None else np.nan), cone_deg=cone, clock_deg=clock, sza_deg=sza,
            mp0=mp0, bs0=bs0, alpha=alpha, r_re=rm, x_mean=xm,
            probe=probe, date=ds, year=year, month=month)
        os.replace(tmp, op)
        return ("ok", eid)
    except Exception as ex:
        return ("err", f"{eid}:{ex}")

if __name__ == "__main__":
    t_start = time.time()
    qual = []
    for f in sorted(os.listdir(STATE_DIR)):
        if f.endswith('.npz'): qual.extend(screen_one(f))
    seen = {}
    for ds, probe, ct, sz, y, mth in qual:
        k = (ds, probe)
        if k not in seen or sz < seen[k][3]: seen[k] = (ds, probe, ct, sz, y, mth)
    qual = sorted(seen.values(), key=lambda v: (v[0], v[1]))
    print(f"qualifying days: {len(qual)}  (WORKERS={WORKERS})", flush=True)
    if LIMIT > 0:
        step = max(1, len(qual) // LIMIT); qual = qual[::step][:LIMIT]
        print(f"LIMIT={LIMIT} -> {len(qual)} sampled", flush=True)
    done = [0]; ok = [0]; cached = [0]; skip = [0]; err = [0]; samples = []; skips = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process, q): q for q in qual}
        for fut in as_completed(futs):
            done[0] += 1
            try:
                tag, info = fut.result()
                if tag == "ok": ok[0] += 1
                elif tag == "cached": cached[0] += 1
                elif tag == "skip": skip[0] += 1; skips.append(info)
                elif tag == "err":
                    err[0] += 1
                    if len(samples) < 12: samples.append(info)
            except Exception as _e:
                err[0] += 1
                if len(samples) < 12: samples.append(f"future:{_e}")
            if done[0] % 200 == 0:
                tprint(f"  [{done[0]}/{len(qual)}] ok={ok[0]} cached={cached[0]} skip={skip[0]} err={err[0]} {time.time()-t_start:.0f}s")
    n_files = len([f for f in os.listdir(SUB_DIR) if f.endswith('.npz')])
    print(f"\nSUBSTRATE DONE {time.time()-t_start:.0f}s. processed={done[0]} recovered_new_ok={ok[0]} cached={cached[0]} skip={skip[0]} err={err[0]} total_npz={n_files}", flush=True)
    if samples: print("SAMPLE ISSUES:", " || ".join(str(s) for s in samples), flush=True)
    # Categorise skips: reason = text after last ':' with digits normalised (sza28 -> szaN).
    # network-ish reasons (no_mom/no_fgm/no_dp) = CDAWeb returned nothing after 3 retries (likely no data, possibly net);
    # filter reasons (szaN/short_mom/no_state/no_pos) = legitimately excluded.
    import re as _re
    from collections import Counter as _Counter
    def _catreason(info):
        r = info.rsplit(":", 1)[-1] if ":" in info else info
        return _re.sub(r"\d+", "N", r)
    _counts = dict(sorted(_Counter(_catreason(s) for s in skips).items(), key=lambda kv: -kv[1]))
    print("SKIP BREAKDOWN:", _counts, flush=True)
    print(f"RECOVERED this run (were transient/network) = {ok[0]}", flush=True)
    _rep = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate_skip_report.txt")
    with open(_rep, "w") as _f:
        _f.write("SKIP REPORT (re-run %s)\n" % time.strftime("%Y-%m-%d %H:%M"))
        _f.write("recovered_this_run (transient/network, now cached) = %d\n" % ok[0])
        _f.write("persistent skips = %d\n" % skip[0])
        _f.write("category counts = %s\n\n" % _counts)
        for s in sorted(skips):
            _f.write(s + "\n")
    print("skip report ->", _rep, flush=True)
