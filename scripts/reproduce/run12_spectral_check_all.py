"""RUN 12 — COMPLETE ESA-spectral check of ALL 107 confirmed candidates.

The pack overstated "107 spectrally confirmed" when only ~12-25 atlas events had spectra. This fetches
the ESA reduced-ion energy spectrogram (th{p}_peir_en_eflux, 32 energies) for EVERY confirmed event and
computes an AUTOMATED sheath-spectrum metric, so each of the 107 gets a per-event spectral flag.

For each event: compute radial s(t) on the substrate grid, interpolate to the ESA spectrum times, then
build the time-median ion spectrum in the near-MP shell (s in [0.05,0.20)) and in the background
(s in [0.6,1.0]). A real shocked-sheath PDL KEEPS the broad ~100s-eV sheath spectral SHAPE as s->0 (just
lower flux); magnetosphere/boundary-layer contamination shifts the peak to keV and/or changes the shape.

Metric:
  E_peak  = energy of max time-median eflux (near vs bg)
  peak_ratio = E_peak_near / E_peak_bg          (~1 sheath PDL; >>1 hot boundary/magnetosphere)
  shape_corr = Pearson r of log10 eflux(near) vs log10 eflux(bg) over channels (high => same population)
  flux_ratio = integral eflux near / bg          (depletion in the spectra)
Status:  SHEATH_CONSISTENT | HOT_BOUNDARY_FLAG | SHAPE_FLAG | AMBIGUOUS_SPEC | NO_DATA
Thresholds are PARAMETERS (printed) — to confirm at D1. Network fetch threaded (I/O releases GIL).

Usage: python run12_spectral_check_all.py [workers]   (default 8)
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from cdasws import CdasWs

SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
CONF = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run9_candidates\confirmed_pdl_candidates.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run12_spectral")
os.makedirs(OUT, exist_ok=True)
NEAR = (0.05, 0.20); BG = (0.6, 1.0); MIN_SPEC = 3
# thresholds (report; confirm at D1)
PEAK_HI = 2.5; PEAK_LO = 0.4; CORR_OK = 0.80; CORR_BAD = 0.60
_cdas = CdasWs()


def e2u(a):
    if hasattr(a, 'dtype') and np.issubdtype(a.dtype, np.datetime64):
        return (a.astype('datetime64[ns]') - np.datetime64('1970-01-01T00:00:00', 'ns')).astype(np.float64) / 1e9
    return np.asarray(a, dtype=np.float64)


def s_of_t(d):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha; r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    return np.where(denom > 0, d_mp / denom, np.nan)


def fetch_spec(probe, t0iso, t1iso):
    pl = probe[-1]; var = f"th{pl}_peir_en_eflux"
    for _ in range(3):
        try:
            st, dd = _cdas.get_data(f"TH{pl.upper()}_L2_ESA", [var], t0iso, t1iso)
            if dd is None or var not in getattr(dd, 'data_vars', {}):
                return None
            spec = np.asarray(dd[var].values, float)
            ec = [c for c in dd.coords if 'epoch' in c.lower()]
            et = e2u(dd.coords[ec[0]].values) if ec else None
            yc = [c for c in dd.coords if 'yaxis' in c.lower()]
            energy = np.asarray(dd.coords[yc[0]].values, float) if yc else None
            if energy is not None and energy.ndim == 2:
                energy = np.nanmedian(energy, axis=0)
            return et, energy, spec
        except Exception:
            pass
    return None


def median_spec(spec, mask, energy):
    """time-median eflux over masked rows; return per-channel median (Ne,) with non-finite/<=0 -> nan."""
    sub = spec[mask]
    if sub.shape[0] < MIN_SPEC:
        return None
    with np.errstate(invalid='ignore'):
        sub = np.where(np.isfinite(sub) & (sub > 0), sub, np.nan)
    med = np.nanmedian(sub, axis=0)
    return med


def classify_spec(peak_ratio, corr):
    if not (np.isfinite(peak_ratio) and np.isfinite(corr)):
        return 'NO_DATA'
    if peak_ratio > PEAK_HI:
        return 'HOT_BOUNDARY_FLAG'
    if corr < CORR_BAD:
        return 'SHAPE_FLAG'
    if PEAK_LO <= peak_ratio <= PEAK_HI and corr >= CORR_OK:
        return 'SHEATH_CONSISTENT'
    return 'AMBIGUOUS_SPEC'


def one(row):
    eid = row['eid']; probe = row['probe']
    p = os.path.join(SUB, eid + ".npz")
    if not os.path.exists(p):
        return dict(eid=eid, probe=probe, status='NO_DATA', note='no_substrate')
    d = np.load(p, allow_pickle=True)
    t = d['t'].astype(float); s = s_of_t(d)
    order = np.argsort(t); t = t[order]; s = s[order]
    sp = fetch_spec(probe, row['t_start'], row['t_end'])
    if sp is None or sp[0] is None or sp[1] is None:
        return dict(eid=eid, probe=probe, status='NO_DATA', note='no_spectra')
    et, energy, spec = sp
    if spec.ndim != 2 or et is None or len(et) < MIN_SPEC:
        return dict(eid=eid, probe=probe, status='NO_DATA', note='bad_spec_shape')
    s_et = np.interp(et, t, s, left=np.nan, right=np.nan)        # s at spectrum times
    near_m = np.isfinite(s_et) & (s_et >= NEAR[0]) & (s_et < NEAR[1])
    bg_m = np.isfinite(s_et) & (s_et >= BG[0]) & (s_et <= BG[1])
    n_near = int(near_m.sum()); n_bg = int(bg_m.sum())
    if n_near < MIN_SPEC or n_bg < MIN_SPEC:
        return dict(eid=eid, probe=probe, status='NO_DATA', note=f'insufficient near={n_near} bg={n_bg}',
                    n_near_spec=n_near, n_bg_spec=n_bg)
    mn = median_spec(spec, near_m, energy); mb = median_spec(spec, bg_m, energy)
    if mn is None or mb is None:
        return dict(eid=eid, probe=probe, status='NO_DATA', note='median_failed', n_near_spec=n_near, n_bg_spec=n_bg)
    emask = np.isfinite(energy) & (energy > 0)
    en = energy[emask]; mn2 = mn[emask]; mb2 = mb[emask]
    fin_n = np.isfinite(mn2) & (mn2 > 0); fin_b = np.isfinite(mb2) & (mb2 > 0)
    if fin_n.sum() < 3 or fin_b.sum() < 3:
        return dict(eid=eid, probe=probe, status='NO_DATA', note='too_few_channels', n_near_spec=n_near, n_bg_spec=n_bg)
    E_peak_near = float(en[fin_n][np.argmax(mn2[fin_n])])
    E_peak_bg = float(en[fin_b][np.argmax(mb2[fin_b])])
    peak_ratio = E_peak_near / E_peak_bg if E_peak_bg > 0 else np.nan
    both = fin_n & fin_b
    if both.sum() >= 3:
        corr = float(np.corrcoef(np.log10(mn2[both]), np.log10(mb2[both]))[0, 1])
    else:
        corr = np.nan
    flux_ratio = float(np.nansum(mn2[fin_n]) / np.nansum(mb2[fin_b])) if np.nansum(mb2[fin_b]) > 0 else np.nan
    status = classify_spec(peak_ratio, corr)
    return dict(eid=eid, probe=probe, status=status, note='', n_near_spec=n_near, n_bg_spec=n_bg,
                E_peak_near=round(E_peak_near, 1), E_peak_bg=round(E_peak_bg, 1),
                peak_ratio=round(peak_ratio, 3) if np.isfinite(peak_ratio) else '',
                shape_corr=round(corr, 3) if np.isfinite(corr) else '',
                flux_ratio=round(flux_ratio, 3) if np.isfinite(flux_ratio) else '')


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    rows = list(csv.DictReader(open(CONF)))
    print(f"spectral check of {len(rows)} confirmed events, {workers} threads", flush=True)
    print(f"thresholds: peak_ratio in [{PEAK_LO},{PEAK_HI}], shape_corr >= {CORR_OK} -> SHEATH_CONSISTENT", flush=True)
    res = [None] * len(rows)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(one, r): i for i, r in enumerate(rows)}
        done = 0
        from concurrent.futures import as_completed
        for fu in as_completed(futs):
            i = futs[fu]
            try:
                res[i] = fu.result()
            except Exception as ex2:
                res[i] = dict(eid=rows[i]['eid'], probe=rows[i]['probe'], status='NO_DATA', note=f'err:{ex2!r}')
            done += 1
            if done % 10 == 0 or done == len(rows):
                print(f"  {done}/{len(rows)} done", flush=True)
    import collections
    cnt = collections.Counter(r['status'] for r in res)
    cols = ['eid', 'probe', 'status', 'n_near_spec', 'n_bg_spec', 'E_peak_near', 'E_peak_bg',
            'peak_ratio', 'shape_corr', 'flux_ratio', 'note']
    path = os.path.join(OUT, "spectral_metrics.csv")
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
        for r in res:
            w.writerow(r)
    summ = [f"RUN 12 — spectral check of all {len(rows)} confirmed candidates",
            f"thresholds: SHEATH_CONSISTENT = peak_ratio in [{PEAK_LO},{PEAK_HI}] AND shape_corr >= {CORR_OK}",
            f"            HOT_BOUNDARY_FLAG = peak_ratio > {PEAK_HI} ; SHAPE_FLAG = shape_corr < {CORR_BAD}",
            "", "status counts: " + str(dict(cnt)), "", f"saved -> {path}"]
    txt = "\n".join(summ)
    print("\n" + txt, flush=True)
    with open(os.path.join(OUT, "SPECTRAL_SUMMARY.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
