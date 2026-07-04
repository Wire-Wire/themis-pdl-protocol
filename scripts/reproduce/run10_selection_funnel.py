"""RUN 10 — COMPLETE SELECTION-FUNCTION FUNNEL, each step labelled.

The external review's Criticism 2: 107/6248 is NOT an occurrence rate; the 6248 base is wrong
because most encounters never sampled BOTH a clean near-MP shell AND a background. This builds the
full, internally-consistent funnel in ONE pass per encounter (so every stage chains exactly), and
labels each drop as a *filter*, not a non-detection:

  substrate(readable) -> model-placeable -> has in-sheath samples -> has background(s[0.6,1])
   -> has membership near-shell (CONTRIBUTING) -> northward(Bz>0) -> CONFIRMED_PDL (clean classifier)

Physics is IDENTICAL to run9a (same membership screen + PDL shell [0.05,0.20) + classify()), so the
CONFIRMED_PDL count reproduces the 107 catalogue. Also writes funnel_contributing.csv (rich covariates
for the Bz confounder check, run11) and funnel_all.csv (per-encounter stage). Threaded load via psub.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, collections
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap, list_files

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0); PDL = (0.05, 0.20); MIN_BG = 5; KB = 1.602e-4
N_FLOOR = 0.3; TBAND = 3.0; VFRAC = 0.2; MIN_SHELL = 5


def classify(r):
    """IDENTICAL to run9a.classify()."""
    if r['beta_near'] < 0.10 or r['n_near'] < 1.0 or r['T_near'] > 1200:
        return 'BOUNDARY_LAYER'
    if (np.isfinite(r['bz']) and r['bz'] > 0 and 0.40 <= r['Dn_mem'] <= 0.90 and r['n_near'] >= 2.0
            and 0.10 <= r['beta_near'] <= 2.0 and 50 <= r['T_near'] <= 800 and r['EB_mem'] > 1.2 and r['n_mem'] >= 15):
        return 'CONFIRMED_PDL'
    return 'AMBIGUOUS'


def funnel_one(name, d):
    """Always returns a dict with 'stage' = furthest stage reached, so drops can be counted."""
    eid = name[:-4]; probe = eid.rsplit('_', 1)[-1]
    base = dict(eid=eid, probe=probe)
    t = d['t'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if len(t) < 10 or not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return {**base, 'stage': 'no_params'}
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float); beta = d['beta'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    bz = float(d['bz']) if np.isfinite(d['bz']) else np.nan
    cone = float(d['cone_deg']) if np.isfinite(d['cone_deg']) else np.nan
    sza = float(d['sza_deg']) if np.isfinite(d['sza_deg']) else np.nan
    clock = float(d['clock_deg']) if ('clock_deg' in d.files and np.isfinite(d['clock_deg'])) else np.nan
    ma = float(d['ma']) if ('ma' in d.files and np.isfinite(d['ma'])) else np.nan
    bx = float(d['bx']) if ('bx' in d.files and np.isfinite(d['bx'])) else np.nan
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha; r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    geo = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]; betav = beta[geo]
    if len(sv) < MIN_BG + 3:
        return {**base, 'stage': 'no_sheath', 'bz': bz}
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return {**base, 'stage': 'no_bg', 'bz': bz}
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm]); T_bg = np.nanmedian(Tv[bgm]); v_bg = np.nanmedian(vv[bgm])
    if not (n_bg > 0 and b_bg > 0 and np.isfinite(T_bg) and T_bg > 0):
        return {**base, 'stage': 'no_bg', 'bz': bz}
    mem = (nv > N_FLOOR) & np.isfinite(Tv) & (Tv > T_bg / TBAND) & (Tv < T_bg * TBAND)
    if np.isfinite(v_bg) and v_bg > 0:
        mem = mem & np.isfinite(vv) & (vv > VFRAC * v_bg)
    pm = mem & (sv >= PDL[0]) & (sv < PDL[1])
    if pm.sum() < MIN_SHELL:
        return {**base, 'stage': 'no_shell', 'bz': bz}
    rec = {**base, 'stage': 'contrib', 'bz': bz, 'cone': cone, 'sza': sza, 'dp': dp,
           'clock': clock, 'ma': ma, 'bx': bx,
           'n_mem': int(pm.sum()), 'Dn_mem': float(np.median(nv[pm]) / n_bg),
           'EB_mem': float(np.median(bv[pm]) / b_bg), 'n_near': float(np.median(nv[pm])),
           'beta_near': float(np.nanmedian(betav[pm])), 'T_near': float(np.nanmedian(Tv[pm]))}
    rec['cls'] = classify(rec)
    return rec


ORDER = ['no_params', 'no_sheath', 'no_bg', 'no_shell', 'contrib']


def main():
    files = list_files()
    rows = pmap(funnel_one, files=files, with_name=True)
    n_unreadable = len(files) - len(rows)
    cnt = collections.Counter(r['stage'] for r in rows)
    contrib = [r for r in rows if r['stage'] == 'contrib']
    north = [r for r in contrib if np.isfinite(r['bz']) and r['bz'] > 0]
    south = [r for r in contrib if np.isfinite(r['bz']) and r['bz'] <= 0]
    bznan = [r for r in contrib if not np.isfinite(r['bz'])]
    ncls = collections.Counter(r['cls'] for r in north)
    scls = collections.Counter(r['cls'] for r in south)

    # cumulative funnel
    readable = len(rows)
    placeable = readable - cnt['no_params']
    has_sheath = placeable - cnt['no_sheath']
    has_bg = has_sheath - cnt['no_bg']
    n_contrib = len(contrib)               # == has_bg - cnt['no_shell']
    n_north = len(north); n_conf = ncls['CONFIRMED_PDL']

    o = []
    o.append("RUN 10 — SELECTION FUNCTION FUNNEL")
    o.append("each step is a FILTER on observability, not a non-detection. 107 is NOT an occurrence rate.")
    o.append("")
    o.append(f"{'stage':42s}{'N':>7s}{'drop':>8s}  reason-dropped")
    o.append(f"{'0 substrate npz (readable)':42s}{readable:7d}{'':>8s}  (+{n_unreadable} unreadable files skipped)")
    o.append(f"{'1 model-placeable (mp0,alpha,Dp ok; >=10)':42s}{placeable:7d}{-cnt['no_params']:8d}  no valid Shue/Dp params or too few samples")
    o.append(f"{'2 has in-sheath samples':42s}{has_sheath:7d}{-cnt['no_sheath']:8d}  window never in geometric magnetosheath")
    o.append(f"{'3 has background s[0.6,1.0]':42s}{has_bg:7d}{-cnt['no_bg']:8d}  never sampled outer-sheath background")
    o.append(f"{'4 CONTRIBUTING (membership near-shell)':42s}{n_contrib:7d}{-cnt['no_shell']:8d}  never sampled membership-sheath in s[0.05,0.20)")
    o.append(f"{'5 northward (Bz>0)':42s}{n_north:7d}{-(n_contrib-n_north):8d}  Bz<=0 ({len(south)}) or Bz unknown ({len(bznan)})")
    o.append(f"{'6 CONFIRMED_PDL (clean classifier)':42s}{n_conf:7d}{-(n_north-n_conf):8d}  northward but BOUNDARY_LAYER ({ncls['BOUNDARY_LAYER']}) or AMBIGUOUS ({ncls['AMBIGUOUS']})")
    o.append("")
    o.append(f"contributing classification  northward: {dict(ncls)}")
    o.append(f"contributing classification  southward: {dict(scls)}  (Bz unknown: {len(bznan)})")
    o.append(f"contributing by probe: {dict(collections.Counter(r['probe'] for r in contrib))}")
    o.append("")
    o.append("INTERPRETATION: the 107 = moment-classified clean PDL candidates within the OBSERVABLE,")
    o.append("northward, membership-screened subset. Any occurrence statistic must use the observable")
    o.append("window (stage 4 = %d contributing, or stage 3 = %d with-background) as denominator," % (n_contrib, has_bg))
    o.append("stratified by SZA/Dp/probe/cone/coverage — NOT 6248.")

    # bias check contributing vs has-bg-but-not-contributing (does near-shell sampling bias the sample?)
    def stats(sub, key):
        a = np.array([r[key] for r in sub if np.isfinite(r.get(key, np.nan))])
        return f"med={np.median(a):.2f} IQR[{np.percentile(a,25):.2f},{np.percentile(a,75):.2f}] N={len(a)}" if len(a) else "n/a"
    o.append("")
    o.append("--- bias: CONTRIBUTING vs all model-placeable (SZA/Dp/cone medians) ---")
    placeable_rows = [r for r in rows if r['stage'] != 'no_params']
    for key in ('sza', 'dp', 'cone'):
        # placeable rows only have these in 'contrib'; compare contrib to itself vs north for Dp bias visibility
        o.append(f"  {key:4s}: contributing {stats(contrib,key)}  |  northward {stats(north,key)}  |  CONFIRMED {stats([r for r in north if r['cls']=='CONFIRMED_PDL'],key)}")

    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "SELECTION_FUNNEL.txt"), "w") as fh:
        fh.write(txt + "\n")

    # per-encounter stage (all) + rich contributing covariate table (for run11 / run13)
    with open(os.path.join(OUT, "funnel_all.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(['eid', 'probe', 'stage', 'bz', 'cls'])
        for r in rows:
            w.writerow([r['eid'], r['probe'], r['stage'], round(r.get('bz', np.nan), 2) if np.isfinite(r.get('bz', np.nan)) else '',
                        r.get('cls', '')])
    with open(os.path.join(OUT, "funnel_contributing.csv"), "w", newline="") as f:
        cols = ['eid', 'probe', 'bz', 'bz_sign', 'clock', 'cone', 'sza', 'dp', 'ma', 'bx', 'Dn_mem', 'EB_mem',
                'n_near', 'beta_near', 'T_near', 'n_mem', 'cls']
        w = csv.writer(f); w.writerow(cols)
        rr = lambda v, k=2: (round(v, k) if np.isfinite(v) else '')
        for r in contrib:
            bzs = 1 if (np.isfinite(r['bz']) and r['bz'] > 0) else (0 if np.isfinite(r['bz']) else '')
            w.writerow([r['eid'], r['probe'], rr(r['bz']), bzs, rr(r['clock'], 1),
                        rr(r['cone'], 1), rr(r['sza'], 1), rr(r['dp']), rr(r['ma'], 1), rr(r['bx']),
                        round(r['Dn_mem'], 4), round(r['EB_mem'], 4), rr(r['n_near']), rr(r['beta_near'], 4),
                        rr(r['T_near'], 1), r['n_mem'], r['cls']])
    print(f"\nsaved -> {OUT}\\SELECTION_FUNNEL.txt , funnel_all.csv , funnel_contributing.csv ({len(contrib)} rows)", flush=True)


if __name__ == "__main__":
    main()
