"""SELECTION FUNCTION: which encounters reach which s-shells, and is the
near-MP-sampling subset biased in SZA / Dp / cone / probe vs the full eligible set?
Threaded over the substrate via psub."""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, collections
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run2_profile_cube")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0); MIN_BG = 5; MIN_SAMP = 3


def s_sheath(d):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    valid = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    return s[valid]


def fnum(d, k):
    try:
        v = float(d[k]); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def sel_one(name, d):
    sv = s_sheath(d)
    probe = name[:-4].rsplit('_', 1)[-1]
    if sv is None or len(sv) == 0:
        return dict(elig=0, probe=probe, sza=np.nan, dp=np.nan, cone=np.nan)
    has_bg = int(((sv >= BG[0]) & (sv <= BG[1])).sum() >= MIN_BG)
    rr = lambda lo, hi: int(((sv >= lo) & (sv < hi)).sum() >= MIN_SAMP)
    clean = rr(0.10, 0.15) or rr(0.15, 0.20)
    return dict(elig=1, has_bg=has_bg, r020=rr(0.0, 0.2), r0010=rr(0.0, 0.10),
                clean=int(clean), contrib=int(has_bg and clean),
                sza=fnum(d, 'sza_deg'), dp=fnum(d, 'dp'), cone=fnum(d, 'cone_deg'), probe=probe)


def main():
    rows = pmap(sel_one, with_name=True)
    elig = [r for r in rows if r.get('elig')]
    withbg = [r for r in elig if r['has_bg']]
    contrib = [r for r in elig if r.get('contrib')]
    o = []
    o.append(f"total substrate npz        : {len(rows)}")
    o.append(f"with sheath samples        : {len(elig)}")
    o.append(f"with background (s[0.6,1]) : {len(withbg)}")
    o.append(f"reach a clean shell [0.10,0.20] : {sum(r['clean'] for r in elig)}")
    o.append(f"CONTRIBUTING (bg AND clean shell): {len(contrib)}")
    o.append(f"  of with-bg, reach s<0.20 : {sum(r['r020'] for r in withbg)}/{len(withbg)}")
    o.append(f"  of with-bg, reach s<0.10 : {sum(r['r0010'] for r in withbg)}/{len(withbg)}")
    o.append("contributing by probe      : " + str(dict(collections.Counter(r['probe'] for r in contrib))))

    def stats(sub, key):
        a = np.array([r[key] for r in sub if np.isfinite(r.get(key, np.nan))])
        return f"med={np.median(a):.2f} IQR[{np.percentile(a,25):.2f},{np.percentile(a,75):.2f}] N={len(a)}" if len(a) else "n/a"

    o.append("--- bias check: contributing vs all-eligible ---")
    for key in ('sza', 'dp', 'cone'):
        o.append(f"  {key:4s}: eligible {stats(elig,key)}  |  contributing {stats(contrib,key)}")
    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "SELECTION_FUNCTION.txt"), "w") as fh:
        fh.write("SELECTION FUNCTION\n\n" + txt + "\n")
    print("\nsaved -> SELECTION_FUNCTION.txt", flush=True)


if __name__ == "__main__":
    main()
