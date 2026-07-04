"""RUN 3 — BIN MULTIVERSE (guard against bin / background / contamination cherry-picking).

For the corrected geometry, sweep:
  - geometry handling : 'sheath' (drop magnetosphere by geometry) vs 'sheath_b' (also drop beta<0.2)
  - background shell   : [0.5,0.8], [0.6,1.0], [0.7,1.0]
  - near shell         : [0.10,0.20], [0.10,0.15], [0.15,0.20]
Report Dn and EB (median over encounters, N>=30) across the whole grid. If the depletion
(Dn<1) and pile-up (EB>1) hold across the grid, they aren't an artifact of one bin choice.
Threaded load via psub.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run3_bin_multiverse")
os.makedirs(OUT, exist_ok=True)
GEOMS = ['sheath', 'sheath_b']
BGS = [(0.5, 0.8), (0.6, 1.0), (0.7, 1.0)]
NEARS = [(0.10, 0.20), (0.10, 0.15), (0.15, 0.20)]
BETA_MAG = 0.2; MIN_BG = 5; MIN_SAMP = 3; MIN_NENC = 30


def multi_one(d):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float); beta = d['beta'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    base = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    out = {}
    for geom in GEOMS:
        mask = base & (beta >= BETA_MAG) if geom == 'sheath_b' else base
        sv = s[mask]; nv = n[mask]; bv = b[mask]
        if len(sv) < MIN_BG + MIN_SAMP:
            continue
        for bgi, (blo, bhi) in enumerate(BGS):
            bgm = (sv >= blo) & (sv <= bhi)
            if bgm.sum() < MIN_BG:
                continue
            n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm])
            if not (n_bg > 0 and b_bg > 0):
                continue
            for ni, (nlo, nhi) in enumerate(NEARS):
                m = (sv >= nlo) & (sv < nhi)
                if m.sum() < MIN_SAMP:
                    continue
                out[(geom, bgi, ni)] = (float(np.median(nv[m]) / n_bg), float(np.median(bv[m]) / b_bg))
    return out or None


def main():
    results = pmap(multi_one)
    print(f"encounters contributing: {len(results)}", flush=True)
    agg = {}
    for r in results:
        for k, (dn, eb) in r.items():
            agg.setdefault(k, ([], []))
            agg[k][0].append(dn); agg[k][1].append(eb)

    def cell(geom, bgi, ni, which):
        v = agg.get((geom, bgi, ni))
        if not v or len(v[0]) < MIN_NENC:
            return "   (N<30)"
        arr = v[0] if which == 'Dn' else v[1]
        return f"{np.median(arr):.3f}(N{len(arr)})"

    lines = []
    for which in ('Dn', 'EB'):
        lines.append(f"\n================  {which}  ================")
        for geom in GEOMS:
            lines.append(f"\n--- geometry={geom} ---")
            lines.append(f"{'near/bg':14s}" + "".join(f"{str(bg):>16s}" for bg in BGS))
            for ni, near in enumerate(NEARS):
                row = f"{str(near):14s}" + "".join(f"{cell(geom,bgi,ni,which):>16s}" for bgi in range(len(BGS)))
                lines.append(row)
    txt = "\n".join(lines)
    print(txt, flush=True)
    with open(os.path.join(OUT, "RUN3_bin_multiverse.txt"), "w") as f:
        f.write("RUN 3 — bin / background / contamination multiverse (sheath_geom; Dn=n/n_bg, EB=|B|/|B|_bg)\n")
        f.write("sheath = drop magnetosphere by geometry; sheath_b = also drop beta<0.2\n")
        f.write(txt + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
