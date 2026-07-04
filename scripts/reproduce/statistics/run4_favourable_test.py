"""RUN 4 — FAVOURABLE-REGIME PDL TEST (pre-registered).

Per user pre-registration (set before seeing results):
  MAIN cut  : Bz>0 (northward) AND cone>45deg AND SZA<20deg AND 1<Dp<4 nPa
  STRICT    : Bz>2 AND cone>60deg AND SZA<15deg            (sensitivity)
  ALL       : no cut                                       (baseline)
Geometry = Jelinek radial; samples = MEMBERSHIP-CONFIRMED magnetosheath (T+v+n+position, Run 7).
Main PDL shell = s in [0.05,0.20); the s<0.05 layer is reported SEPARATELY (contamination).
A descriptive Bz x cone x SZA grid is also reported.

LANGUAGE DISCIPLINE: positives are a "membership-verified PDL-relevant favourable-regime
profile", NOT a "confirmed PDL occurrence" (that requires Run 8 spectra/atlas).
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

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run4_favourable")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0)
SHELLS = [('contam[0.00,0.05)', 0.00, 0.05), ('PDL[0.05,0.20)', 0.05, 0.20),
          ('[0.05,0.10)', 0.05, 0.10), ('[0.10,0.15)', 0.10, 0.15), ('[0.15,0.20)', 0.15, 0.20),
          ('[0.20,0.40)', 0.20, 0.40)]
MIN_SAMP = 3; MIN_BG = 5; MIN_NENC = 30
KB = 1.602e-4; N_FLOOR = 0.3; TBAND = 3.0; VFRAC = 0.2


def one(d):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    bz = float(d['bz']) if np.isfinite(d['bz']) else np.nan
    cone = float(d['cone_deg']) if np.isfinite(d['cone_deg']) else np.nan
    sza = float(d['sza_deg']) if np.isfinite(d['sza_deg']) else np.nan
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    geo = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]
    if len(sv) < MIN_BG + MIN_SAMP:
        return None
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return None
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm]); T_bg = np.nanmedian(Tv[bgm]); v_bg = np.nanmedian(vv[bgm])
    if not (n_bg > 0 and b_bg > 0 and np.isfinite(T_bg) and T_bg > 0):
        return None
    member = (nv > N_FLOOR) & np.isfinite(Tv) & (Tv > T_bg / TBAND) & (Tv < T_bg * TBAND)
    if np.isfinite(v_bg) and v_bg > 0:
        member = member & np.isfinite(vv) & (vv > VFRAC * v_bg)
    rec = {}
    for name, lo, hi in SHELLS:
        m = (sv >= lo) & (sv < hi)
        if m.sum() < MIN_SAMP:
            continue
        e = dict(Dn_all=float(np.median(nv[m]) / n_bg), memfrac=float(np.mean(member[m])))
        mm = m & member
        if mm.sum() >= MIN_SAMP:
            e['Dn_mem'] = float(np.median(nv[mm]) / n_bg); e['EB_mem'] = float(np.median(bv[mm]) / b_bg)
        rec[name] = e
    return (rec, bz, cone, sza, dp)


def passes(tag, bz, cone, sza, dp):
    if tag == 'ALL':
        return True
    ok = np.isfinite(bz) and np.isfinite(cone) and np.isfinite(sza) and np.isfinite(dp)
    if not ok:
        return False
    if tag == 'MAIN':
        return bz > 0 and cone > 45 and sza < 20 and 1.0 < dp < 4.0
    if tag == 'STRICT':
        return bz > 2 and cone > 60 and sza < 15
    return False


def agg_table(res, tag):
    sub = [r for r in res if passes(tag, r[1], r[2], r[3], r[4])]
    lines = [f"--- cut={tag}  (encounters={len(sub)}) ---",
             f"{'shell':18s} {'N':>5s} {'Dn_all':>7s} {'mem%':>6s} {'Dn_mem':>7s} {'EB_mem':>7s}"]
    for name, lo, hi in SHELLS:
        dnm = [r[0][name]['Dn_mem'] for r in sub if name in r[0] and 'Dn_mem' in r[0][name]]
        ebm = [r[0][name]['EB_mem'] for r in sub if name in r[0] and 'EB_mem' in r[0][name]]
        dna = [r[0][name]['Dn_all'] for r in sub if name in r[0]]
        mfr = [r[0][name]['memfrac'] for r in sub if name in r[0]]
        n = len(dnm)
        flag = "" if n >= MIN_NENC else " (N<30)"
        if not dna:
            lines.append(f"{name:18s}   (no data)"); continue
        dn_mem = f"{np.median(dnm):.3f}" if dnm else "  -  "
        eb_mem = f"{np.median(ebm):.3f}" if ebm else "  -  "
        lines.append(f"{name:18s} {n:5d} {np.median(dna):7.3f} {np.median(mfr):6.0%} {dn_mem:>7s} {eb_mem:>7s}{flag}")
    return "\n".join(lines)


def main():
    res = pmap(one)
    print(f"encounters contributing: {len(res)}", flush=True)
    out = ["RUN 4 — membership-verified PDL-relevant FAVOURABLE-REGIME profile (NOT a confirmed PDL occurrence)",
           "geometry=Jelinek radial; samples=membership-confirmed sheath; main PDL shell s in [0.05,0.20); s<0.05 = contamination layer", ""]
    for tag in ('ALL', 'MAIN', 'STRICT'):
        out.append(agg_table(res, tag)); out.append("")
    # descriptive grid: Dn_mem & EB_mem in PDL[0.05,0.20) by Bz-sign x cone x SZA
    out.append("--- DESCRIPTIVE GRID: PDL[0.05,0.20) Dn_mem / EB_mem  (cells with N>=30 only) ---")
    out.append(f"{'Bz':>6s} {'cone':>8s} {'SZA':>8s} {'N':>5s} {'Dn_mem':>7s} {'EB_mem':>7s}")
    for bzlab, bzf in [('Bz<0', lambda v: v < 0), ('Bz>0', lambda v: v > 0)]:
        for clab, cf in [('cone<45', lambda v: v < 45), ('cone>45', lambda v: v >= 45)]:
            for slab, sf in [('SZA<15', lambda v: v < 15), ('SZA15-30', lambda v: v >= 15)]:
                cell = [r for r in res if np.isfinite(r[1]) and np.isfinite(r[2]) and np.isfinite(r[3])
                        and bzf(r[1]) and cf(r[2]) and sf(r[3])]
                dnm = [r[0]['PDL[0.05,0.20)']['Dn_mem'] for r in cell
                       if 'PDL[0.05,0.20)' in r[0] and 'Dn_mem' in r[0]['PDL[0.05,0.20)']]
                ebm = [r[0]['PDL[0.05,0.20)']['EB_mem'] for r in cell
                       if 'PDL[0.05,0.20)' in r[0] and 'EB_mem' in r[0]['PDL[0.05,0.20)']]
                if len(dnm) >= MIN_NENC:
                    out.append(f"{bzlab:>6s} {clab:>8s} {slab:>8s} {len(dnm):5d} {np.median(dnm):7.3f} {np.median(ebm):7.3f}")
    txt = "\n".join(out)
    print(txt, flush=True)
    with open(os.path.join(OUT, "RUN4_favourable.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
