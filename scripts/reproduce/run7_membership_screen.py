"""RUN 7 — MAGNETOSHEATH-MEMBERSHIP SCREEN (temperature + velocity + density + position).

The near-MP density depletion lives in low-beta plasma, and beta can't tell a real sheath
PDL from magnetosphere/boundary-layer leakage (both low-beta). Temperature and flow CAN:
  - magnetosheath (incl. a PDL): SHOCKED solar wind -> T ~ 100s eV (~ background sheath T), still flowing.
  - magnetosphere: hot plasma-sheet (keV) or near-vacuum, flow collapsed.

Per sample: T[eV] = p_th[nPa] / (n[cm^-3] * 1.602e-4);  v = vmag[km/s].
Per encounter, characterise the confirmed-sheath BACKGROUND (s in [0.6,1.0]): T_bg, v_bg.
Then for each shell report (a) DIAGNOSTIC: how sheath-like the near-MP plasma is (T/T_bg, v/v_bg),
and (b) the depletion among MEMBERSHIP-CONFIRMED sheath samples only.

Membership (sample is magnetosheath) = geometry-sheath AND n>N_FLOOR AND T in [T_bg/TBAND, T_bg*TBAND]
AND v > VFRAC*v_bg.  Threaded load via psub.
(NOTE: energy spectra are not in the substrate; T is the moment proxy. Spectra need a re-fetch -> Run 8 atlas.)
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run7_membership")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0)
EDGES = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00]
MIN_SAMP_BIN = 3; MIN_BG = 5; MIN_NENC = 30; NB = len(EDGES) - 1
KB = 1.602e-4          # nPa per (cm^-3 * eV)
N_FLOOR = 0.3          # cm^-3
TBAND = 3.0            # T within [T_bg/3, T_bg*3]
VFRAC = 0.2            # v > 0.2 * v_bg


def one(d):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
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
    if len(sv) < MIN_BG + MIN_SAMP_BIN:
        return None
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return None
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm])
    T_bg = np.nanmedian(Tv[bgm]); v_bg = np.nanmedian(vv[bgm])
    if not (n_bg > 0 and b_bg > 0 and np.isfinite(T_bg) and T_bg > 0):
        return None
    member = (nv > N_FLOOR) & np.isfinite(Tv) & (Tv > T_bg / TBAND) & (Tv < T_bg * TBAND)
    if np.isfinite(v_bg) and v_bg > 0:
        member = member & np.isfinite(vv) & (vv > VFRAC * v_bg)
    out = {}
    for i in range(NB):
        lo, hi = EDGES[i], EDGES[i + 1]
        m = (sv >= lo) & (sv <= hi) if i == NB - 1 else (sv >= lo) & (sv < hi)
        if m.sum() < MIN_SAMP_BIN:
            continue
        lowb = m & (Tv < T_bg)  # placeholder not used; keep simple
        rec = dict(Dn_all=float(np.median(nv[m]) / n_bg),
                   Tratio=float(np.nanmedian(Tv[m]) / T_bg) if np.isfinite(np.nanmedian(Tv[m])) else np.nan,
                   vratio=(float(np.nanmedian(vv[m]) / v_bg) if (np.isfinite(v_bg) and v_bg > 0 and np.isfinite(np.nanmedian(vv[m]))) else np.nan),
                   memfrac=float(np.mean(member[m])))
        mm = m & member
        if mm.sum() >= MIN_SAMP_BIN:
            rec['Dn_mem'] = float(np.median(nv[mm]) / n_bg)
            rec['EB_mem'] = float(np.median(bv[mm]) / b_bg)
        out[i] = rec
    return out


def main():
    res = pmap(one)
    print(f"encounters contributing: {len(res)}", flush=True)
    agg = {i: {k: [] for k in ('Dn_all', 'Tratio', 'vratio', 'memfrac', 'Dn_mem', 'EB_mem')} for i in range(NB)}
    for r in res:
        for i, rec in r.items():
            for k, val in rec.items():
                if k in agg[i] and val is not None and np.isfinite(val):
                    agg[i][k].append(val)

    def med(i, k):
        v = agg[i][k]
        return np.median(v) if len(v) >= MIN_NENC else np.nan

    lines = ["RUN 7 — magnetosheath-membership screen (T + v + n + position)",
             f"membership: geometry-sheath & n>{N_FLOOR} & T in [T_bg/{TBAND:.0f},T_bg*{TBAND:.0f}] & v>{VFRAC}*v_bg",
             "",
             f"{'shell':13s} {'N':>5s} {'Dn_all':>7s} {'T/Tbg':>7s} {'v/vbg':>7s} {'mem%':>6s} {'Dn_mem':>7s} {'EB_mem':>7s}"]
    for i in range(NB):
        n_enc = len(agg[i]['Dn_all'])
        if n_enc < MIN_NENC:
            lines.append(f"[{EDGES[i]:.2f},{EDGES[i+1]:.2f}]  {n_enc:5d}   (N<{MIN_NENC})"); continue
        lines.append(f"[{EDGES[i]:.2f},{EDGES[i+1]:.2f}]  {n_enc:5d} {med(i,'Dn_all'):7.3f} {med(i,'Tratio'):7.2f} "
                     f"{med(i,'vratio'):7.2f} {med(i,'memfrac'):6.0%} {med(i,'Dn_mem'):7.3f} {med(i,'EB_mem'):7.3f}")
    txt = "\n".join(lines)
    print(txt, flush=True)
    with open(os.path.join(OUT, "RUN7_membership.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
