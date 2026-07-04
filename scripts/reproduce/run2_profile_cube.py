"""RUN 2 — RADIAL_JELINEK_PROFILE_CUBE (the centerpiece).

Builds the radial profile of magnetosheath plasma vs normalised position s
(Shue MP + Jelinek-2012 BS geometry) from the frozen substrate, to answer:
does a real near-magnetopause depletion/pileup profile survive under correct
geometry — and is the deep `s<0.1` drop a real PDL or magnetosphere contamination?

Two geometries are compared per constraint-A transparency:
  CLIP        : s clipped to [0,1] (OLD behaviour) -> magnetospheric samples
                (r < model MP) pile into the near bin = the contamination artefact.
  SHEATH_GEOM : s UNCLIPPED; keep only samples strictly between model MP and BS
                (d_mp>0 & d_bs>0) -> magnetosphere dropped, not dumped at s=0.

Per s-shell: Dn (n/n_bg), EB (|B|/|B|_bg), median beta, and magnetosphere-fraction
(frac samples with beta<BETA_MAG) as a contamination flag. Background = s[0.6,1.0].
Aggregate shell reported only if N_enc>=30 (constraint B). Single pass over the substrate.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, glob, csv
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap, list_files

SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run2_profile_cube")
os.makedirs(OUT, exist_ok=True)

BG = (0.6, 1.0)
EDGES = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00]
MIN_SAMP_BIN = 3
MIN_BG = 5
MIN_NENC = 30
BETA_MAG = 0.2


def s_unclipped(x, y, z, mp0, alpha, dp):
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp
    d_bs = r_bs - r
    denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    return s, d_mp, d_bs


def encounter_profile(d, mode):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float); beta = d['beta'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    s, d_mp, d_bs = s_unclipped(x, y, z, mp0, alpha, dp)
    base = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0)
    if mode == 'clip':
        s = np.where(np.isfinite(s), np.clip(s, 0.0, 1.0), np.nan)
        valid = base & np.isfinite(s)
    else:
        valid = base & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[valid]; nv = n[valid]; bv = b[valid]; betav = beta[valid]
    if len(sv) < (MIN_BG + MIN_SAMP_BIN):
        return None
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return None
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm])
    if not (n_bg > 0 and b_bg > 0):
        return None
    out = {}
    nb = len(EDGES) - 1
    for i in range(nb):
        lo, hi = EDGES[i], EDGES[i + 1]
        m = (sv >= lo) & (sv <= hi) if i == nb - 1 else (sv >= lo) & (sv < hi)
        if m.sum() < MIN_SAMP_BIN:
            continue
        bm = betav[m]
        out[i] = dict(Dn=float(np.median(nv[m]) / n_bg), EB=float(np.median(bv[m]) / b_bg),
                      beta=float(np.nanmedian(bm)) if np.isfinite(bm).any() else np.nan,
                      magfrac=float(np.mean(bm < BETA_MAG)) if np.isfinite(bm).any() else np.nan)
    return out


def favourable(d):
    try:
        bz = float(d['bz']); cone = float(d['cone_deg']); sza = float(d['sza_deg']); dp = float(d['dp'])
    except Exception:
        return False
    return (np.isfinite(bz) and bz > 0 and np.isfinite(cone) and cone > 45
            and np.isfinite(sza) and sza < 20 and np.isfinite(dp) and 1.0 <= dp <= 4.0)


def rows_from(perbin):
    nb = len(EDGES) - 1
    rows = []
    for i in range(nb):
        vals = perbin[i]; nenc = len(vals)
        if nenc < MIN_NENC:
            rows.append([EDGES[i], EDGES[i + 1], nenc, None, None, None, None]); continue
        bvals = [v['beta'] for v in vals if np.isfinite(v['beta'])]
        mvals = [v['magfrac'] for v in vals if np.isfinite(v['magfrac'])]
        rows.append([EDGES[i], EDGES[i + 1], nenc,
                     float(np.median([v['Dn'] for v in vals])),
                     float(np.median([v['EB'] for v in vals])),
                     float(np.median(bvals)) if bvals else float('nan'),
                     float(np.median(mvals)) if mvals else float('nan')])
    return rows


def fmt(rows):
    out = [f"  {'s-shell':13s} {'N_enc':>6s} {'Dn':>7s} {'EB':>7s} {'beta':>7s} {'mag_frac':>8s}"]
    for lo, hi, nenc, Dn, EB, beta, magf in rows:
        if Dn is None:
            out.append(f"  [{lo:.2f},{hi:.2f}]  {nenc:6d}   (N<{MIN_NENC})")
        else:
            out.append(f"  [{lo:.2f},{hi:.2f}]  {nenc:6d} {Dn:7.3f} {EB:7.3f} {beta:7.2f} {magf:8.1%}")
    return "\n".join(out)


def save_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['s_lo', 's_hi', 'N_enc', 'Dn', 'EB', 'beta_med', 'mag_frac'])
        for r in rows:
            w.writerow([r[0], r[1], r[2]] + ['' if x is None else round(x, 4) for x in r[3:]])


def profile_one(d):
    pc = encounter_profile(d, 'clip')
    ps = encounter_profile(d, 'sheath_geom')
    if pc is None and ps is None:
        return None
    return (pc, ps, favourable(d))


def main():
    files = list_files()
    print(f"substrate encounters: {len(files)}", flush=True)
    nb = len(EDGES) - 1
    acc = {k: {i: [] for i in range(nb)} for k in ('clip', 'sheath', 'fav')}
    used = {'clip': 0, 'sheath': 0, 'fav': 0}
    for (pc, ps, isfav) in pmap(profile_one, files=files):   # threaded over the substrate (GIL-released np.load + numpy)
        if pc:
            used['clip'] += 1
            for i, v in pc.items(): acc['clip'][i].append(v)
        if ps:
            used['sheath'] += 1
            if isfav: used['fav'] += 1
            for i, v in ps.items():
                acc['sheath'][i].append(v)
                if isfav: acc['fav'][i].append(v)
    report = []
    titles = {'clip': "OVERALL — geometry=CLIP (old; shows magnetosphere contamination)",
              'sheath': "OVERALL — geometry=SHEATH_GEOM (corrected: magnetosphere dropped)",
              'fav': "FAVOURABLE conds — SHEATH_GEOM (Bz>0, cone>45, sza<20, 1<=Dp<=4)"}
    for k in ('clip', 'sheath', 'fav'):
        rows = rows_from(acc[k])
        save_csv(os.path.join(OUT, f"profile_{k}.csv"), rows)
        hdr = f"\n=== {titles[k]}  (N_enc used={used[k]}) ==="
        print(hdr, flush=True); print(fmt(rows), flush=True); report.append(hdr + "\n" + fmt(rows))
    with open(os.path.join(OUT, "RUN2_SUMMARY.txt"), 'w') as f:
        f.write("RUN 2 — RADIAL_JELINEK_PROFILE_CUBE\n")
        f.write("Dn=n/n_bg, EB=|B|/|B|_bg, bg=s[0.6,1.0]; mag_frac=frac samples beta<%.2f (contamination flag)\n\n" % BETA_MAG)
        f.write("\n".join(report) + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
