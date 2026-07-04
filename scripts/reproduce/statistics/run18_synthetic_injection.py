"""RUN 18 — synthetic-injection forward test of the dynamic-pressure ordering.

The sceptical review objection: the observed Dn-vs-Dp ordering may be the coordinate's own Dp dependence
(s is built from Dp-parameterised Shue MP + Jelinek BS), not a physical depth change. This forward test
discriminates artefact from physics: assign SYNTHETIC plasma (which is Dp-INDEPENDENT in PHYSICAL terms)
to every REAL sample using REAL geometry / REAL OMNI Dp / REAL Shue-Jelinek boundaries, then push it
through the IDENTICAL s-shell pipeline and see whether a Dp ordering appears spuriously.

Scenarios (per user spec):
  null            : n=1, |B|=1 everywhere (uniform; no depletion)            -> expect Dn=EB=1 flat
  abs-700/1000/1300: depletion of FIXED absolute thickness L km just outside the model MP, FIXED amplitude
                     (n=DEPL, |B|=ENH within 0<dist<L km; else 1), NOT Dp-dependent
  norm-layer      : depletion FIXED in s in [0.10,0.20), FIXED amplitude, NOT Dp-dependent (control)
Each measured EXACTLY as the real pipeline: Dn = median(n_syn in s[0.05,0.20)) / median(n_syn in s[0.6,1.0]);
EB likewise; aggregated over contributing encounters and by the OBSERVED Dp tertiles (edges 1.91, 3.14).
Artefact fraction = synthetic (Dn_high - Dn_low) / observed (1.158 - 0.711). Threaded load via psub.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, collections
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run18_synthetic")
os.makedirs(OUT, exist_ok=True)
RE_KM = 6371.0
BG = (0.6, 1.0); PDL = (0.05, 0.20); MIN_BG = 5; MIN_SHELL = 5
DEPL = 0.5   # fixed density depletion amplitude (n -> 0.5 of background) inside the layer
ENH = 2.0    # fixed field enhancement (|B| -> 2x background) inside the layer
DP_EDGES = (1.91, 3.14)      # observed tertile edges (run16), for like-for-like binning
OBS_DN = (0.711, 0.928, 1.158)   # observed Dn by Dp tertile (run16)
OBS_EB = (2.257, 2.077, 1.633)
ABS_L = {'abs-700': 700.0, 'abs-1000': 1000.0, 'abs-1300': 1300.0}


SCENARIOS = ['null', 'abs-700', 'abs-1000', 'abs-1300', 'norm-layer', 'dpcorr+', 'dpcorr-']


def layer_mask(scen, s, dist_km, s_true_plus, s_true_minus):
    """Boolean 'in synthetic depletion layer' for each scenario (measured-s for fixed scenarios;
    TRUE-s for the Dp-correlated-boundary-error scenarios)."""
    if scen == 'null':
        return np.zeros_like(s, dtype=bool)
    if scen == 'norm-layer':
        return (s >= 0.10) & (s < 0.20)
    if scen == 'dpcorr+':
        return (s_true_plus >= 0.10) & (s_true_plus < 0.20)
    if scen == 'dpcorr-':
        return (s_true_minus >= 0.10) & (s_true_minus < 0.20)
    L = ABS_L[scen]
    return (dist_km > 0) & (dist_km < L)


def one(name, d):
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    geo = (d_mp > 0) & (d_bs > 0) & np.isfinite(denom) & (denom > 0)
    if geo.sum() < MIN_BG + MIN_SHELL:
        return None
    dmp_g = d_mp[geo]; den_g = denom[geo]
    s = dmp_g / den_g
    dist_km = dmp_g * RE_KM
    # Dp-correlated magnetopause mis-placement (worst case ~ +/-0.5 R_E across the Dp range), measured with
    # the model MP but layer defined in the TRUE (offset) coordinate -> tests a Dp-correlated boundary error.
    off = 0.5 * float(np.clip((dp - 2.4) / 2.0, -1.2, 1.2))     # R_E, increases with Dp
    s_true_plus = (dmp_g - off) / np.where((den_g - off) > 0, den_g - off, np.nan)
    s_true_minus = (dmp_g + off) / np.where((den_g + off) > 0, den_g + off, np.nan)
    bgm = (s >= BG[0]) & (s <= BG[1]); shm = (s >= PDL[0]) & (s < PDL[1])
    if bgm.sum() < MIN_BG or shm.sum() < MIN_SHELL:
        return None
    rec = {'dp': dp, 'scen': {}}
    for scen in SCENARIOS:
        m = layer_mask(scen, s, dist_km, s_true_plus, s_true_minus)
        n_syn = np.where(m, DEPL, 1.0); b_syn = np.where(m, ENH, 1.0)
        n_bg = np.median(n_syn[bgm]); b_bg = np.median(b_syn[bgm])
        if n_bg <= 0 or b_bg <= 0:
            continue
        rec['scen'][scen] = (float(np.median(n_syn[shm]) / n_bg), float(np.median(b_syn[shm]) / b_bg))
    return rec


def tert(rows, scen, idx):
    """median Dn(idx=0) or EB(idx=1) per Dp tertile for a scenario."""
    sub = [(r['dp'], r['scen'][scen][idx]) for r in rows if scen in r['scen'] and np.isfinite(r['dp'])]
    lo = [v for dp, v in sub if dp <= DP_EDGES[0]]
    mid = [v for dp, v in sub if DP_EDGES[0] < dp <= DP_EDGES[1]]
    hi = [v for dp, v in sub if dp > DP_EDGES[1]]
    m = lambda a: float(np.median(a)) if a else np.nan
    return m(lo), m(mid), m(hi), len(sub)


def main():
    rows = pmap(one, with_name=True)
    o = ["RUN 18 — synthetic-injection forward test",
         f"synthetic depletion amplitude n={DEPL}, field enhancement |B|={ENH}; Dp tertile edges {DP_EDGES}",
         f"contributing synthetic encounters: {len(rows)}",
         "",
         f"OBSERVED (real data):           Dn  {OBS_DN[0]:.3f} / {OBS_DN[1]:.3f} / {OBS_DN[2]:.3f}   "
         f"EB  {OBS_EB[0]:.3f} / {OBS_EB[1]:.3f} / {OBS_EB[2]:.3f}   (slope Dn = {OBS_DN[2]-OBS_DN[0]:+.3f})",
         ""]
    obs_slope = OBS_DN[2] - OBS_DN[0]
    o.append(f"{'scenario':12s} {'N':>5s}  {'Dn low/mid/high':>22s} {'Dn slope':>9s} {'art.frac':>9s}   {'EB low/mid/high':>22s}")
    for scen in ['null', 'norm-layer', 'abs-700', 'abs-1000', 'abs-1300', 'dpcorr+', 'dpcorr-']:
        dl, dm, dh, n = tert(rows, scen, 0)
        el, em, eh, _ = tert(rows, scen, 1)
        slope = dh - dl
        af = slope / obs_slope if np.isfinite(slope) and obs_slope != 0 else np.nan
        o.append(f"{scen:12s} {n:5d}  {f'{dl:.3f}/{dm:.3f}/{dh:.3f}':>22s} {slope:+9.3f} {af:>9.2f}   {f'{el:.3f}/{em:.3f}/{eh:.3f}':>22s}")
    o.append("")
    o.append("READING: art.frac = synthetic Dn slope / observed Dn slope (+0.447).")
    o.append("  ~1  => the observed Dp ordering is reproduced by a Dp-INDEPENDENT layer => coordinate/sampling artefact.")
    o.append("  ~0 or NEGATIVE => a Dp-independent layer does NOT produce the observed (rising) trend => observed ordering is")
    o.append("        physical-leaning (the artefact channel has small or opposite-sign effect); bound = |art.frac|.")
    o.append("  norm-layer is the control: a layer fixed in s must be Dp-flat (art.frac ~ 0).")
    o.append("  dpcorr+/- inject a WORST-CASE Dp-CORRELATED magnetopause error (~+/-0.5 R_E across the Dp range) on an")
    o.append("        otherwise Dp-flat layer: if these stay ~flat, even a Dp-correlated boundary error cannot manufacture")
    o.append("        the observed ordering (the channel the Devil's Advocate named is excluded).")
    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "SYNTHETIC_INJECTION.txt"), "w") as f:
        f.write(txt + "\n")
    # csv of per-encounter for audit
    with open(os.path.join(OUT, "synthetic_per_encounter.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(['dp'] + [f'{s}_Dn' for s in SCENARIOS])
        for r in rows:
            w.writerow([round(r['dp'], 3)] + [round(r['scen'][s][0], 4) if s in r['scen'] else '' for s in SCENARIOS])
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
