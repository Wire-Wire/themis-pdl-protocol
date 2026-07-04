"""run31a + run31c (no network).

31a - representativeness of the 60 spectrally supported candidates: probe / year / geometry
spread, concentration (Herfindahl), and a leave-one-probe-out / leave-one-year-out robustness
check on the two headline A-vs-NN axes (cone angle, background density) to show neither finding
hinges on a single probe or year.

31c - theta_Bn vs cone-angle collinearity, analytic: the Jelinek bow-shock normal tilt from the
Sun-Earth line as a function of polar angle. For a near-subsolar archive this bounds |theta_Bn -
cone|, and because A and NN share the same SZA distribution the tilt is a near-common offset, so
the A-vs-NN cone difference maps into theta_Bn essentially unchanged. Retires the "global cone vs
local shock geometry" question without a streamline assumption.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import sys, csv
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
R28 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions\run28_per_encounter.csv")
with open(R28, newline="") as f:
    recs = list(csv.DictReader(f))


def col(rows, c, num=True):
    out = []
    for r in rows:
        v = r[c]
        out.append(float(v) if num and v not in ("", "nan") else v)
    return np.array(out)


A = [r for r in recs if r["grp"] == "A"]
NN = [r for r in recs if r["grp"] == "NONCAND" and r["north"] == "True"]
for r in recs:
    r["probe"] = r["eid"].split("_")[1]
    r["year"] = int(r["eid"][:4])

L = []
def log(s=""):
    print(s, flush=True); L.append(s)

log("RUN31a - REPRESENTATIVENESS OF THE 60 SPECTRALLY SUPPORTED CANDIDATES")
probes = [r["eid"].split("_")[1] for r in A]
years = [int(r["eid"][:4]) for r in A]
up, cp = np.unique(probes, return_counts=True)
uy, cy = np.unique(years, return_counts=True)
log(f"  N=60; probes: " + ", ".join(f"{p}:{c}" for p, c in zip(up, cp)))
log(f"    top-probe share {100*cp.max()/60:.0f}%  Herfindahl {np.sum((cp/60)**2):.2f} (even={1/len(up):.2f})")
log(f"    note: THB/THC left for lunar ARTEMIS orbit in 2010, so the Earth-dayside load is carried by THA/THD/THE")
log(f"  years: {len(uy)} distinct, {uy.min()}-{uy.max()}; top-year share {100*cy.max()/60:.0f}%; "
    f"top-3 {100*np.sort(cy)[::-1][:3].sum()/60:.0f}%  Herfindahl {np.sum((cy/60)**2):.2f} (even={1/len(uy):.2f})")
for c in ("sza", "cone", "mp0"):
    a = col(A, c)
    log(f"  {c:4s}: median {np.median(a):.1f}  IQR [{np.percentile(a,25):.1f},{np.percentile(a,75):.1f}]  "
        f"range [{a.min():.1f},{a.max():.1f}]")
log(f"  near-subsolar: {100*np.mean(col(A,'sza')<30):.0f}% have SZA<30 deg")
log(f"  SZA distributions A vs NN (for the theta_Bn common-offset argument): "
    f"A median {np.median(col(A,'sza')):.1f}, NN median {np.median(col(NN,'sza')):.1f}")
log("")

# leave-one-out robustness of the two headline A-vs-NN axes
log("RUN31a - LEAVE-ONE-OUT robustness of the headline A-vs-NN axes")
for axis in ("cone", "n_bg"):
    full = stats.mannwhitneyu(col(A, axis), col(NN, axis), alternative="two-sided").pvalue
    dfull = np.median(col(A, axis)) - np.median(col(NN, axis))
    log(f"  {axis}: full A-vs-NN dMed={dfull:+.2f}, p={full:.4g}")
    # drop each probe (track the LEAST significant = largest p across drops)
    worst_p = (None, -1.0, 0.0)
    for p in up:
        Ai = [r for r in A if r["eid"].split("_")[1] != p]
        Ni = [r for r in NN if r["eid"].split("_")[1] != p]
        if len(Ai) < 10:
            continue
        d = np.median(col(Ai, axis)) - np.median(col(Ni, axis))
        pv = stats.mannwhitneyu(col(Ai, axis), col(Ni, axis), alternative="two-sided").pvalue
        if pv > worst_p[1]:
            worst_p = (p, pv, d)
    log(f"    worst single-probe drop: omit {worst_p[0]} -> dMed={worst_p[2]:+.2f}, p={worst_p[1]:.4g}")
    # drop each year (track the largest p across drops)
    worst_y = (None, -1.0, 0.0)
    for y in uy:
        Ai = [r for r in A if int(r["eid"][:4]) != y]
        Ni = [r for r in NN if int(r["eid"][:4]) != y]
        if len(Ai) < 10:
            continue
        d = np.median(col(Ai, axis)) - np.median(col(Ni, axis))
        pv = stats.mannwhitneyu(col(Ai, axis), col(Ni, axis), alternative="two-sided").pvalue
        if pv > worst_y[1]:
            worst_y = (y, pv, d)
    log(f"    worst single-year drop:  omit {worst_y[0]} -> dMed={worst_y[2]:+.2f}, p={worst_y[1]:.4g}")
log("")

# 31c: Jelinek bow-shock normal tilt vs polar angle (analytic)
log("RUN31c - theta_Bn vs cone: Jelinek shock-normal tilt from the Sun-Earth line (analytic)")
R0, EPS, LAM = 15.02, 6.55, 1.17
def rbs(theta, dp=2.4):
    return 2 * R0 * dp ** (-1.0 / EPS) / (np.cos(theta) + np.sqrt(np.cos(theta)**2 + LAM**2 * np.sin(theta)**2))
def normal_tilt_deg(theta_deg, dp=2.4):
    th = np.radians(theta_deg); h = 1e-5
    r = rbs(th, dp)
    rp = (rbs(th + h, dp) - rbs(th - h, dp)) / (2 * h)
    # position and tangent in 2D (X = r cos th, R = r sin th); outward normal angle from +X axis
    px, py = r*np.cos(th), r*np.sin(th)
    tx, ty = rp*np.cos(th) - r*np.sin(th), rp*np.sin(th) + r*np.cos(th)
    nx, ny = ty, -tx                      # rotate tangent -90 -> outward (points away from origin)
    if nx*px + ny*py < 0:
        nx, ny = -nx, -ny
    return np.degrees(np.arccos(np.clip(nx/np.hypot(nx, ny), -1, 1)))
for thd in (0, 10, 14.8, 20, 26.5, 30, 40):
    log(f"  polar angle {thd:5.1f} deg -> shock-normal tilt from Sun-Earth line {normal_tilt_deg(thd):5.1f} deg")
tmax = normal_tilt_deg(26.5)
log(f"  => over the candidates' SZA range (all <30 deg, median 14.8), the shock normal lies within "
    f"~{tmax:.0f} deg of the Sun-Earth line.")
log(f"  => |theta_Bn - cone| is bounded by that tilt; and since A and NN share the SZA distribution")
log(f"     (medians 14.8 vs {np.median(col(NN,'sza')):.1f} deg) the tilt is a near-common offset, so the")
log(f"     A-vs-NN cone difference (69.3 vs 63.3) maps into theta_Bn essentially unchanged.")

import os
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run31_candidate_context")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "RUN31a_RUN31c_repr_thetabn.txt"), "w", encoding="utf-8", newline="") as f:
    f.write("\n".join(L) + "\n")
print("\nsaved -> run31_candidate_context/RUN31a_RUN31c_repr_thetabn.txt", flush=True)
