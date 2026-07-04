"""run31b analysis - does the A-vs-NN cone preference survive on the upstream-steady subset?

Uses run31_cone_stability.csv (within-window cone(t) median + IQR, re-fetched) merged with run28.
- validation: re-fetched within-window cone median vs the stored encounter-mean cone (run28).
- steady subset: encounters with small within-window cone IQR (well-defined IMF direction);
  re-test the A-vs-NN cone difference there. Survival => the cone preference is not an artefact
  of OMNI direction instability.
- also: are candidates themselves found under steadier or less-steady upstream IMF?
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import sys, csv
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run31_candidate_context")
CS = OUT + r"\run31_cone_stability.csv"
with open(CS, newline="") as f:
    rows = [r for r in csv.DictReader(f)]

def fnum(x):
    try:
        v = float(x); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan

for r in rows:
    r["cone28"] = fnum(r["cone28"]); r["cone_omni_med"] = fnum(r["cone_omni_med"]); r["cone_iqr"] = fnum(r["cone_iqr"])
rows = [r for r in rows if np.isfinite(r["cone_iqr"]) and np.isfinite(r["cone28"])]
A = [r for r in rows if r["grp"] == "A"]
NN = [r for r in rows if r["grp"] == "NONCAND"]

L = []
def log(s=""):
    print(s, flush=True); L.append(s)

log("RUN31b - OMNI-STABLE CONE RE-CHECK")
log(f"  coverage: {len(rows)} northward encounters with re-fetched cone(t)  (A={len(A)}, NN={len(NN)})")
# validation of the fetch against the stored encounter-mean cone
d = np.array([abs(r["cone28"] - r["cone_omni_med"]) for r in rows if np.isfinite(r["cone_omni_med"])])
log(f"  fetch validation |stored cone - re-fetched within-window median cone|: median {np.median(d):.1f} deg, 90th {np.percentile(d,90):.1f} deg")
log("")

ac = np.array([r["cone28"] for r in A]); nc = np.array([r["cone28"] for r in NN])
full_p = stats.mannwhitneyu(ac, nc, alternative="two-sided").pvalue
log(f"  full A-vs-NN cone: {np.median(ac):.1f} vs {np.median(nc):.1f}  dMed {np.median(ac)-np.median(nc):+.1f}  p={full_p:.4g}")
log("")

iqr_all = np.array([r["cone_iqr"] for r in rows])
log(f"  within-window cone IQR over the universe: median {np.median(iqr_all):.1f} deg, IQR [{np.percentile(iqr_all,25):.1f},{np.percentile(iqr_all,75):.1f}]")
# are candidates found under steadier IMF?
ai = np.array([r["cone_iqr"] for r in A]); ni = np.array([r["cone_iqr"] for r in NN])
log(f"  A vs NN cone-IQR (is upstream steadier for candidates?): {np.median(ai):.1f} vs {np.median(ni):.1f}, "
    f"MWU p={stats.mannwhitneyu(ai, ni, alternative='two-sided').pvalue:.3g}")
log("")

log("  re-test the cone preference on upstream-STEADY subsets (small within-window cone IQR):")
log(f"  {'subset':28s} {'N(A)':>5s} {'N(NN)':>6s} {'A cone':>7s} {'NN cone':>8s} {'dMed':>6s} {'p':>9s}")
for label, thr in [("all (no stability cut)", np.inf),
                   ("cone IQR <= universe median", float(np.median(iqr_all))),
                   ("cone IQR <= 15 deg", 15.0),
                   ("cone IQR <= 10 deg", 10.0)]:
    As = [r for r in A if r["cone_iqr"] <= thr]; Ns = [r for r in NN if r["cone_iqr"] <= thr]
    if len(As) < 8 or len(Ns) < 8:
        log(f"  {label:28s} {len(As):>5d} {len(Ns):>6d}   (too few)")
        continue
    a = np.array([r["cone28"] for r in As]); n = np.array([r["cone28"] for r in Ns])
    p = stats.mannwhitneyu(a, n, alternative="two-sided").pvalue
    log(f"  {label:28s} {len(As):>5d} {len(Ns):>6d} {np.median(a):>7.1f} {np.median(n):>8.1f} {np.median(a)-np.median(n):>+6.1f} {p:>9.4g}")
log("")
log("  Reading: if the A-vs-NN cone difference stays ~+6 deg and nominally significant on the")
log("  steady-IMF subset, the quasi-perpendicular preference is not an artefact of OMNI direction")
log("  instability; it is carried by encounters with a well-defined upstream IMF orientation.")

import os
with open(os.path.join(OUT, "RUN31b_CONE_STABILITY.txt"), "w", encoding="utf-8", newline="") as f:
    f.write("\n".join(L) + "\n")
print("\nsaved -> run31_candidate_context/RUN31b_CONE_STABILITY.txt", flush=True)
