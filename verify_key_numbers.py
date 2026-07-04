"""30-second verification: re-derive the dissertation's headline numbers from committed_outputs/ alone.
No spacecraft data needed. Exit code 0 = all PASS."""
import csv, os, sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, "committed_outputs")


def rows(*p):
    with open(os.path.join(C, *p), newline="") as f:
        return list(csv.DictReader(f))


checks = []
# 1. the coordinate result: paired fixed-axis vs radial (the headline)
pr = rows("run23_paired", "paired_1d_vs_radial.csv")
d1 = np.array([float(r["Dn_1d"]) for r in pr]); dr = np.array([float(r["Dn_radial"]) for r in pr])
w = stats.wilcoxon(dr - d1, alternative="two-sided")
checks += [("paired N = 672", len(pr) == 672),
           ("fixed-axis median D_n = 0.689", abs(np.median(d1) - 0.689) < 0.001),
           ("radial median D_n = 1.114", abs(np.median(dr) - 1.114) < 0.001),
           ("median paired shift = +0.184", abs(np.median(dr - d1) - 0.184) < 0.001),
           ("Wilcoxon p ~ 4e-44", 1e-45 < w.pvalue < 1e-42)]
# 2. funnel + population
fc = rows("run10_selection", "funnel_contributing.csv")
dn = np.array([float(r["Dn_mem"]) for r in fc]); eb = np.array([float(r["EB_mem"]) for r in fc])
checks += [("contributing N = 661", len(fc) == 661),
           ("northward = 332", sum(r["bz_sign"] == "1" for r in fc) == 332),
           ("population D_n = 0.92", abs(np.median(dn) - 0.917) < 0.005),
           ("population E_B = 1.97", abs(np.median(eb) - 1.970) < 0.005),
           ("D_n vs unity Wilcoxon p = 0.013", abs(stats.wilcoxon(dn - 1).pvalue - 0.013) < 0.005)]
# 3. spectral tiers
vs = rows("run13_validation", "validation_status.csv")
tier = [r["spec_status"] for r in vs]
checks += [("candidates = 107", len(vs) == 107),
           ("spectrally supported = 60", tier.count("SHEATH_CONSISTENT") == 60),
           ("hot magnetospheric/boundary-layer = 28, shape = 4, borderline = 15",
            (tier.count("HOT_BOUNDARY_FLAG"), tier.count("SHAPE_FLAG"), tier.count("AMBIGUOUS_SPEC")) == (28, 4, 15))]
# 4. the 60-set conditional contrast
a = [r for r in vs if r["spec_status"] == "SHEATH_CONSISTENT"]
checks += [("60-set D_n = 0.68", abs(np.median([float(r["Dn_mem"]) for r in a]) - 0.676) < 0.005),
           ("60-set E_B = 2.36", abs(np.median([float(r["EB_mem"]) for r in a]) - 2.36) < 0.01)]


# 5. environment + stability (runs 27, 29 addendum, 30, 31)
def txt(*p):
    import io
    return io.open(os.path.join(C, *p), encoding="utf-8").read()

t27 = txt("run27_boundary_motion_qc", "RUN27_BOUNDARY_MOTION_QC.txt")
t29 = txt("run29_themis_omni_ratio", "RUN29_THEMIS_OMNI_RATIO.txt")
t30 = txt("run30_selection_null", "RUN30_SELECTION_NULL.txt")
t31 = txt("run31_candidate_context", "RUN31a_RUN31c_repr_thetabn.txt")
t31b = txt("run31_candidate_context", "RUN31b_CONE_STABILITY.txt")
import re as _re
checks += [
    ("run27: within-window MP motion median 1.45 R_E", "1.45" in t27),
    ("run29 addendum: A-NN background ratio +0.83, p=3.9e-05",
     ("+0.83" in t29) and ("3.92e-05" in t29 or "3.9e-05" in t29)),
    ("run30: raw density gap INSIDE the selection null (both strata)",
     _re.search(r"raw n_bg gap.*P\[null>=obs\]=0\.70", t30) is not None and
     _re.search(r"raw n_bg gap.*P\[null>=obs\]=0\.59", t30) is not None),
    ("run30: cone gap EXCEEDS the null (0.0070 / 0.0242)",
     "P[null>=obs]=0.0070" in t30 and "P[null>=obs]=0.0242" in t30),
    ("run31a: 60 candidates span 5 probes, 14 years",
     "years: 14 distinct" in t31 and all(k in t31 for k in ["tha:", "thb:", "thc:", "thd:", "the:"])),
    ("run31a: cone preference survives every leave-one-out (worst p = 0.044)",
     "p=0.04361" in t31 and "p=0.02066" in t31),
    ("run31b: preference strengthens on steady upstream (+11.2 deg, p = 0.020)",
     "+11.2" in t31b and "0.01986" in t31b),
]

fail = 0
print("VERIFYING DISSERTATION HEADLINE NUMBERS FROM committed_outputs/\n")
for label, ok in checks:
    print(("  PASS  " if ok else "  FAIL  ") + label)
    fail += (not ok)
print(f"\n{len(checks) - fail}/{len(checks)} checks passed")
sys.exit(1 if fail else 0)
