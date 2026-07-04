import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
# run29 addendum: the A-vs-NN test on the background compression ratio (was only run A-vs-REJ).
# Appends to the committed run29 output so the number is traceable (verified_only).
import io, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
RUN = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run29_themis_omni_ratio")
df = pd.read_csv(RUN + r"\run29_per_encounter.csv")

# groups as defined in run28/run29: A = spectrally supported; NN = northward non-candidates
A = df[(df.grp == "A")]
NN = df[(df.grp == "NONCAND") & (df.north == True)]


def med_iqr(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return np.median(x), np.percentile(x, 25), np.percentile(x, 75), len(x)


def boot_dmed(a, b, n=20000, seed=29):
    rng = np.random.default_rng(seed)
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    b = np.asarray(b, float); b = b[np.isfinite(b)]
    d = np.empty(n)
    for i in range(n):
        d[i] = np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b)))
    return np.percentile(d, 2.5), np.percentile(d, 97.5)


lines = ["", "=== ADDENDUM (committed post-hoc): A-vs-NN on the OMNI-normalised ratios ===",
         "(the section-2 table above only tested A-vs-REJ; the dense-sheath finding of run28 in",
         " ratio form, n_bg/n_sw, removes solar-wind density variation; Dp and M_A are A-vs-NN null",
         " in run28, so the normalisation does not import a known group difference)"]
for col, name in [("ratio_bg", "background n_bg / n_sw"), ("ratio_near", "near-shell n / n_sw"),
                  ("ratio_v", "background v_bg / v_sw")]:
    a = A[col].astype(float); a = a[np.isfinite(a)]
    b = NN[col].astype(float); b = b[np.isfinite(b)]
    ma, q1a, q3a, na = med_iqr(a)
    mb, q1b, q3b, nb = med_iqr(b)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    lo, hi = boot_dmed(a, b)
    lines.append(f"  {name:26s}  A {ma:5.2f} [{q1a:.2f},{q3a:.2f}] (N={na})   "
                 f"NN {mb:5.2f} [{q1b:.2f},{q3b:.2f}] (N={nb})   MWU(A-NN) p={p:.3g}   dMed[95%CI] {ma-mb:+.2f} [{lo:+.2f},{hi:+.2f}]")
lines.append("Reading: descriptive, selection-limited (same discipline as run28); a significant background-")
lines.append("ratio excess in A would state the run28 dense-sheath preference in compression-ratio form")
lines.append("(stronger-than-average shock compression), independent of absolute solar-wind density.")
txt = "\n".join(lines) + "\n"
print(txt)
with io.open(RUN + r"\RUN29_THEMIS_OMNI_RATIO.txt", "a", encoding="utf-8", newline="") as f:
    f.write(txt)
print("APPENDED to RUN29_THEMIS_OMNI_RATIO.txt")
