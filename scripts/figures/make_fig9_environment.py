import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
# Figure 9 - candidate environments against the northward baseline (manuscript Section 7).
# Panel (a): A-vs-NN median differences on the clean axes, scaled by the baseline IQR, with
#            bootstrap 95% CIs and the run30 denominator-null 95% band on the six tested axes.
# Panel (b): A-vs-REJ Mann-Whitney separation (-log10 p) on environmental axes vs the three
#            spectral discriminants.
# Self-checking: recomputed medians/p-values are asserted against the committed run28/run29/run30
# outputs before the figure is written.
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
R28 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions\run28_per_encounter.csv")
R29 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run29_themis_omni_ratio\run29_per_encounter.csv")
R30 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run30_selection_null\run30_null_distributions.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\figures\fig9_environment.png")

d = pd.read_csv(R28).merge(pd.read_csv(R29)[["eid", "ratio_bg"]], on="eid", how="left")
nulls = pd.read_csv(R30)
A = d[d.grp == "A"]
REJ = d[d.grp == "REJ"]
NN = d[(d.grp == "NONCAND") & (d.north == True)]

# ---- cross-checks against committed values (run28/run29 txt) ----
ck = [
    ("A=60", len(A) == 60), ("REJ=32", len(REJ) == 32), ("NN=225", len(NN) == 225),
    ("A n_bg ~31.1", abs(A.n_bg.median() - 31.1) < 0.15),
    ("NN n_bg ~25.2", abs(NN.n_bg.median() - 25.2) < 0.15),
    ("A cone ~69.3", abs(A.cone.median() - 69.3) < 0.15),
    ("NN cone ~63.3", abs(NN.cone.median() - 63.3) < 0.15),
    ("A ratio ~4.47", abs(A.ratio_bg.median() - 4.47) < 0.03),
    ("NN ratio ~3.64", abs(NN.ratio_bg.median() - 3.64) < 0.03),
]
for name, ok in ck:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
assert all(ok for _, ok in ck), "cross-checks failed - inputs do not reproduce committed values"

# ---- panel (a): clean-axis forest, A - NN, scaled by NN IQR ----
AXES = [  # (column, label, null-csv key or None)
    ("dp",       "dynamic pressure $D_p$",        None),
    ("cone",     "IMF cone angle",                "cone"),
    ("shear",    "|clock| from due north",        None),
    ("bz",       "IMF $B_z$",                     None),
    ("ma",       "Alfvén Mach $M_A$",        None),
    ("sza",      "solar-zenith angle",            "sza"),
    ("mp0",      "model MP standoff",             None),
    ("thick",    "model sheath thickness",        None),
    ("n_bg",     "background density $n_{bg}$",   "raw"),
    ("b_bg",     "background $|B|_{bg}$",         None),
    ("t_bg",     "background $T_{bg}$",           "t_bg"),
    ("v_bg",     "background flow $v_{bg}$",      "v_bg"),
    ("ratio_bg", "background $n_{bg}/n_{sw}$",    "rat"),
]
rng = np.random.default_rng(9)


def boot_dmed(a, b, n=10000):
    a = a[np.isfinite(a)].to_numpy(); b = b[np.isfinite(b)].to_numpy()
    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b)))
    return np.percentile(out, [2.5, 97.5])


# observed pass-vs-fail screen gaps (the quantity the run30 null bounds; RUN30 txt)
PASSGAP = {"raw": +5.76, "rat": +0.74, "cone": +6.05, "sza": +0.30, "t_bg": -33.76, "v_bg": -18.43}

rows = []
for col, lab, nk in AXES:
    a, b = A[col].dropna(), NN[col].dropna()
    iqr = b.quantile(0.75) - b.quantile(0.25)
    dmed = a.median() - b.median()
    lo, hi = boot_dmed(a, b)
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    band = gap = None
    if nk is not None:  # run30 null band, envelope of strata variants A and B
        pool = np.concatenate([nulls[f"nullA_{nk}"].dropna(), nulls[f"nullB_{nk}"].dropna()])
        band = np.percentile(pool, [2.5, 97.5]) / iqr
        gap = PASSGAP[nk] / iqr
    rows.append(dict(lab=lab, d=dmed / iqr, lo=lo / iqr, hi=hi / iqr, p=p, band=band, gap=gap,
                     posthoc=(col == "ratio_bg")))
    print(f"  {col:9s} dMed={dmed:+.2f} scaled={dmed/iqr:+.2f} p={p:.3g}")

# ---- panel (b): A-vs-REJ separation, environment vs spectra ----
ENV_B = [("n_bg", "background $n_{bg}$"), ("ratio_bg", "background $n_{bg}/n_{sw}$"),
         ("eb", "field contrast $E_B$"), ("r_nb", "per-event $n$–$|B|$ corr."),
         ("dp", "dynamic pressure $D_p$"), ("cone", "IMF cone angle"), ("sza", "solar-zenith angle")]
SPEC_B = [("peak", "peak-energy ratio"), ("shape", "shape correlation"), ("flux", "flux ratio")]
committed_b = {"n_bg": 0.191, "ratio_bg": 0.271, "eb": 0.0801, "r_nb": 0.806, "dp": 0.127,
               "cone": 0.0459, "sza": 0.032, "peak": 1.27e-10, "shape": 2.23e-11, "flux": 1.73e-12}
bars = []
for col, lab in ENV_B + SPEC_B:
    a, r = A[col].dropna(), REJ[col].dropna()
    p = stats.mannwhitneyu(a, r, alternative="two-sided").pvalue
    ref = committed_b[col]
    ok = abs(np.log10(p) - np.log10(ref)) < 0.15
    print(f"  A-REJ {col:9s} p={p:.3g} (committed {ref:.3g}) [{'PASS' if ok else 'FAIL'}]")
    assert ok, f"A-vs-REJ p for {col} does not reproduce the committed value"
    bars.append(dict(lab=lab, p=p, spec=(col in dict(SPEC_B))))

# ---- draw ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 6.4), width_ratios=[1.15, 1.0])
ys = np.arange(len(rows))[::-1]
for y, r in zip(ys, rows):
    if r["band"] is not None:
        ax1.barh(y, r["band"][1] - r["band"][0], left=r["band"][0], height=0.62,
                 color="0.82", zorder=1, label="_")
    ax1.plot([r["lo"], r["hi"]], [y, y], color="#1f4e79", lw=2, zorder=3)
    star = "**" if r["p"] < 0.0042 else ("*" if r["p"] < 0.05 else "")
    ax1.plot(r["d"], y, "o", color="#c0392b" if star else "#1f4e79", ms=7, zorder=4)
    if r["gap"] is not None:
        ax1.plot(r["gap"], y, "D", mfc="none", mec="#2c3e50", ms=7, mew=1.4, zorder=5)
    if star:
        ax1.annotate(star, (r["d"], y + 0.18), ha="center", color="#c0392b", fontsize=11)
ax1.axhline(ys[-1] + 0.5, color="0.55", lw=0.8, ls=":")  # separate the post-hoc ratio row
ax1.axvline(0, color="0.4", lw=0.8)
ax1.set_yticks(ys); ax1.set_yticklabels([r["lab"] for r in rows], fontsize=9)
ax1.set_xlabel("median difference, spectrally supported − baseline\n(units of baseline IQR)", fontsize=9)
ax1.set_title("(a) Where are candidates found? A (60) vs northward baseline (225)\n"
              "grey band: 95% range inducible by the density selection alone (run30 null)", fontsize=9.5)

yb = np.arange(len(bars))[::-1]
cols = ["#7f8c8d" if not b["spec"] else "#c0392b" for b in bars]
vals = [-np.log10(b["p"]) for b in bars]
ax2.barh(yb, vals, color=cols, height=0.62)
ax2.axvline(-np.log10(0.05), color="0.45", lw=0.9, ls=":")
ax2.axvline(-np.log10(0.0042), color="0.25", lw=0.9, ls="--")
ax2.text(-np.log10(0.05), len(bars) - 0.32, " p=0.05", fontsize=7.5, color="0.4")
ax2.text(-np.log10(0.0042), len(bars) - 0.32, " Bonferroni", fontsize=7.5, color="0.25")
ax2.set_yticks(yb); ax2.set_yticklabels([b["lab"] for b in bars], fontsize=9)
ax2.set_xlabel("Mann–Whitney separation of supported vs rejected, −log$_{10}$ p", fontsize=9)
ax2.set_title("(b) What validates candidates? environment (grey) vs ion spectra (red)\n"
              "A (60) vs spectrally rejected (32)", fontsize=9.5)
fig.suptitle("Figure 9 — candidate environments: discovery conditions in a selection-limited sample, not occurrence rates",
             fontsize=10.5, y=0.995)
fig.tight_layout(rect=(0, 0.045, 1, 0.96))
fig.text(0.055, 0.012, "** p<0.0042 (12-test Bonferroni; the $n_{bg}/n_{sw}$ row below the dotted rule is a post-hoc run29 axis) · "
                       "* nominal · ◇ observed pass-vs-fail screen gap (the quantity the null bounds)",
         fontsize=7.6, color="0.3")
fig.savefig(OUT, dpi=200)
print("written:", OUT)
