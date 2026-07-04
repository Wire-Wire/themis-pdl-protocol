"""run28 — candidate-vs-baseline characterisation.

Review question: "look at the solar-wind conditions
during your candidate PDLs; compare and contrast the subset with the set as a whole;
what is the difference between a PDL and the normal magnetosheath; under what conditions
does it typically occur; is your mean density significantly different at the two-sigma level?"

Populations (committed tables; nothing re-selected here):
  A    = 60 spectrally supported candidates  (run13 spec_status == SHEATH_CONSISTENT)
  REJ  = 32 spectrally rejected              (HOT_BOUNDARY_FLAG 28 + SHAPE_FLAG 4)
  BORD = 15 borderline                       (AMBIGUOUS_SPEC)
  NN   = northward contributing NON-candidates (B_z>0, not among the 107)  — PRIMARY baseline
         (B_z>0 was a candidate-selection condition, so the baseline is IMF-matched)
  NA   = all contributing non-candidates (661 - 107)                        — secondary

Axis discipline (selection-limited scope rules):
  CLEAN axes  (not selection inputs): D_p, cone, |clock|-from-north, M_A, SZA, model sheath
              thickness / MP standoff, background n/|B|/T/v  -> A-vs-baseline reportable
              (descriptive regime of where candidates are FOUND; never causal; D_p note below)
  CIRCULAR    (selection inputs: 0.40<=D_n<=0.90, beta in [0.1,2], 50<=T<=800 eV, E_B>1.2,
              n_near>=2, B_z>0): reported ONLY as A-vs-REJ (both passed the identical screen)
  DISCRIMINANT (peak-energy ratio / shape corr / flux ratio): the separation itself (A vs REJ)

Stats: two-sided Mann-Whitney U + median difference with bootstrap 95% CI ("two-sigma");
~12 clean-axis tests -> exploratory; Bonferroni threshold quoted alongside nominal p.

Inputs (committed): run10 funnel_contributing.csv, run13 validation_status.csv,
run19 contrast_per_encounter.csv, frozen substrate (background state via run7 load_enc).
Output: runs/run28_candidate_conditions/ (TXT + per-encounter CSV + 2 figures).
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from run_hardening_checks import load_enc
from psub import pmap

RUNS = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs")
OUT = os.path.join(RUNS, "run28_candidate_conditions")
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(28)
NBOOT = 5000


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def shear_from_north(clock_deg):
    """Angular distance of the IMF clock angle from due north (0 = north, 180 = south)."""
    return abs(((clock_deg + 180.0) % 360.0) - 180.0)


# ---------- 1. committed tables ----------
contrib = read_csv(os.path.join(RUNS, "run10_selection", "funnel_contributing.csv"))   # 661
valid = read_csv(os.path.join(RUNS, "run13_validation", "validation_status.csv"))      # 107
contrast_tbl = read_csv(os.path.join(RUNS, "run19_contrast_checks", "contrast_per_encounter.csv"))      # r_nB
spec = {r["eid"]: r for r in valid}
rnb = {r["eid"]: fnum(r.get("r_nB", "")) for r in contrast_tbl}

# ---------- 2. substrate pass: background state + geometry scalars ----------
def bg_rec(name, d):
    eid = name[:-4]
    e = load_enc(d)
    if e is None:
        return None
    mp0 = float(d["mp0"]); bs0 = float(d["bs0"])
    return (eid, dict(n_bg=e["n_bg"], b_bg=e["b_bg"], t_bg=e["T_bg"], v_bg=e["v_bg"],
                      mp0=mp0, bs0=bs0, thick=bs0 - mp0))

bg = dict(pmap(bg_rec, with_name=True))

# ---------- 3. join (the canonical 661 contributing set) ----------
recs, missing_bg = [], 0
for r in contrib:
    eid = r["eid"]
    v = spec.get(eid)
    if v is None:
        grp = "NONCAND"
    elif v["spec_status"] == "SHEATH_CONSISTENT":
        grp = "A"
    elif v["spec_status"] in ("HOT_BOUNDARY_FLAG", "SHAPE_FLAG"):
        grp = "REJ"
    else:
        grp = "BORD"
    b = bg.get(eid)
    if b is None:
        missing_bg += 1
        b = dict(n_bg=np.nan, b_bg=np.nan, t_bg=np.nan, v_bg=np.nan,
                 mp0=np.nan, bs0=np.nan, thick=np.nan)
    rec = dict(eid=eid, grp=grp, north=(r["bz_sign"] == "1"),
               dp=fnum(r["dp"]), cone=fnum(r["cone"]), sza=fnum(r["sza"]),
               ma=fnum(r["ma"]), bz=fnum(r["bz"]),
               shear=shear_from_north(fnum(r["clock"])),
               dn=fnum(r["Dn_mem"]), eb=fnum(r["EB_mem"]),
               beta_near=fnum(r["beta_near"]), t_near=fnum(r["T_near"]),
               n_near=fnum(r["n_near"]), r_nb=rnb.get(eid, np.nan))
    rec.update(b)
    if v is not None:
        rec.update(peak=fnum(v["peak_ratio"]), shape=fnum(v["shape_corr"]),
                   flux=fnum(v["flux_ratio"]))
    else:
        rec.update(peak=np.nan, shape=np.nan, flux=np.nan)
    recs.append(rec)

A = [r for r in recs if r["grp"] == "A"]
REJ = [r for r in recs if r["grp"] == "REJ"]
BORD = [r for r in recs if r["grp"] == "BORD"]
NN = [r for r in recs if r["grp"] == "NONCAND" and r["north"]]
NA = [r for r in recs if r["grp"] == "NONCAND"]


def arr(rs, k):
    a = np.array([r[k] for r in rs], float)
    return a[np.isfinite(a)]


def mwu_p(a, b):
    if len(a) < 5 or len(b) < 5:
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def boot_diff_ci(a, b):
    if len(a) < 5 or len(b) < 5:
        return (np.nan, np.nan)
    da = np.median(RNG.choice(a, (NBOOT, len(a))), axis=1)
    db = np.median(RNG.choice(b, (NBOOT, len(b))), axis=1)
    lo, hi = np.percentile(da - db, [2.5, 97.5])
    return float(lo), float(hi)


def fm(a, dig=3):
    if len(a) == 0:
        return "--"
    return f"{np.median(a):.{dig}g} [{np.percentile(a, 25):.{dig}g},{np.percentile(a, 75):.{dig}g}]"


CLEAN = [
    ("dp", "D_p (nPa)",
     "indirect selection: choosing depleted D_n prefers low-D_p via the Sec.7 model-coordinate ordering, itself degenerate with a Dp-correlated boundary error (Sec.7.4) -> descriptive only"),
    ("cone", "IMF cone angle (deg)", ""),
    ("shear", "|clock| from due north (deg)", "within-northward comparison (B_z>0 is a selection condition)"),
    ("bz", "IMF B_z (nT)", "within-northward comparison"),
    ("ma", "Alfven Mach number M_A", ""),
    ("sza", "solar-zenith angle (deg)", ""),
    ("mp0", "model MP standoff (R_E)", "D_p/B_z-induced (Shue parameters)"),
    ("thick", "model subsolar sheath thickness BS0-MP0 (R_E)", "D_p-induced (both standoffs scale with D_p)"),
    ("n_bg", "background sheath density (cm^-3)", "partially induced: floor n_bg >= n_near/0.90 >= 2.2 from the candidate cuts"),
    ("b_bg", "background sheath |B| (nT)", ""),
    ("t_bg", "background sheath T (eV)", ""),
    ("v_bg", "background sheath flow (km/s)", ""),
]
CIRC = [("dn", "D_n (near/bg density)"), ("eb", "E_B (near/bg field)"),
        ("beta_near", "near-shell beta"), ("t_near", "near-shell T (eV)"),
        ("n_near", "near-shell density (cm^-3)"), ("r_nb", "per-event n-|B| corr r_nB")]
DISC = [("peak", "peak-energy ratio"), ("shape", "shape correlation"), ("flux", "flux ratio")]

NTEST = len(CLEAN)
BONF = 0.05 / NTEST

L = []
L.append("RUN28 — CANDIDATE-VS-BASELINE CHARACTERISATION")
L.append(f"groups: A(spectrally supported)={len(A)}  REJ(spectral false pos)={len(REJ)}  "
         f"BORD={len(BORD)}  NN(northward non-candidate baseline)={len(NN)}  NA(all non-candidates)={len(NA)}")
L.append(f"contributing total={len(recs)} (committed run10 set); substrate background join missing for {missing_bg}")
L.append("primary test: A vs NN (IMF-matched baseline), two-sided Mann-Whitney U;")
L.append(f"median difference with bootstrap 95% CI; {NTEST} clean-axis tests -> Bonferroni p<{BONF:.4f} (nominal p also shown)")
L.append("")

# ---------- CLEAN AXES ----------
L.append("=== CLEAN AXES (not candidate-selection inputs) — where are candidates FOUND? (descriptive, selection-limited) ===")
hdr = f"{'axis':44s} {'A median[IQR]':>24s} {'NN median[IQR]':>24s} {'NA med':>8s} {'REJ med':>8s} {'p(A-NN)':>9s} {'dMed[95%CI]':>26s} sig"
L.append(hdr)
forest = []
for k, label, note in CLEAN:
    a, nn = arr(A, k), arr(NN, k)
    na, rj = arr(NA, k), arr(REJ, k)
    p = mwu_p(a, nn)
    lo, hi = boot_diff_ci(a, nn)
    dmed = np.median(a) - np.median(nn) if len(a) and len(nn) else np.nan
    sig = "**" if (np.isfinite(p) and p < BONF) else ("*" if (np.isfinite(p) and p < 0.05) else "")
    L.append(f"{label:44s} {fm(a):>24s} {fm(nn):>24s} {np.median(na) if len(na) else np.nan:>8.3g} "
             f"{np.median(rj) if len(rj) else np.nan:>8.3g} {p:>9.3g} "
             f"{dmed:>+9.3g} [{lo:>+.3g},{hi:>+.3g}] {sig}")
    if note:
        L.append(f"{'':44s} note: {note}")
    iqr_nn = np.percentile(nn, 75) - np.percentile(nn, 25) if len(nn) else np.nan
    forest.append((label, dmed, lo, hi, p, iqr_nn))
L.append("  sig: * nominal p<0.05;  ** survives Bonferroni p<%.4f" % BONF)
L.append("")

# ---------- CIRCULAR AXES ----------
L.append("=== CIRCULAR AXES (candidate-selection inputs) — A vs REJ ONLY (both passed the identical moment screen);")
L.append("    A-vs-baseline differences on these axes are BY CONSTRUCTION and are NOT findings ===")
L.append(f"{'axis':36s} {'A median[IQR]':>24s} {'REJ median[IQR]':>24s} {'p(A-REJ)':>10s}   (NN median, context only)")
for k, label in CIRC:
    a, rj, nn = arr(A, k), arr(REJ, k), arr(NN, k)
    p = mwu_p(a, rj)
    L.append(f"{label:36s} {fm(a):>24s} {fm(rj):>24s} {p:>10.3g}   ({np.median(nn) if len(nn) else np.nan:.3g})")
L.append("")

# ---------- DISCRIMINANT ----------
L.append("=== SPECTRAL DISCRIMINANT (the separation itself; A vs REJ) ===")
for k, label in DISC:
    a, rj = arr(A, k), arr(REJ, k)
    L.append(f"{label:24s} A {fm(a):>22s}   REJ {fm(rj):>22s}   MWU p={mwu_p(a, rj):.3g}")
L.append("")

# ---------- THE REVIEWER'S 2-SIGMA QUESTIONS ----------
L.append("=== THE TWO-SIGMA QUESTIONS, ANSWERED DIRECTLY ===")
dn_all = arr(recs, "dn")
w = stats.wilcoxon(dn_all - 1.0, alternative="two-sided")
L.append(f"(1) population D_n vs unity (contributing, N={len(dn_all)}): median={np.median(dn_all):.3f}, "
         f"Wilcoxon two-sided p={w.pvalue:.3g} (~{abs(stats.norm.ppf(w.pvalue/2)):.1f} sigma nominal)")
L.append("    BUT not boundary-robust: a +/-0.5 R_E MP offset moves the median across [0.55,1.14] (Appendix C) — sign not boundary-proof.")
eb_all = arr(recs, "eb")
frac = float(np.mean(eb_all > 1))
bt = stats.binomtest(int((eb_all > 1).sum()), len(eb_all), 0.5, alternative="two-sided")
L.append(f"(2) field enhancement: {frac:.0%} of contributing have E_B>1 (sign test p={bt.pvalue:.2g}) — far beyond 2 sigma,")
L.append("    but generic to near-boundary draping/pile-up (E_B indistinguishable A-vs-REJ, below) — robustly MEASURED, not PDL-diagnostic.")
a_dn, rj_dn = arr(A, "dn"), arr(REJ, "dn")
a_eb, rj_eb = arr(A, "eb"), arr(REJ, "eb")
L.append(f"(3) within the SAME moment screen (legitimate, not by construction): D_n A vs REJ = "
         f"{np.median(a_dn):.3f} vs {np.median(rj_dn):.3f} (MWU p={mwu_p(a_dn, rj_dn):.3g}); "
         f"E_B A vs REJ = {np.median(a_eb):.2f} vs {np.median(rj_eb):.2f} (MWU p={mwu_p(a_eb, rj_eb):.3g})")
L.append(f"(4) A vs NN on D_n: {np.median(a_dn):.2f} vs {np.median(arr(NN, 'dn')):.2f} — BY CONSTRUCTION")
L.append("    (the candidate screen requires 0.40<=D_n<=0.90); quoted only to answer the question, never as a finding.")
L.append("")

# ---------- CROSS-CHECKS vs committed runs ----------
L.append("=== CROSS-CHECKS vs committed outputs (pipeline-consistency gate) ===")
ck = []
ck.append(("A D_n median ~0.68 (Table 3)", abs(np.median(a_dn) - 0.68) < 0.01))
ck.append(("A E_B median ~2.36 (Table 3)", abs(np.median(a_eb) - 2.36) < 0.02))
ck.append(("E_B A-vs-REJ MWU p ~0.08 (run22)", abs(mwu_p(a_eb, rj_eb) - 0.08) < 0.04))
a_r, rj_r = arr(A, "r_nb"), arr(REJ, "r_nb")
ck.append(("r_nB A-vs-REJ MWU p ~0.80 (run21)", abs(mwu_p(a_r, rj_r) - 0.80) < 0.15))
ck.append(("all 60 have cone >= 45 deg (App C)", float(np.min(arr(A, "cone"))) >= 45.0))
ck.append(("population D_n Wilcoxon p ~0.013 (run21)", abs(w.pvalue - 0.013) < 0.006))
for label, ok in ck:
    L.append(f"  [{'PASS' if ok else 'FAIL'}] {label}")
L.append("")

# ---------- READING ----------
L.append("=== READING (selection-limited, descriptive) ===")
L.append("These are the conditions under which clean candidates are FOUND in this selection-limited sample —")
L.append("descriptive distributions, not occurrence rates and not physical drivers (the D_p axis is degenerate")
L.append("with a D_p-correlated boundary error, Sec.7.4; an IMF-axis difference must be read beside the Sec.7.3")
L.append("attenuation result). The PDL-specific separation remains spectral: peak-energy ratio 1.0 vs ~24.")

txt = "\n".join(L)
print(txt, flush=True)
with open(os.path.join(OUT, "RUN28_CANDIDATE_CONDITIONS.txt"), "w") as f:
    f.write(txt + "\n")

# ---------- per-encounter CSV (provenance) ----------
cols = ["eid", "grp", "north", "dp", "cone", "shear", "bz", "ma", "sza", "mp0", "bs0", "thick",
        "n_bg", "b_bg", "t_bg", "v_bg", "dn", "eb", "beta_near", "t_near", "n_near", "r_nb",
        "peak", "shape", "flux"]
with open(os.path.join(OUT, "run28_per_encounter.csv"), "w", newline="") as f:
    wcsv = csv.writer(f)
    wcsv.writerow(cols)
    for r in recs:
        wcsv.writerow([r.get(c, "") for c in cols])

# ---------- figures ----------
plt.rcParams.update({"font.family": "serif", "font.size": 9})

# Fig 1: upstream / clean-axis distributions, A vs NN
panels = [("dp", "D_p (nPa)", (0, 8)), ("cone", "cone angle (deg)", (0, 90)),
          ("shear", "|clock| from north (deg)", (0, 180)), ("ma", "M_A", (0, 25)),
          ("t_bg", "background T (eV)", (0, 600)), ("thick", "sheath thickness (R_E)", (1, 6))]
fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6))
for ax, (k, label, xlim) in zip(axes.flat, panels):
    a, nn = arr(A, k), arr(NN, k)
    bins = np.linspace(xlim[0], xlim[1], 25)
    ax.hist(nn, bins=bins, density=True, color="#4C72B0", alpha=0.45,
            label=f"northward non-candidates (N={len(nn)})")
    ax.hist(a, bins=bins, density=True, histtype="step", lw=2.0, color="#C44E52",
            label=f"spectrally supported (N={len(a)})")
    ax.axvline(np.median(nn), color="#4C72B0", ls="--", lw=1)
    ax.axvline(np.median(a), color="#C44E52", ls="-", lw=1.4)
    p = mwu_p(a, nn)
    ax.set_xlabel(label)
    ax.set_title(f"MWU p = {p:.3g}", fontsize=8)
    ax.set_xlim(*xlim)
axes.flat[0].set_ylabel("density")
axes.flat[3].set_ylabel("density")
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=8)
fig.suptitle("Where the spectrally supported candidates are FOUND — descriptive, selection-limited "
             "(not occurrence rates; D_p axis degenerate with boundary error)", fontsize=9.5)
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
fig.savefig(os.path.join(OUT, "fig_run28_conditions.png"), dpi=170)

# Fig 2: forest plot, robust standardized median differences (A - NN)
fig2, ax2 = plt.subplots(figsize=(7.6, 0.45 * len(forest) + 1.6))
ys = np.arange(len(forest))[::-1]
for y, (label, d, lo, hi, p, iqr) in zip(ys, forest):
    if not (np.isfinite(d) and np.isfinite(iqr) and iqr > 0):
        continue
    eff, elo, ehi = d / iqr, lo / iqr, hi / iqr
    sig = np.isfinite(p) and p < 0.05
    bonf = np.isfinite(p) and p < BONF
    col = "#C44E52" if bonf else ("#DD8452" if sig else "0.45")
    ax2.errorbar(eff, y, xerr=[[eff - elo], [ehi - eff]], fmt="o", color=col,
                 capsize=3, ms=5, lw=1.6)
    ax2.text(ax2.get_xlim()[1], y, f" p={p:.2g}", va="center", fontsize=7.5)
ax2.axvline(0, color="0.3", lw=1, ls="--")
ax2.set_yticks(ys)
ax2.set_yticklabels([f[0] for f in forest], fontsize=8.5)
ax2.set_xlabel("median difference (A − northward non-candidates), in units of baseline IQR")
ax2.set_title("Candidate-vs-baseline differences on clean (non-selection) axes\n"
              "red = survives Bonferroni; orange = nominal p<0.05; grey = n.s.", fontsize=9)
fig2.tight_layout()
fig2.savefig(os.path.join(OUT, "fig_run28_forest.png"), dpi=170)

print("\nsaved -> run28_candidate_conditions/RUN28_CANDIDATE_CONDITIONS.txt + run28_per_encounter.csv "
      "+ fig_run28_conditions.png + fig_run28_forest.png", flush=True)
