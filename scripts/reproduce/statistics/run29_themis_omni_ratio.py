"""run29 — THEMIS-to-OMNI ratio check (a reviewer suggestion: "people use in-situ-to-OMNI ratios
to classify region, rather than only boundary models").

Lite implementation of the suggested data-driven region check, on the contributing set:
  ratio_bg   = background sheath density / OMNI upstream density  (shock compression ~3-4x
               for real shocked sheath -> validates the background placement, answering
               "how do you know the far end is outer magnetosheath?")
  ratio_near = near-shell (membership-screened) density / OMNI upstream density
  ratio_v    = background sheath flow / OMNI flow speed (subsolar sheath ~0.2-0.4)
Groups as in run28 (A = 60 spectrally supported, REJ = 32 spectral false positives,
NN = northward non-candidates). KEY QUESTION: does the suggested OMNI-ratio classifier
separate the spectrally rejected contamination (REJ) from the clean candidates (A)?
Whatever the answer, it is a result: if it does NOT separate them, the ion spectra remain
the only working discriminator (the paper's central thesis, now tested against the
reviewer's own suggested alternative).

Inputs: omni_cache/{eid}.npz + run28_per_encounter.csv. Output: runs/run29_themis_omni_ratio/.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv
import numpy as np
from scipy import stats

CACHE = P(r"H:\0mssl\review\repair\option3\omni_cache")
R28 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions\run28_per_encounter.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run29_themis_omni_ratio")
os.makedirs(OUT, exist_ok=True)


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


with open(R28, newline="") as f:
    recs = list(csv.DictReader(f))

rows, nocache = [], 0
for r in recs:
    eid = r["eid"]
    p = os.path.join(CACHE, eid + ".npz")
    if not os.path.exists(p):
        nocache += 1
        continue
    try:
        o = np.load(p)
        nsw = np.array(o["proton_density"], float)
        vsw = np.array(o["flow_speed"], float)
    except Exception:
        nocache += 1
        continue
    nsw = nsw[np.isfinite(nsw) & (nsw > 0)]
    vsw = vsw[np.isfinite(vsw) & (vsw > 0)]
    if len(nsw) < 60:
        continue
    n_sw = float(np.median(nsw))
    v_sw = float(np.median(vsw)) if len(vsw) >= 60 else np.nan
    n_bg, n_near, v_bg = fnum(r["n_bg"]), fnum(r["n_near"]), fnum(r["v_bg"])
    rows.append(dict(
        eid=eid, grp=r["grp"], north=r["north"],
        n_sw=n_sw, v_sw=v_sw,
        ratio_bg=n_bg / n_sw if np.isfinite(n_bg) else np.nan,
        ratio_near=n_near / n_sw if np.isfinite(n_near) else np.nan,
        ratio_v=v_bg / v_sw if (np.isfinite(v_bg) and np.isfinite(v_sw)) else np.nan,
    ))


def arr(rs, k):
    a = np.array([x[k] for x in rs], float)
    return a[np.isfinite(a)]


def fm(a, dig=3):
    if len(a) == 0:
        return "--"
    return f"{np.median(a):.{dig}g} [{np.percentile(a,25):.{dig}g},{np.percentile(a,75):.{dig}g}]"


def mwu(a, b):
    if len(a) < 5 or len(b) < 5:
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


A = [x for x in rows if x["grp"] == "A"]
REJ = [x for x in rows if x["grp"] == "REJ"]
NN = [x for x in rows if x["grp"] == "NONCAND" and x["north"] == "True"]
ALL = rows

L = []
L.append("RUN29 — THEMIS-TO-OMNI RATIO CHECK (reviewer-suggested cross-check)")
L.append(f"encounters with usable OMNI density: {len(rows)} of {len(recs)} (missing/short cache: {len(recs)-len(rows)})")
L.append("")
L.append("=== (1) BACKGROUND VALIDATION — is the background really shocked sheath? ===")
rb = arr(ALL, "ratio_bg")
L.append(f"  n_bg / n_sw (all contributing): {fm(rb)}   (gas-dynamic shock compression ~3-4x expected)")
L.append(f"  fraction in [2, 6]: {float(np.mean((rb >= 2) & (rb <= 6))):.0%};  fraction < 1.5: {float(np.mean(rb < 1.5)):.0%}")
rv = arr(ALL, "ratio_v")
L.append(f"  v_bg / v_sw (all contributing): {fm(rv)}   (decelerated subsolar sheath flow expected ~0.2-0.4)")
L.append("")
L.append("=== (2) DOES THE OMNI-RATIO CLASSIFIER SEPARATE THE SPECTRAL FALSE POSITIVES? ===")
for k, lab in (("ratio_near", "near-shell n / n_sw"), ("ratio_bg", "background n / n_sw")):
    a, rj, nn = arr(A, k), arr(REJ, k), arr(NN, k)
    L.append(f"  {lab:24s} A {fm(a):>22s}   REJ {fm(rj):>22s}   NN {fm(nn):>20s}   MWU(A-REJ) p={mwu(a, rj):.3g}")
L.append("")
L.append("Reading: (1) calibrates the suggested classifier on the unambiguous case (background vs solar")
L.append("wind); (2) is the decisive test — if A and REJ overlap in the in-situ-to-OMNI MOMENT-density")
L.append("ratio, the suggested classifier inherits the same moment blindness the membership screen has,")
L.append("and the ion energy SPECTRA remain the only discriminator that separates the contamination")
L.append("(peak-energy ratio 1.0 vs ~24). Either outcome is reportable; neither changes the three headline results.")

txt = "\n".join(L)
print(txt, flush=True)
with open(os.path.join(OUT, "RUN29_THEMIS_OMNI_RATIO.txt"), "w") as f:
    f.write(txt + "\n")
with open(os.path.join(OUT, "run29_per_encounter.csv"), "w", newline="") as f:
    w = csv.writer(f)
    cols = ["eid", "grp", "north", "n_sw", "v_sw", "ratio_bg", "ratio_near", "ratio_v"]
    w.writerow(cols)
    for x in rows:
        w.writerow([x[c] for c in cols])
print("\nsaved -> run29_themis_omni_ratio/RUN29_THEMIS_OMNI_RATIO.txt + run29_per_encounter.csv", flush=True)
