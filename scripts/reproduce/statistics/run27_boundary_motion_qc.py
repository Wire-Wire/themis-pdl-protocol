"""run27 — within-window boundary-motion QC (boundary-motion quality control).

Motivating question: "check that your magnetopause and bow-shock positions do not change
significantly during the interval — e.g. where the magnetopause was originally may later
be the bow-shock position; exclude such events."

Context (discovered while preparing this run): the frozen substrate stores ONE set of
boundary parameters per encounter (scalar mp0/bs0/alpha from encounter-averaged OMNI), as
the manuscript states (Sec.3.2 "per encounter"). Within-window boundary motion is therefore
UNMODELLED in the pipeline, and this QC bounds its impact using the re-fetched 1-min OMNI:
  Shue-1998 subsolar standoff r0(t) = (10.22 + 1.29 tanh(0.184 (Bz+8.14))) Dp^(-1/6.6)
  Jelinek-2012 subsolar bow shock R_BS(t) = 15.02 Dp^(-1/6.55)
Metrics per encounter: range of r0(t), range of R_BS(t), region-swap flag
(max r0(t) >= min R_BS(t)), upstream coverage. Exclusion sweep at dMP > 0.5/1.0/1.5 R_E:
impact on the contributing count, the population D_n/E_B medians, and the 60-candidate set.

Cross-checks: median r0(t) vs the stored encounter scalar mp0; OMNI-median Dp vs stored dp.
Inputs: omni_cache/{eid}.npz (fetch_omni_contributing.py) + run28_per_encounter.csv (committed).
Output: runs/run27_boundary_motion_qc/ (TXT + per-encounter CSV).
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv
import numpy as np

sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
CACHE = P(r"H:\0mssl\review\repair\option3\omni_cache")
SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
R28 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions\run28_per_encounter.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run27_boundary_motion_qc")
os.makedirs(OUT, exist_ok=True)


def shue_r0(dp, bz):
    return (10.22 + 1.29 * np.tanh(0.184 * (bz + 8.14))) * dp ** (-1.0 / 6.6)


def jel_rbs(dp):
    return 15.02 * dp ** (-1.0 / 6.55)


with open(R28, newline="") as f:
    recs = list(csv.DictReader(f))


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


rows, nocache, lowcov = [], 0, 0
for r in recs:
    eid = r["eid"]
    p = os.path.join(CACHE, eid + ".npz")
    if not os.path.exists(p):
        nocache += 1
        continue
    try:
        o = np.load(p)
        dp = np.array(o["Pressure"], float)
        bz = np.array(o["BZ_GSM"], float)
    except Exception:
        nocache += 1
        continue
    fin = np.isfinite(dp) & np.isfinite(bz) & (dp > 0)
    if fin.sum() < 60:  # need >= 1 h of valid upstream minutes
        lowcov += 1
        continue
    r0 = shue_r0(dp[fin], bz[fin])
    rbs = jel_rbs(dp[fin])
    d = np.load(os.path.join(SUB, eid + ".npz"), allow_pickle=True)
    rows.append(dict(
        eid=eid, grp=r["grp"], dn=fnum(r["dn"]), eb=fnum(r["eb"]),
        cov=float(fin.mean()),
        mp_med=float(np.median(r0)), mp_rng=float(np.max(r0) - np.min(r0)),
        bs_rng=float(np.max(rbs) - np.min(rbs)),
        swap=bool(np.max(r0) >= np.min(rbs)),
        mp0_stored=float(d["mp0"]), dp_stored=float(d["dp"]),
        dp_omni_med=float(np.median(dp[fin])),
    ))

L = []
L.append("RUN27 — WITHIN-WINDOW BOUNDARY-MOTION QC")
L.append(f"contributing with usable 1-min OMNI: {len(rows)} of {len(recs)} "
         f"(no/failed cache: {nocache}; <1 h upstream coverage: {lowcov})")
L.append("")
L.append("NOTE (pipeline fact, stated in manuscript Sec.3.2): boundary parameters are evaluated PER")
L.append("ENCOUNTER from encounter-averaged OMNI; within-window boundary motion is unmodelled. This QC")
L.append("bounds its impact with re-fetched 1-min OMNI.")
L.append("")

# cross-checks
mdiff = np.array([abs(x["mp_med"] - x["mp0_stored"]) for x in rows])
ddiff = np.array([abs(x["dp_omni_med"] - x["dp_stored"]) / max(x["dp_stored"], 1e-3) for x in rows])
L.append("=== CROSS-CHECKS (fetch + formula consistency) ===")
L.append(f"  |median Shue r0(t) - stored mp0|: median {np.median(mdiff):.3f} R_E, 90th pct {np.percentile(mdiff, 90):.3f} R_E")
L.append(f"  |OMNI-median Dp - stored dp| / stored: median {np.median(ddiff):.1%}, 90th pct {np.percentile(ddiff, 90):.1%}")
L.append("")

mp_rng = np.array([x["mp_rng"] for x in rows])
bs_rng = np.array([x["bs_rng"] for x in rows])
swaps = np.array([x["swap"] for x in rows])
L.append("=== WITHIN-WINDOW BOUNDARY MOTION (model subsolar standoffs, 1-min OMNI) ===")
L.append(f"  MP standoff range dMP: median {np.median(mp_rng):.2f} R_E, IQR [{np.percentile(mp_rng,25):.2f},{np.percentile(mp_rng,75):.2f}], 90th {np.percentile(mp_rng,90):.2f}, max {np.max(mp_rng):.2f}")
L.append(f"  BS standoff range dBS: median {np.median(bs_rng):.2f} R_E, 90th {np.percentile(bs_rng,90):.2f}")
L.append(f"  REGION-SWAP encounters (max r0(t) >= min R_BS(t), the reviewer's scenario): {int(swaps.sum())} of {len(rows)}")
L.append("")

dn_all = np.array([x["dn"] for x in rows])
eb_all = np.array([x["eb"] for x in rows])
isA = np.array([x["grp"] == "A" for x in rows])
L.append("=== EXCLUSION SWEEP — drop encounters with dMP above threshold; do the headline numbers move? ===")
L.append(f"{'cut':18s} {'N kept':>7s} {'D_n med':>9s} {'E_B med':>9s} {'60-set kept':>12s} {'D_n(A) med':>11s}")
base_dn, base_eb = np.median(dn_all), np.median(eb_all)
L.append(f"{'none (baseline)':18s} {len(rows):>7d} {base_dn:>9.3f} {base_eb:>9.3f} {int(isA.sum()):>12d} {np.median(dn_all[isA]):>11.3f}")
for cut in (1.5, 1.0, 0.5):
    keep = mp_rng <= cut
    kA = keep & isA
    L.append(f"{'dMP <= ' + str(cut) + ' R_E':18s} {int(keep.sum()):>7d} {np.median(dn_all[keep]):>9.3f} "
             f"{np.median(eb_all[keep]):>9.3f} {int(kA.sum()):>12d} {np.median(dn_all[kA]) if kA.sum() else float('nan'):>11.3f}")
keep = ~swaps
kA = keep & isA
L.append(f"{'no region-swap':18s} {int(keep.sum()):>7d} {np.median(dn_all[keep]):>9.3f} "
         f"{np.median(eb_all[keep]):>9.3f} {int(kA.sum()):>12d} {np.median(dn_all[kA]) if kA.sum() else float('nan'):>11.3f}")
L.append("")
L.append("Reading: if the D_n/E_B medians and the candidate count are stable as strongly-moving-boundary")
L.append("encounters are removed, within-window boundary motion is not driving the contrast; the swap")
L.append("count answers the reviewer's specific scenario directly.")

txt = "\n".join(L)
print(txt, flush=True)
with open(os.path.join(OUT, "RUN27_BOUNDARY_MOTION_QC.txt"), "w") as f:
    f.write(txt + "\n")
with open(os.path.join(OUT, "run27_per_encounter.csv"), "w", newline="") as f:
    w = csv.writer(f)
    cols = ["eid", "grp", "dn", "eb", "cov", "mp_med", "mp_rng", "bs_rng", "swap", "mp0_stored", "dp_stored", "dp_omni_med"]
    w.writerow(cols)
    for x in rows:
        w.writerow([x[c] for c in cols])
print("\nsaved -> run27_boundary_motion_qc/RUN27_BOUNDARY_MOTION_QC.txt + run27_per_encounter.csv", flush=True)
