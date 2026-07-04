#!/usr/bin/env python
"""
run24 — Pi-2024 reconciliation: does the fixed-axis->radial paired shift in D_n
grow with solar-zenith angle (SZA, off-subsolar) and vanish near the subsolar line?

Joins the committed 672-row paired table (run23) to the frozen per-encounter
substrate (which carries per-sample sza_deg) and bins the paired shift
  dDn = Dn_radial - Dn_1d
by each encounter's median SZA. numbers-traceable: SZA is the committed substrate field.

Geometric expectation: fixed-axis and radial coincide at the subsolar point
(SZA->0) and diverge toward the flanks, so |dDn| should be ~0 at low SZA and
grow with SZA. If borne out, this reconciles with Pi et al. (2024), whose
construction comparison is restricted to near-subsolar geometry.
"""
import os, csv
import numpy as np

HERE = os.path.dirname(__file__)
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "01_CURRENT__rebuild"))
RUNS = os.path.join(ROOT, "runs")
SUB  = os.path.join(ROOT, "substrate")
PAIRED = os.path.join(RUNS, "run23_paired", "paired_1d_vs_radial.csv")
OUTDIR = os.path.join(RUNS, "run24_sza_reconciliation")
os.makedirs(OUTDIR, exist_ok=True)

rng = np.random.default_rng(24)

def boot_ci(a, nboot=5000):
    a = np.asarray(a, float)
    if len(a) < 3:
        return (np.nan, np.nan)
    meds = [np.median(rng.choice(a, len(a), replace=True)) for _ in range(nboot)]
    return (float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5)))

rows = list(csv.DictReader(open(PAIRED)))
out = []
missing = 0
for r in rows:
    eid = r["encounter_id"]
    p = os.path.join(SUB, eid + ".npz")
    if not os.path.exists(p):
        missing += 1
        continue
    d = np.load(p, allow_pickle=True)
    sza = np.asarray(d["sza_deg"], float)
    sza = sza[np.isfinite(sza)]
    if sza.size == 0:
        missing += 1
        continue
    out.append(dict(eid=eid,
                    sza_med=float(np.median(sza)),
                    Dn_1d=float(r["Dn_1d"]),
                    Dn_radial=float(r["Dn_radial"]),
                    dDn=float(r["Dn_radial"]) - float(r["Dn_1d"]),
                    EB_1d=float(r["EB_1d"]),
                    EB_radial=float(r["EB_radial"])))

print(f"joined {len(out)} / {len(rows)} paired encounters to substrate ({missing} missing)")

# write per-encounter joined table
with open(os.path.join(OUTDIR, "sza_shift.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["eid","sza_med","Dn_1d","Dn_radial","dDn","EB_1d","EB_radial"])
    w.writeheader()
    for o in out:
        w.writerow(o)

sza = np.array([o["sza_med"] for o in out])
dDn = np.array([o["dDn"] for o in out])
d1  = np.array([o["Dn_1d"] for o in out])
dr  = np.array([o["Dn_radial"] for o in out])

edges = [0, 10, 20, 30, 45, 60, 90]
lines = []
lines.append("run24 — paired D_n shift (radial - fixed-axis) vs solar-zenith angle")
lines.append(f"joined {len(out)}/{len(rows)} encounters; SZA = median per-encounter substrate sza_deg")
lines.append("")
lines.append(f"{'SZA bin (deg)':>14} | {'N':>4} | {'med SZA':>7} | {'med dDn':>8} | {'95% CI':>18} | {'med Dn_1d':>9} | {'med Dn_rad':>10}")
lines.append("-"*92)
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (sza >= lo) & (sza < hi)
    if m.sum() == 0:
        continue
    ci = boot_ci(dDn[m])
    lines.append(f"{lo:>6}-{hi:<7} | {m.sum():>4} | {np.median(sza[m]):>7.1f} | {np.median(dDn[m]):>8.3f} | "
                 f"[{ci[0]:>+7.3f},{ci[1]:>+7.3f}] | {np.median(d1[m]):>9.3f} | {np.median(dr[m]):>10.3f}")
lines.append("-"*92)
# Spearman of dDn vs SZA (monotonic growth test)
from scipy.stats import spearmanr
rho, pval = spearmanr(sza, dDn)
lines.append(f"Spearman(SZA, dDn) = {rho:+.3f}  (p = {pval:.2e})   [positive => shift grows off-subsolar]")
# near-subsolar vs flank
near = dDn[sza < 20]; flank = dDn[sza >= 45]
lines.append(f"near-subsolar (SZA<20, N={near.size}): median dDn = {np.median(near):+.3f}" if near.size else "near-subsolar: none")
lines.append(f"flank (SZA>=45, N={flank.size}): median dDn = {np.median(flank):+.3f}" if flank.size else "flank: none")
lines.append(f"overall: median dDn = {np.median(dDn):+.3f}, median Dn_1d={np.median(d1):.3f}, median Dn_radial={np.median(dr):.3f}")

txt = "\n".join(lines)
open(os.path.join(OUTDIR, "SZA_RECONCILIATION.txt"), "w").write(txt + "\n")
print(txt)
print("\nwrote", os.path.join(OUTDIR, "sza_shift.csv"), "and SZA_RECONCILIATION.txt")
