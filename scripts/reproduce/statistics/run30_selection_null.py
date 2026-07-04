import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
# RUN30 - denominator-selection null for the dense-background finding (external-review Q1/Part A).
# Question: how much of the observed dense-background preference of moment-gate passers could be
# manufactured by the D_n = n_near/n_bg gate ALONE, with no physical association?
# Null model: within Dp-tertile x cone-tertile strata of the northward contributing universe,
# permute the near-shell density across encounters (breaking any physical near<->background link),
# replay the density part of the candidate gate (0.40 <= D_n <= 0.90 AND n_near >= 2.0), and
# record the pass-vs-fail background-density gap. The independence null OVERSTATES the mechanical
# induction (real near/background densities are positively coupled within a sheath, which weakens
# the gate's dependence on the background level), so "observed > null" is a CONSERVATIVE verdict.
import io, sys, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
R28 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions\run28_per_encounter.csv")
R29 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run29_themis_omni_ratio\run29_per_encounter.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run30_selection_null")
os.makedirs(OUT, exist_ok=True)

d28 = pd.read_csv(R28)
d29 = pd.read_csv(R29)[["eid", "n_sw", "ratio_bg", "ratio_near"]]
df = d28.merge(d29, on="eid", how="left")

L = []


def log(s=""):
    print(s)
    L.append(s)


log("RUN30 - DENOMINATOR-SELECTION NULL")
log("gate-lite = the density part of the candidate screen: 0.40 <= D_n <= 0.90 AND n_near >= 2.0 cm^-3")
log("universe  = northward contributing encounters (the population the candidate gate operates on)")
log("")

# ---- universe: northward, finite densities
U = df[(df.north == True) & np.isfinite(df.n_near) & np.isfinite(df.n_bg)].copy()
U_omni = U[np.isfinite(U.ratio_bg) & np.isfinite(U.n_sw)].copy()
log(f"universe N={len(U)} (northward contributing); with usable OMNI density N={len(U_omni)}")

# ---- cross-checks: reproduce the committed group counts/medians from this same table
A = U[U.grp == "A"]
NN = U[U.grp == "NONCAND"]
ck = [
    ("A count == 60", len(A) == 60),
    ("NN count == 225", len(NN) == 225),
    ("A n_bg median ~31.1", abs(A.n_bg.median() - 31.1) < 0.15),
    ("NN n_bg median ~25.2", abs(NN.n_bg.median() - 25.2) < 0.15),
]
Ao = U_omni[U_omni.grp == "A"]; NNo = U_omni[U_omni.grp == "NONCAND"]
ck += [("A ratio_bg median ~4.47", abs(Ao.ratio_bg.median() - 4.47) < 0.03),
       ("NN ratio_bg median ~3.64", abs(NNo.ratio_bg.median() - 3.64) < 0.03)]
for name, ok in ck:
    log(f"  [{'PASS' if ok else 'FAIL'}] {name}")
assert all(ok for _, ok in ck), "cross-checks failed - inputs do not reproduce committed values"
log("")

# ---- observed gate-lite split on the REAL pairing
dn_real = U.n_near / U.n_bg
gate_real = (dn_real >= 0.40) & (dn_real <= 0.90) & (U.n_near >= 2.0)
cand_overlap = U.loc[gate_real, "grp"].isin(["A", "REJ", "BORD"]).mean()
in_gate_frac = U.loc[U.grp.isin(["A", "REJ", "BORD"]), :].index.isin(U.index[gate_real]).mean()
log("=== OBSERVED (real pairing) ===")
log(f"gate-lite passers: {int(gate_real.sum())} of {len(U)} "
    f"(committed moment candidates are 107 = A60+REJ32+BORD15; gate-lite is the density part only)")
log(f"  fraction of committed A/REJ/BORD that pass gate-lite: {in_gate_frac:.2%}")
log(f"  fraction of gate-lite passers that are committed candidates: {cand_overlap:.2%}")


def gap(tab, mask, col):
    a = tab.loc[mask, col].dropna(); b = tab.loc[~mask, col].dropna()
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan, np.nan
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    return a.median() - b.median(), a.median(), p


obs_raw, med_raw, p_raw = gap(U, gate_real, "n_bg")
gate_real_o = gate_real & np.isfinite(U.ratio_bg)
m_o = gate_real_o[U_omni.index]
obs_rat, med_rat, p_rat = gap(U_omni, m_o, "ratio_bg")
log(f"  pass-vs-fail background density gap : dMed = {obs_raw:+.2f} cm^-3 (MWU p={p_raw:.2g})")
log(f"  pass-vs-fail OMNI-ratio gap         : dMed = {obs_rat:+.2f}        (MWU p={p_rat:.2g})")
log(f"  (committed A-vs-NN references: +5.85 cm^-3 [run28]; +0.83 ratio [run29 addendum])")
log("")

# ---- strata for the permutation: Dp tertile x cone tertile within the universe
U["dp_t"] = pd.qcut(U.dp, 3, labels=False, duplicates="drop")
U["cone_t"] = pd.qcut(U.cone, 3, labels=False, duplicates="drop")
U["stratum"] = U.dp_t.astype(str) + "_" + U.cone_t.astype(str)
log(f"strata: Dp tertile x cone tertile -> {U.stratum.nunique()} cells, "
    f"sizes {sorted(U.groupby('stratum').size().tolist())}")

# ---- gate-drag exposure: which clean axes correlate with n_bg in the universe?
log("")
log("=== GATE-DRAG EXPOSURE (Spearman vs n_bg in the northward universe) ===")
log("(a density-selecting gate can indirectly induce a preference on any axis correlated with n_bg)")
for c in ("cone", "sza", "t_bg", "v_bg", "dp", "ma"):
    ok = np.isfinite(U[c]) & np.isfinite(U.n_bg)
    r, p = stats.spearmanr(U.loc[ok, c], U.loc[ok, "n_bg"])
    log(f"  {c:5s} vs n_bg: rho={r:+.3f} (p={p:.2g})")

# ---- observed gate-pass gaps on the dragged axes
AX = ["cone", "sza", "t_bg", "v_bg"]
obs_ax = {}
for c in AX:
    obs_ax[c], _, p_c = gap(U, gate_real, c)
    log(f"  observed pass-vs-fail {c} gap: {obs_ax[c]:+.2f} (MWU p={p_c:.2g})")


def run_null(strata_col, nperm, seed):
    rng = np.random.default_rng(seed)
    strat = U[strata_col].to_numpy()
    idx_by_s = {s: np.where(strat == s)[0] for s in np.unique(strat)}
    n_near = U.n_near.to_numpy(); n_bg = U.n_bg.to_numpy()
    ratio_bg = U.ratio_bg.to_numpy()
    axes = {c: U[c].to_numpy() for c in AX}
    res = {k: [] for k in ["raw", "rat"] + AX}
    for k in range(nperm):
        perm = np.arange(len(U))
        for s, ix in idx_by_s.items():
            perm[ix] = rng.permutation(ix)
        nn_p = n_near[perm]                  # near-shell density re-paired within stratum
        dn_p = nn_p / n_bg
        g = (dn_p >= 0.40) & (dn_p <= 0.90) & (nn_p >= 2.0)
        if g.sum() < 5 or (~g).sum() < 5:
            continue
        res["raw"].append(np.median(n_bg[g]) - np.median(n_bg[~g]))
        r_ok = np.isfinite(ratio_bg)
        res["rat"].append(np.median(ratio_bg[g & r_ok]) - np.median(ratio_bg[(~g) & r_ok]))
        for c in AX:
            ok = np.isfinite(axes[c])
            res[c].append(np.median(axes[c][g & ok]) - np.median(axes[c][(~g) & ok]))
    return {k: np.asarray(v) for k, v in res.items()}


def report(res, obs_map):
    for name, key, obs in obs_map:
        nl = res[key]
        lo, hi = np.percentile(nl, [2.5, 97.5])
        pct = (nl < obs).mean() * 100
        p_exceed = (nl >= obs).mean()
        log(f"  {name:22s}: null median {np.median(nl):+.2f}, 95% band [{lo:+.2f},{hi:+.2f}]  "
            f"observed {obs:+.2f} -> {pct:.1f}th pct (P[null>=obs]={p_exceed:.4f})")


# ---- permutation null, main variant (Dp x cone strata)
NPERM, SEED = 5000, 30
res_main = run_null("stratum", NPERM, SEED)
log("")
log(f"=== NULL A (independence within Dp x cone tertile strata; {NPERM} permutations, seed={SEED}) ===")
obs_map = [("raw n_bg gap (cm^-3)", "raw", obs_raw), ("OMNI-ratio gap", "rat", obs_rat)] + \
          [(f"{c} gap (dragged axis)", c, obs_ax[c]) for c in AX]
report(res_main, obs_map)
log("  (cone strata partially pin the cone composition; the cone row here is a LOWER bound on drag)")

# ---- variant B: Dp-only strata (frees the cone/SZA drag channel)
U["stratum_dp"] = U.dp_t.astype(str)
res_dp = run_null("stratum_dp", NPERM, SEED + 1)
log("")
log(f"=== NULL B (Dp-tertile strata only - full drag channel open; {NPERM} permutations) ===")
report(res_dp, obs_map)

null_raw, null_rat = res_main["raw"], res_main["rat"]
log("")
log("=== READING (discipline: descriptive, selection-limited) ===")
log("The independence null quantifies the WORST-CASE mechanical induction of a dense-background")
log("preference by the D_n gate alone (real near/background coupling weakens it). Compare the")
log("observed gate-pass gap with the null band:")
log("  - observed ABOVE the null 97.5th pct -> the dense-background association exceeds what the")
log("    density-ratio gate can manufacture; report as a robust selection-limited association.")
log("  - observed INSIDE the null band -> the association is within mechanical-induction reach;")
log("    report it only as 'consistent with the candidate definition' (downgrade the wording).")
pd.DataFrame({("nullA_" + k): pd.Series(v) for k, v in res_main.items()} |
             {("nullB_" + k): pd.Series(v) for k, v in res_dp.items()}).to_csv(
    os.path.join(OUT, "run30_null_distributions.csv"), index=False)
io.open(os.path.join(OUT, "RUN30_SELECTION_NULL.txt"), "w", encoding="utf-8", newline="").write("\n".join(L) + "\n")
print("\nwritten:", OUT)
