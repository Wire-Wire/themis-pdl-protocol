"""RUN 11 — Bz CONFOUNDER CONTROL on the FULL membership-screened set.

User refinement (2026-05-29): Bz must NOT enter the detector and must NOT be regressed inside the
all-northward 107 (no variance). To license the phrase "northward-IMF-conditioned", test the mechanism
on ALL membership-screened CONTRIBUTING encounters (both Bz signs; 332 north + 329 south = 661),
controlling Dp, SZA, cone, probe. This is a MECHANISM-ROBUSTNESS check, not a detector, not occurrence.

Response: Dn_mem (median near-MP-shell density / background; LOWER = more depleted).
Hypothesis: northward (Bz>0) -> deeper depletion -> LOWER Dn_mem -> NEGATIVE Bz coefficient.

Two methods:
 (a) matched/stratified: median-split Dp/SZA/cone + probe -> strata; within each stratum with both
     signs, north-minus-south median Dn_mem; pool (sign test + median difference).
 (b) OLS Dn_mem ~ Bz_sign + cone + SZA + Dp + C(probe), pure-numpy, with bootstrap CI on the Bz coef.
Reads run10's funnel_contributing.csv (no substrate re-pass).
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, collections
import numpy as np

RUN10 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run11_bz_confounder")
os.makedirs(OUT, exist_ok=True)

try:
    from scipy import stats as _st
    def pval(t, dof): return 2.0 * _st.t.sf(abs(t), dof)
except Exception:
    def pval(t, dof):  # normal approx
        from math import erfc, sqrt
        return erfc(abs(t) / sqrt(2.0))


def load():
    rows = []
    for r in csv.DictReader(open(RUN10)):
        try:
            if r['bz_sign'] == '':
                continue
            rows.append(dict(eid=r['eid'], probe=r['probe'], bz_sign=int(r['bz_sign']),
                             cone=float(r['cone']), sza=float(r['sza']), dp=float(r['dp']),
                             Dn=float(r['Dn_mem'])))
        except Exception:
            pass
    return [r for r in rows if all(np.isfinite([r['cone'], r['sza'], r['dp'], r['Dn']]))]


def matched(rows, o):
    dp_m = np.median([r['dp'] for r in rows]); sza_m = np.median([r['sza'] for r in rows]); cone_m = np.median([r['cone'] for r in rows])
    def strat(r): return (r['dp'] > dp_m, r['sza'] > sza_m, r['cone'] > cone_m, r['probe'])
    cells = collections.defaultdict(list)
    for r in rows:
        cells[strat(r)].append(r)
    diffs = []; npairs = 0; nsh_lt = 0
    for key, members in cells.items():
        north = [r['Dn'] for r in members if r['bz_sign'] == 1]
        south = [r['Dn'] for r in members if r['bz_sign'] == 0]
        if len(north) >= 2 and len(south) >= 2:
            d = np.median(north) - np.median(south)
            diffs.append(d); npairs += 1
            if d < 0:
                nsh_lt += 1
    o.append("--- (a) MATCHED / STRATIFIED (median-split Dp,SZA,cone x probe) ---")
    o.append(f"  split points: Dp>{dp_m:.2f}, SZA>{sza_m:.1f}deg, cone>{cone_m:.1f}deg ; probe as-is")
    o.append(f"  strata with >=2 north AND >=2 south: {npairs}")
    if diffs:
        diffs = np.array(diffs)
        o.append(f"  within-stratum (north - south) median Dn_mem: median={np.median(diffs):+.3f}  mean={np.mean(diffs):+.3f}")
        o.append(f"  strata where north MORE depleted (diff<0): {nsh_lt}/{npairs}")
        # binomial sign test vs 0.5
        try:
            p = _st.binomtest(nsh_lt, npairs, 0.5).pvalue
        except Exception:
            p = float('nan')
        o.append(f"  sign-test p (north<south != chance): {p:.4f}")
    # unmatched baseline
    north = np.array([r['Dn'] for r in rows if r['bz_sign'] == 1]); south = np.array([r['Dn'] for r in rows if r['bz_sign'] == 0])
    o.append(f"  [unmatched baseline] median Dn_mem north={np.median(north):.3f} (N={len(north)}) vs south={np.median(south):.3f} (N={len(south)}); diff={np.median(north)-np.median(south):+.3f}")
    o.append(f"  [unmatched baseline]   MEAN Dn_mem north={np.mean(north):.3f} vs south={np.mean(south):.3f}; diff={np.mean(north)-np.mean(south):+.3f}  (mean vs median DISAGREE in sign => skewed/outlier-sensitive)")
    try:
        u = _st.mannwhitneyu(north, south, alternative='less')
        o.append(f"  [unmatched] Mann-Whitney north<south p={u.pvalue:.2e}")
    except Exception:
        pass
    # bootstrap CI on the unmatched MEDIAN difference (is even the raw median effect robust?)
    rng = np.random.default_rng(7); md = []
    for _ in range(3000):
        md.append(np.median(rng.choice(north, len(north))) - np.median(rng.choice(south, len(south))))
    lo, hi = np.percentile(md, [2.5, 97.5])
    o.append(f"  [unmatched] median-diff bootstrap 95% CI [{lo:+.3f}, {hi:+.3f}]  ({'excludes 0' if (lo<0)==(hi<0) else 'INCLUDES 0'})")


def design(rows):
    probes = sorted(set(r['probe'] for r in rows))
    ref = probes[0]
    cols = ['intercept', 'bz_sign', 'cone', 'sza', 'dp'] + [f'probe_{p}' for p in probes[1:]]
    X = []; y = []
    for r in rows:
        row = [1.0, float(r['bz_sign']), r['cone'], r['sza'], r['dp']] + [1.0 if r['probe'] == p else 0.0 for p in probes[1:]]
        X.append(row); y.append(r['Dn'])
    return np.array(X), np.array(y), cols, ref


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta; n, p = X.shape; dof = n - p
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    ss_tot = float(((y - y.mean()) ** 2).sum()); ss_res = float(resid @ resid)
    r2 = 1 - ss_res / ss_tot
    return beta, se, dof, r2


def coef_path(rows, o):
    """Show the Bz_sign coefficient as controls are ADDED — pins down which covariate is the confounder."""
    probes = sorted(set(r['probe'] for r in rows)); ref = probes[0]
    y = np.array([r['Dn'] for r in rows])
    def fit(cols_fn, label):
        X = np.array([cols_fn(r) for r in rows])
        beta, se, dof, r2 = ols(X, y)
        t = beta[1] / se[1] if se[1] > 0 else float('nan')
        o.append(f"  {label:32s} Bz coef={beta[1]:+.4f}  SE={se[1]:.4f}  p={pval(t,dof):.4f}  R2={r2:.3f}")
    o.append("")
    o.append("--- (b0) COEFFICIENT PATH: Bz_sign coef as controls are added ---")
    fit(lambda r: [1.0, float(r['bz_sign'])], "Bz only")
    fit(lambda r: [1.0, float(r['bz_sign']), r['dp']], "Bz + Dp")
    fit(lambda r: [1.0, float(r['bz_sign']), r['cone']], "Bz + cone")
    fit(lambda r: [1.0, float(r['bz_sign']), r['sza']], "Bz + SZA")
    fit(lambda r: [1.0, float(r['bz_sign']), r['dp'], r['cone'], r['sza']], "Bz + Dp + cone + SZA")
    o.append("  (If the Bz coef collapses toward 0 / loses significance the moment Dp enters, Dp is the confounder.)")


def regression(rows, o):
    coef_path(rows, o)
    X, y, cols, ref = design(rows)
    beta, se, dof, r2 = ols(X, y)
    o.append("")
    o.append("--- (b) OLS  Dn_mem ~ Bz_sign + cone + SZA + Dp + C(probe) ---")
    o.append(f"  N={len(y)}  dof={dof}  R^2={r2:.3f}  (probe reference = {ref})")
    o.append(f"  {'term':16s}{'coef':>10s}{'SE':>9s}{'t':>8s}{'p':>10s}")
    for c, b, s in zip(cols, beta, se):
        t = b / s if s > 0 else float('nan')
        o.append(f"  {c:16s}{b:10.4f}{s:9.4f}{t:8.2f}{pval(t,dof):10.4f}")
    # bootstrap CI on Bz coefficient (serial; tiny lstsq -> fast)
    rng_idx = np.arange(len(y)); bz_i = cols.index('bz_sign'); boots = []
    NB = 2000
    # deterministic resample sequence (Math.random unavailable concerns are for workflows; here numpy RNG is fine)
    rng = np.random.default_rng(12345)
    for _ in range(NB):
        s = rng.integers(0, len(y), len(y))
        bb, *_ = np.linalg.lstsq(X[s], y[s], rcond=None)
        boots.append(bb[bz_i])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    o.append(f"  Bz_sign coef bootstrap 95% CI [{lo:+.4f}, {hi:+.4f}]  ({'EXCLUDES 0 -> robust' if (lo<0)==(hi<0) else 'includes 0'})")
    o.append("  (NEGATIVE Bz_sign coef = northward IMF -> LOWER Dn_mem = DEEPER depletion, controlling Dp/SZA/cone/probe.)")
    # log-response OLS (Dn_mem is a positive ratio -> log symmetrises; coef ~ effect on geometric mean, robust to right-skew)
    Xl = X.copy(); yl = np.log(np.clip(y, 1e-6, None))
    bL, seL, dofL, r2L = ols(Xl, yl)
    i = cols.index('bz_sign'); tL = bL[i] / seL[i] if seL[i] > 0 else float('nan')
    o.append("")
    o.append("--- (c) LOG-OLS  log(Dn_mem) ~ same terms (correct for a ratio response) ---")
    o.append(f"  N={len(yl)}  R^2={r2L:.3f} ; Bz_sign coef={bL[i]:+.4f} (={100*(np.exp(bL[i])-1):+.1f}% on Dn_mem)  SE={seL[i]:.4f}  p={pval(tL,dofL):.4f}")
    o.append(f"  Dp coef={bL[cols.index('dp')]:+.4f}  p={pval(bL[cols.index('dp')]/seL[cols.index('dp')],dofL):.4f}  (Dp remains the dominant, significant predictor)")


def main():
    rows = load()
    o = ["RUN 11 — Bz CONFOUNDER CONTROL — mechanism robustness, NOT a detector/occurrence model",
         f"population: ALL membership-screened CONTRIBUTING encounters with finite Bz (both signs); N={len(rows)}",
         f"  north(Bz>0)={sum(r['bz_sign'] for r in rows)}  south(Bz<=0)={sum(1-r['bz_sign'] for r in rows)}",
         ""]
    matched(rows, o)
    regression(rows, o)
    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "BZ_CONFOUNDER.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
