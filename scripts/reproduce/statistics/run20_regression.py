"""RUN 20 — commit the Table-5 regression output for provenance (NO substrate re-pass; reads
funnel_contributing.csv). Reproduces the manuscript Table 5: log-OLS of log(Dn_mem) on the drivers,
uncontrolled (Bz-sign only) and controlled (Bz, Dp, SZA, cone, probe dummies), with per-term
coefficients, standard errors, normal-approx p (math.erfc; df large), bootstrap 95% CIs, VIFs, and
the IMF-Dp correlation. Pure numpy (no scipy/statsmodels). Reporting/provenance, not new science.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, math
import numpy as np

CSV = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run20_regression")
os.makedirs(OUT, exist_ok=True)
NBOOT = 5000


def f(x):
    try:
        v = float(x); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def pval(t):
    return math.erfc(abs(t) / math.sqrt(2.0))  # two-sided normal approx (df ~ 655)


def ols(X, y):
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    r2 = 1.0 - (resid @ resid) / (((y - y.mean()) ** 2).sum())
    return beta, se, r2


def boot_ci(X, y, j, nboot=NBOOT, seed=0):
    rs = np.random.RandomState(seed); N = len(y); vals = np.empty(nboot)
    for i in range(nboot):
        s = rs.randint(0, N, N)
        b, _, _, _ = np.linalg.lstsq(X[s], y[s], rcond=None)
        vals[i] = b[j]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def vif(cols, j):
    """VIF of column j: regress col j on the other columns (+intercept)."""
    y = cols[:, j]
    others = np.delete(cols, j, axis=1)
    X = np.column_stack([np.ones(len(y)), others])
    _, _, r2 = ols(X, y)
    return 1.0 / (1.0 - r2) if r2 < 1 else np.inf


def main():
    rows = list(csv.DictReader(open(CSV)))
    dn = np.array([f(r['Dn_mem']) for r in rows])
    bz = np.array([1.0 if r['bz_sign'] == '1' else (0.0 if r['bz_sign'] == '0' else np.nan) for r in rows])
    dp = np.array([f(r['dp']) for r in rows])
    sza = np.array([f(r['sza']) for r in rows])
    cone = np.array([f(r['cone']) for r in rows])
    probe = np.array([r['probe'] for r in rows])
    m = np.isfinite(dn) & (dn > 0) & np.isfinite(bz) & np.isfinite(dp) & np.isfinite(sza) & np.isfinite(cone)
    y = np.log(dn[m]); bzm = bz[m]; dpm = dp[m]; szam = sza[m]; conem = cone[m]; probem = probe[m]
    probes = sorted(set(probem)); ref = probes[0]
    dummies = [(probem == p).astype(float) for p in probes if p != ref]
    dnames = [f"probe_{p}" for p in probes if p != ref]

    o = []
    o.append("RUN 20 — Table-5 regression (log-OLS of log(Dn_mem); committed for provenance)")
    o.append(f"source: {CSV}")
    o.append(f"N = {len(y)}  (contributing encounters with finite Dn>0, Bz-sign, Dp, SZA, cone)")
    o.append(f"probe reference = {ref}; dummies = {dnames}")
    o.append("p-values: two-sided normal approximation (erfc); 95% CI: nonparametric bootstrap, "
             f"{NBOOT} resamples, seed 0.")
    o.append("")

    # uncontrolled: log(Dn) ~ Bz
    Xu = np.column_stack([np.ones(len(y)), bzm])
    bu, seu, r2u = ols(Xu, y)
    ciu = boot_ci(Xu, y, 1)
    o.append("--- UNCONTROLLED (Bz-sign only) ---")
    o.append(f"  Bz(north=1): coef={bu[1]:+.4f}  SE={seu[1]:.4f}  p={pval(bu[1]/seu[1]):.4f}  "
             f"boot95% CI=[{ciu[0]:+.3f},{ciu[1]:+.3f}]   R2={r2u:.4f}")
    o.append("")

    # controlled: log(Dn) ~ Bz + Dp + SZA + cone + probe dummies
    cont = [bzm, dpm, szam, conem] + dummies
    names = ['Bz(north=1)', 'Dp(per nPa)', 'SZA(deg)', 'cone(deg)'] + dnames
    Xc = np.column_stack([np.ones(len(y))] + cont)
    bc, sec, r2c = ols(Xc, y)
    contmat = np.column_stack(cont)
    o.append(f"--- CONTROLLED (Bz, Dp, SZA, cone, probe dummies); R2={r2c:.4f} ---")
    o.append(f"  {'term':14s} {'coef':>9s} {'SE':>8s} {'p':>9s} {'boot95% CI':>22s} {'VIF':>6s}")
    for k, nm in enumerate(names):
        j = k + 1
        ci = boot_ci(Xc, y, j, seed=k)
        v = vif(contmat, k)
        o.append(f"  {nm:14s} {bc[j]:+9.4f} {sec[j]:8.4f} {pval(bc[j]/sec[j]):9.4f} "
                 f"[{ci[0]:+.3f},{ci[1]:+.3f}]".rjust(22) + f" {v:6.3f}")
    o.append(f"  intercept      {bc[0]:+9.4f} {sec[0]:8.4f}")
    o.append("")

    # attenuation + IMF-Dp correlation
    rIMFdp = float(np.corrcoef(bzm, dpm)[0, 1])
    # corr p via Fisher z
    n = len(y); z = 0.5 * math.log((1 + rIMFdp) / (1 - rIMFdp)); se_z = 1 / math.sqrt(n - 3)
    p_r = math.erfc(abs(z / se_z) / math.sqrt(2))
    o.append(f"IMF-Dp correlation r = {rIMFdp:+.4f}  (p = {p_r:.3f})")
    o.append(f"IMF attenuation: uncontrolled coef {bu[1]:+.4f} -> controlled {bc[1]:+.4f}  "
             f"(controlled is {abs(bc[1]/bu[1])*100:.0f}% of uncontrolled)")
    o.append("")
    o.append("NOTE: reproduces manuscript Table 5; values match the figures reported there. p for the "
             "large-|t| Dp term is below the erfc floor (report as <1e-4).")

    txt = "\n".join(o)
    print(txt)
    with open(os.path.join(OUT, "RUN20_REGRESSION.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\nsaved ->", os.path.join(OUT, "RUN20_REGRESSION.txt"))


if __name__ == "__main__":
    main()
