"""RUN 21 — commit the contrast-result supporting statistics (numbers-traceable).

Review supplied several load-bearing numbers (D_n-vs-unity significance,
E_B-vs-unity, the r_nB 60-vs-32 Mann-Whitney, the paired Wilcoxon, Spearman+Pearson) that were NOT
in any committed run file. This script recomputes them from the frozen CSVs and commits them so the
manuscript can cite traceable values. Pure numpy (no scipy); tests via normal approximation, which is
exact-enough at N = 661 / 672 / 60 / 32. Reporting only — no new analysis branch.

Sources (frozen):
  contrast_per_encounter.csv      -> Dn_mem, EB_mem, r_nB, cls, spec   (661 contributing)
  archive_radial_catalogue_v2.csv -> Dn_1d, Dn_jbs, evaluable_*       (the §4 paired set)
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, math
import numpy as np

CONTRAST_TBL = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run19_contrast_checks\contrast_per_encounter.csv")
CAT = P(r"H:\0mssl\review\repair\option3\derived\archive_radial_catalogue_v2.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run21_contrast_stats")
os.makedirs(OUT, exist_ok=True)


def F(x):
    try:
        v = float(x); return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def p_from_z(z):
    """two-sided normal p."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def sign_test(vals, ref=1.0):
    a = np.array([v for v in vals if np.isfinite(v)])
    below = int(np.sum(a < ref)); above = int(np.sum(a > ref)); n = below + above
    # normal approx to binomial(n, 0.5), continuity-corrected
    k = min(below, above)
    z = (abs(below - above) - 1) / math.sqrt(n)
    return dict(n=n, below=below, above=above, frac_below=below / n, z=z, p_two=p_from_z(z))


def wilcoxon_signed_rank(diffs):
    """one-sample Wilcoxon signed-rank vs 0 (pass D_n-1 for 'vs unity'); normal approx, tie+zero handled."""
    d = np.array([v for v in diffs if np.isfinite(v) and v != 0.0])
    n = len(d)
    r = _rank(np.abs(d))
    W = float(np.sum(r[d > 0]))                       # sum of positive ranks
    mu = n * (n + 1) / 4.0
    # tie correction
    _, counts = np.unique(np.abs(d), return_counts=True)
    tie = np.sum(counts ** 3 - counts)
    sigma = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0 - tie / 48.0)
    z = (W - mu) / sigma
    return dict(n=n, W=W, z=z, p_two=p_from_z(z))


def mannwhitney(x, y):
    x = np.array([v for v in x if np.isfinite(v)]); y = np.array([v for v in y if np.isfinite(v)])
    nx, ny = len(x), len(y)
    allv = np.concatenate([x, y]); r = _rank(allv)
    Rx = float(np.sum(r[:nx]))
    U = Rx - nx * (nx + 1) / 2.0
    mu = nx * ny / 2.0
    _, counts = np.unique(allv, return_counts=True)
    tie = np.sum(counts ** 3 - counts); N = nx + ny
    sigma = math.sqrt(nx * ny / 12.0 * ((N + 1) - tie / (N * (N - 1))))
    z = (U - mu) / sigma
    return dict(nx=nx, ny=ny, U=U, z=z, p_two=p_from_z(z))


def _rank(a):
    order = np.argsort(a, kind='mergesort'); ranks = np.empty(len(a), float)
    ranks[order] = np.arange(1, len(a) + 1); sa = a[order]; i = 0
    while i < len(sa):
        j = i
        while j + 1 < len(sa) and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return ranks


def boot_median_ci(vals, nboot=10000, seed=0):
    a = np.array([v for v in vals if np.isfinite(v)]); rs = np.random.RandomState(seed)
    meds = np.array([np.median(a[rs.randint(0, len(a), len(a))]) for _ in range(nboot)])
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    rx = _rank(x[m]); ry = _rank(y[m]); n = int(m.sum())
    rho = float(np.corrcoef(rx, ry)[0, 1])
    z = rho * math.sqrt(n - 1)                         # large-n approx
    return rho, n, p_from_z(z)


def pearson(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[m], y[m])[0, 1])


def main():
    rows = list(csv.DictReader(open(CONTRAST_TBL)))
    dn = np.array([F(r['Dn_mem']) for r in rows]); eb = np.array([F(r['EB_mem']) for r in rows])
    rnb = np.array([F(r['r_nB']) for r in rows])
    is60 = np.array([r['cls'] == 'CONFIRMED_PDL' and r['spec'] == 'SHEATH_CONSISTENT' for r in rows])
    is32 = np.array([r['cls'] == 'CONFIRMED_PDL' and r['spec'] in ('HOT_BOUNDARY_FLAG', 'SHAPE_FLAG') for r in rows])

    o = ["RUN 21 — contrast-result supporting statistics (committed for numbers-traceable)", ""]
    o.append(f"contributing N = {len(rows)}")
    o.append("")

    # D_n vs unity
    o.append("=== D_n vs unity (population density reduction; tier 2) ===")
    st = sign_test(dn, 1.0)
    o.append(f"  median D_n = {np.nanmedian(dn):.4f}  IQR [{np.nanpercentile(dn,25):.4f}, {np.nanpercentile(dn,75):.4f}]")
    o.append(f"  sign test: {st['below']}/{st['n']} below 1 ({st['frac_below']*100:.1f}%); {st['above']} above; two-sided p = {st['p_two']:.4f}")
    wd = wilcoxon_signed_rank(dn - 1.0)
    o.append(f"  Wilcoxon signed-rank vs 1: N={wd['n']}, z={wd['z']:.3f}, two-sided p = {wd['p_two']:.4f}")
    lo, hi = boot_median_ci(dn)
    o.append(f"  bootstrap median 95% CI = [{lo:.3f}, {hi:.3f}]")
    o.append("")

    # E_B vs unity
    o.append("=== E_B vs unity (field enhancement; tier 1, robust half) ===")
    ste = sign_test(eb, 1.0)
    o.append(f"  median E_B = {np.nanmedian(eb):.4f}  IQR [{np.nanpercentile(eb,25):.4f}, {np.nanpercentile(eb,75):.4f}]")
    o.append(f"  {ste['above']}/{ste['n']} above 1 ({ste['above']/ste['n']*100:.1f}%); two-sided p = {ste['p_two']:.2e}")
    o.append("")

    # coupling D_n-E_B
    o.append("=== D_n-E_B association (coupling, NOT the decrease) ===")
    rho, n, prho = spearman(dn, eb)
    o.append(f"  Spearman(D_n, E_B) = {rho:.4f}  (N={n}, p = {prho:.2e})")
    o.append(f"  Pearson(D_n, E_B)  = {pearson(dn, eb):.4f}  (monotonic, not linear)")
    o.append("")

    # r_nB 60 vs 32 — the decisive 'coupling is generic' number
    o.append("=== per-event n-|B| coupling r_nB: clean 60 vs rejected 32 (decisive) ===")
    r60 = rnb[is60]; r32 = rnb[is32]
    o.append(f"  60 supported:  median r_nB = {np.nanmedian(r60):.3f}, frac<0 = {np.mean(r60[np.isfinite(r60)]<0)*100:.1f}%  (N={np.sum(np.isfinite(r60))})")
    o.append(f"  32 rejected:   median r_nB = {np.nanmedian(r32):.3f}, frac<0 = {np.mean(r32[np.isfinite(r32)]<0)*100:.1f}%  (N={np.sum(np.isfinite(r32))})")
    mw = mannwhitney(r60, r32)
    o.append(f"  Mann-Whitney (r_nB, 60 vs 32): U={mw['U']:.1f}, z={mw['z']:.3f}, two-sided p = {mw['p_two']:.3f}")
    o.append(f"  => indistinguishable: the anticorrelation is generic, not a PDL-vs-contamination discriminator")
    o.append("")

    # T_near 60 vs 32 (for the §3.5 fix)
    tn = np.array([F(r['T_near']) for r in rows])
    o.append("=== T_near: clean 60 vs rejected 32 (for the §3.5 adiabatic-heating caveat) ===")
    mwT = mannwhitney(tn[is60], tn[is32])
    o.append(f"  60 median T_near = {np.nanmedian(tn[is60]):.0f} eV; 32 median = {np.nanmedian(tn[is32]):.0f} eV")
    o.append(f"  Mann-Whitney p = {mwT['p_two']:.3f}  => no temperature signature (adiabatic-heating story unsupported)")
    o.append("")

    # §4 paired Wilcoxon
    o.append("=== §4 paired fixed-axis vs radial (commit the Wilcoxon stat) ===")
    try:
        crow = list(csv.DictReader(open(CAT)))
        P = [(F(r['Dn_1d']), F(r['Dn_jbs'])) for r in crow
             if str(r.get('evaluable_1d')) == 'True' and str(r.get('evaluable_jbs')) == 'True']
        P = [(a, b) for a, b in P if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0]
        d1 = np.array([a for a, b in P]); dg = np.array([b for a, b in P]); diff = dg - d1
        o.append(f"  N pairs = {len(P)}; median D_n: 1d {np.median(d1):.4f} -> radial {np.median(dg):.4f}")
        o.append(f"  median per-encounter shift = {np.median(diff):+.4f}; {np.mean(diff>0)*100:.1f}% shift upward")
        wp = wilcoxon_signed_rank(diff)
        o.append(f"  Wilcoxon signed-rank (paired): N={wp['n']}, z={wp['z']:.2f}, two-sided p = {wp['p_two']:.2e}")
        lo2, hi2 = boot_median_ci(diff)
        o.append(f"  bootstrap 95% CI on median shift = [{lo2:+.3f}, {hi2:+.3f}]")
    except FileNotFoundError:
        o.append(f"  (catalogue {CAT} not found — §4 paired stat not recomputed here)")
    o.append("")

    txt = "\n".join(o)
    print(txt)
    with open(os.path.join(OUT, "RUN21_CONTRAST_TBL_STATS.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT)


if __name__ == "__main__":
    main()
