#!/usr/bin/env python
"""
RUN 22 - supporting statistics (committed for numbers-traceable).

Commits the two facets of a consolidated review objection, both
recomputed from the frozen per-encounter substrate (no new selection):

  FACET 1 - the near-MP FIELD enhancement is GENERIC, not PDL-diagnostic:
            E_B does not distinguish the spectrally-clean PDL set from the
            spectrally-rejected contamination, and is present without depletion.

  FACET 2 - the tier-2 density verdict ("statistically below unity") quotes
            only SAMPLING uncertainty; the boundary-placement SYSTEMATIC
            (Appendix C / RUN17) dwarfs it and flips the sub-unity sign.

Source data : runs/run19_contrast_checks/contrast_per_encounter.csv  (N=661)
Cross-refs  : runs/run21_contrast_stats/RUN21_CONTRAST_TBL_STATS.txt (tier-2 sampling stats)
              runs/run17_robustness/ROBUSTNESS.txt           (boundary band)
              runs/run2_profile_cube/RUN2_SUMMARY.txt         (sheath-wide E_B profile)
Method      : medians + Mann-Whitney U (normal approx w/ tie correction, erfc
              two-sided p) - same pure-stdlib discipline as run21. No scipy needed.
"""
import csv, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
CSV  = os.path.join(HERE, "..", "..", "01_CURRENT__rebuild", "runs",
                    "run19_contrast_checks", "contrast_per_encounter.csv")

def median(xs):
    xs = sorted(float(x) for x in xs)
    n = len(xs)
    if n == 0: return float("nan")
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m-1] + xs[m])

def mannwhitney_two_sided(a, b):
    """U for sample a vs b, normal approx with tie correction, two-sided erfc p."""
    a = [float(x) for x in a]; b = [float(x) for x in b]
    n1, n2 = len(a), len(b)
    combined = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    # average ranks (1-based), handling ties
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j + 1 < len(combined) and combined[j+1][0] == combined[i][0]:
            j += 1
        avg = 0.5 * ((i + 1) + (j + 1))
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    R1 = sum(r for r, (_, g) in zip(ranks, combined) if g == 0)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    U2 = n1 * n2 - U1
    U  = min(U1, U2)
    mu = n1 * n2 / 2.0
    # tie correction
    from collections import Counter
    cnt = Counter(v for v, _ in combined)
    N = n1 + n2
    tie = sum(t**3 - t for t in cnt.values())
    var = (n1 * n2 / 12.0) * ((N + 1) - tie / (N * (N - 1)))
    z = (U - mu) / math.sqrt(var) if var > 0 else 0.0
    p = math.erfc(abs(z) / math.sqrt(2))
    return U1, z, p

rows = list(csv.DictReader(open(CSV)))
EB   = lambda r: float(r["EB_mem"])
DN   = lambda r: float(r["Dn_mem"])

clean  = [r for r in rows if r["spec"] == "SHEATH_CONSISTENT"]                 # 60
rej    = [r for r in rows if r["spec"] in ("HOT_BOUNDARY_FLAG", "SHAPE_FLAG")] # 32
hot    = [r for r in rows if r["spec"] == "HOT_BOUNDARY_FLAG"]                 # 28
nondep = [r for r in rows if DN(r) > 1.0]                                      # 293

print("RUN 22 - supporting statistics (field-generic + tier-2 systematic)")
print(f"contributing N = {len(rows)}\n")

print("=== FACET 1: the field enhancement E_B is GENERIC (not a PDL discriminator) ===")
print(f"  clean 60   (spec=SHEATH_CONSISTENT): median E_B = {median(EB(r) for r in clean):.3f}")
print(f"  rejected32 (HOT_BOUNDARY+SHAPE_FLAG): median E_B = {median(EB(r) for r in rej):.3f}")
print(f"  hot 28     (HOT_BOUNDARY_FLAG)      : median E_B = {median(EB(r) for r in hot):.3f}"
      f"  (>1 in {sum(EB(r) > 1 for r in hot)}/{len(hot)} = {sum(EB(r) > 1 for r in hot)/len(hot):.0%})")
U1, z, p = mannwhitney_two_sided([EB(r) for r in clean], [EB(r) for r in rej])
print(f"  Mann-Whitney E_B clean-60 vs rejected-32: U={U1:.1f}, z={z:.3f}, two-sided p = {p:.3f}")
print(f"    => indistinguishable (p>0.05), and the median trends HIGHER in the rejected")
print(f"       contamination -> the field rise does NOT discriminate PDL from contamination.")
print(f"  non-depleted encounters (D_n>1, i.e. NO density depletion): N = {len(nondep)}")
print(f"    median E_B = {median(EB(r) for r in nondep):.3f}; "
      f"E_B>1 in {sum(EB(r) > 1 for r in nondep)}/{len(nondep)} = "
      f"{sum(EB(r) > 1 for r in nondep)/len(nondep):.0%}")
print(f"    => the field is enhanced even where there is no depletion at all.")
print(f"  (cf. RUN2 SHEATH_GEOM profile: E_B already 1.67 at mid-sheath s in [0.20,0.30]")
print(f"       where D_n=1.04 - the field rises smoothly across the WHOLE sheath: generic")
print(f"       compression/draping, not a localised depletion-layer feature.)")
print(f"  PARALLEL to the anticorrelation result (r_nB clean -0.81 vs rejected -0.80,")
print(f"  MWU p=0.81): BOTH legs of the n-|B| coupling are generic. Calling the field the")
print(f"  'robust half of the [depletion-layer] signature' headlines a non-PDL-specific quantity.\n")

print("=== FACET 2: tier-2 'statistically below unity' omits the dominant systematic ===")
dn_med = median(DN(r) for r in rows)
ci_lo, ci_hi = 0.871, 0.972            # RUN21 bootstrap median 95% CI (sampling, fixed boundary)
ci_half = 0.5 * (ci_hi - ci_lo)
print(f"  population median D_n = {dn_med:.4f}")
print(f"  SAMPLING uncertainty (RUN21, fixed boundary): Wilcoxon-vs-1 two-sided p=0.013;")
print(f"    bootstrap median 95% CI [{ci_lo}, {ci_hi}] -> half-width = {ci_half:.4f}")
print(f"  SYSTEMATIC uncertainty (RUN17 / Appendix C, correlated boundary offset):")
print(f"    magnetopause -0.5 R_E -> median D_n = 0.550   (excursion {0.550-dn_med:+.3f})")
print(f"    magnetopause +0.5 R_E -> median D_n = 1.135   (excursion {1.135-dn_med:+.3f}; >1 = NO decrease)")
print(f"    bow shock   +10%      -> median D_n = 1.452")
print(f"  systematic excursion (to +0.5 R_E) / sampling CI half-width = "
      f"{abs(1.135-dn_med)/ci_half:.1f}x")
print(f"  => the p=0.013 / CI[0.87,0.97] are conditional on an UNBIASED boundary; a plausible")
print(f"     correlated offset (Staples 2020: model under-tracks OUTWARD at high D_p) moves the")
print(f"     median to 1.135 (>1), erasing the sub-unity verdict. The systematic is ~{abs(1.135-dn_med)/ci_half:.0f}x the")
print(f"     sampling CI and is DIRECTIONAL. The fine-shell profile (S5) shares the same")
print(f"     coordinate, so it does not provide an independent escape.")
