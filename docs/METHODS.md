# Methods — what is computed, and why

Distilled from the dissertation (§2–§3, §7–§8); every threshold below is tested in the committed sensitivity runs.

## 1. The normalised magnetosheath coordinate

Each sample gets `s = d_MP / (d_MP + d_BS)` — fractional position between the model magnetopause (s=0) and bow shock (s=1). **Why:** the sheath is thin and its boundaries move with the solar wind, so absolute position is meaningless across encounters [Pi24, MdW24]. **Decisive choice:** distances are measured **radially** through the spacecraft, not along the Sun–Earth axis — the fixed-axis chord misassigns off-subsolar samples and manufactures an apparent depletion (the paired test, `run23`: D_n 0.689→1.114, p=4×10⁻⁴⁴).

- **Magnetopause:** Shue et al. 1998 [Shue98]: `r(θ) = r0 (2/(1+cosθ))^α`, with `r0(Dp, Bz)` and `α(Dp, Bz)` from encounter-averaged OMNI. **Why Shue:** the community-standard dayside MP with explicit Dp/Bz response.
- **Bow shock:** Jelínek et al. 2012 [Jel12]: `r_BS(θ) = 2 R0 Dp^(−1/ε) / (cosθ + sqrt(cos²θ + λ² sin²θ))`, R0=15.02 R_E, ε=6.55, λ=1.17. **Why Jelínek:** fitted and validated on THEMIS crossings of the same era, and an *independent* surface (not a scaled copy of the MP), so the coordinate can respond to genuine sheath-thickness changes. It is dynamic-pressure-only by construction (weak Mach sensitivity for M>4).
- Boundary parameters are evaluated **once per encounter** (encounter-averaged OMNI); within-window boundary motion is unmodelled and is quantified in `run27` (median ΔMP 1.45 R_E → motivates the statistical-smearing reading).

## 2. Contrast metrics

Per encounter: background = median over s∈[0.6, 1.0] (outer sheath; insensitive choice, run3/App A); near-shell = s∈[0.05, 0.20). `D_n = n_near/n_bg`, `E_B = B_near/B_bg`. **Why within-encounter ratios:** each encounter is its own control, removing order-of-magnitude upstream variability. **Why s<0.05 is excluded:** boundary mispositioning concentrates magnetospheric contamination there.

## 3. Magnetosheath-membership screen

A near-shell sample counts only if: density > 0.3 cm⁻³, temperature within ×3 of the encounter's background T, flow > 0.2 × background flow. **Why these:** the magnetosphere/boundary layer is hot, tenuous and slow; shocked sheath is dense, sheath-temperature and flowing. Thresholds are deliberately permissive and tested in `run25`: sweeping the density floor 0.1→2.0 cm⁻³ moves D_n only 0.92→0.97 (floor non-binding below 1.0).

## 4. Per-event ion-spectral validation

Moments can be fooled by mixed populations; the ESA reduced-ion energy spectrum (32 log-spaced channels) is not. Three metrics per candidate (near-shell vs background time-median spectra): **peak-energy ratio** (sheath-like in [0.4, 2.5]; magnetospheric contamination shifts the peak to keV → ratio ≫1), **log-shape correlation** (≥0.8 sheath-like), **flux ratio** (the depletion itself). Result: 60 of 107 supported / 28 hot magnetospheric–boundary-layer / 4 shape-mismatch / 15 borderline; the populations are separated by an order of magnitude (peak ratio 1.0 vs ~24), so the classification is threshold-insensitive. **Why this carries the inference:** both the field enhancement and the n–|B| anticorrelation are generic near-boundary signatures (indistinguishable between clean and rejected candidates, `run21/run22`); only the spectra separate contamination — confirmed against the in-situ-to-OMNI moment-ratio alternative in `run29`.

## 5. Statistics

Paired two-sided Wilcoxon signed-rank for the coordinate test; Mann–Whitney U for group comparisons; bootstrap 95% CIs on medians (encounter = independent unit); log-OLS with VIF for confounder control (`run20`); forward synthetic injection for identifiability (`run18`: the Dp ordering is degenerate with a Dp-correlated boundary error ≤0.6 R_E, hence never claimed as a driver). Selection discipline: axes used to select candidates are never compared candidate-vs-baseline as findings (`run28` clean/circular split).

## 6. Selection-preserving environmental null and stability (runs 28–31)

"Under what conditions are candidates found?" is answerable only against what the candidate selection *itself* can induce. **Design** (`run30`): the density component of the candidate screen (0.40 ≤ D_n ≤ 0.90, n_near ≥ 2 cm⁻³) is replayed on 5,000 permutations of the near-shell densities across the northward universe, within strata (D_p × cone tertiles, and D_p tertiles alone as the worst-case drag channel); the observed candidate-vs-baseline gap on each axis is then compared with the null band. **Why:** a density-selecting gate can mechanically induce a preference on any axis correlated with background density — the null quantifies that reach. **Result:** the dense-background gap sits *inside* the null (reported as a property of the candidate definition), while the quasi-perpendicular cone-angle gap *exceeds* it (exceedance 0.007–0.024). **Stability** (`run31`): the 60 supported candidates span 5 probes and 14 years; the cone preference survives every leave-one-probe-out and leave-one-year-out; a re-fetched within-window IMF check shows the preference is not an OMNI-direction artefact and *strengthens* on upstream-steady subsets (+11.2° at within-window cone IQR ≤ 15°, p = 0.020); and because the archive is near-subsolar (SZA < 30°), the analytic Jelínek shock-normal tilt (≤ 18°) makes cone angle ≈ θ_Bn without a streamline reconstruction. Discipline throughout: these are **conditions of candidate discovery in a selection-limited sample**, never occurrence rates or demonstrated drivers.

## References

| Key | Reference | Used for |
|---|---|---|
| Shue98 | Shue et al., JGR 103, 17691 (1998). doi:10.1029/98JA01103 | magnetopause model |
| Jel12 | Jelínek, Němeček & Šafránková, JGR 117, A05208 (2012). doi:10.1029/2011JA017252 | bow-shock model |
| ZW76 | Zwan & Wolf, JGR 81, 1636 (1976). doi:10.1029/JA081i010p01636 | PDL theory (pile-up + drainage) |
| Cro79 | Crooker et al., JGR 84, 869 (1979). doi:10.1029/JA084iA03p00869 | early PDL observations |
| Phan94 | Phan et al., JGR 99, 121 (1994). doi:10.1029/93JA02444 | low-shear PDL climatology |
| Sou08 | Soucek, Lucek & Dandouras, JGR 113, A04203 (2008). doi:10.1029/2007JA012649 | mirror-mode n–B anticorrelation (non-uniqueness) |
| SK92 | Southwood & Kivelson, JGR 97, 2873 (1992). doi:10.1029/91JA02446 | slow-mode/pressure-balance context |
| Sta20 | Staples et al., JGR 125, e2019JA027289 (2020). doi:10.1029/2019JA027289 | MP-model error under compression (Dp-correlated error physical basis) |
| Pi24 | Pi et al., Front. Astron. Space Sci. 11, 1401078 (2024). doi:10.3389/fspas.2024.1401078 | normalised-coordinate sheath profiles (THEMIS) |
| MdW24 | Michotte de Welle et al., Front. Astron. Space Sci. 11, 1427791 (2024). doi:10.3389/fspas.2024.1427791 | sheath n/B vs IMF orientation (cone-angle dependence) |
| Ang08 | Angelopoulos, Space Sci. Rev. 141, 5 (2008). doi:10.1007/s11214-008-9336-1 | THEMIS mission |
| Aus08 | Auster et al., Space Sci. Rev. 141, 235 (2008). doi:10.1007/s11214-008-9365-9 | FGM |
| McF08 | McFadden et al., Space Sci. Rev. 141, 277 (2008). doi:10.1007/s11214-008-9440-2 | ESA |
| KP05 | King & Papitashvili, JGR 110, A02104 (2005). doi:10.1029/2004JA010649 | OMNI |
