# Pipeline map — which script produces which number

> Scripts live in `scripts/reproduce/` (significance/robustness runs in `statistics/`), figures in `scripts/figures/`; the reusable core is `src/pdl_protocol/`. Fuller per-script detail (reads/writes): [SCRIPT_REFERENCE.md](SCRIPT_REFERENCE.md). Dissertation anchors use the **final** section/figure numbering; repo figure *filenames* keep the pre-renumbering scheme (map in the README).

Chain order (each stage reads the frozen substrate and/or earlier committed outputs; nothing is recomputed silently):

| # | Script | What it does | Committed output | Dissertation anchor |
|---|---|---|---|---|
| 0 | `build_substrate.py` | Builds 6,248 per-encounter NPZ time series (6-h window per probe per dayside day; n, \|B\|, β, p_th, v + per-encounter OMNI drivers and Shue/Jelínek parameters) from CDAWeb | `data/substrate/` (rebuilt locally) | §2.3 |
| 1 | `psub.py`, `radial_models.py`, `run_hardening_checks.py` | Shared loaders: threaded substrate mapper; Jelínek bow shock; membership screen + shell contrasts | — | §3.1–3.4 |
| 2 | `run2_profile_cube.py` (+`run2b`, `run2c`) | Radial fine-shell profile of D_n/E_B/β vs s, with systematic + random boundary-error bands | `run2_profile_cube/` | §5, Fig 3; App C |
| 3 | `run10_selection_funnel.py` | The observability funnel 6,248 → 3,869 → 1,187 → 661 → 332 → 107 → 60 | `run10_selection/` | §6, Table 2, Fig 9 |
| 4 | `run8a_select_candidates.py` → `run9a_classify_candidates.py` | Conservative absolute moment screen on northward encounters → 107 candidates (legacy "confirmed" labels = "spectrally supported" in the dissertation) | `run8_atlas/`, `run9_candidates/` | §6.2 |
| 5 | `run12_spectral_check_all.py` → `run13_validation_status.py` | Per-event ESA ion-spectral validation of all 107 → 60/28/4/15 tiers | `run12_spectral/`, `run13_validation/` | §3.5, §6.3, Fig 5, App B |
| 6 | `run16_dp_ordering.py` | The D_p tertile ordering of the contrast | `run16_dp_ordering/` | §8.1, Table 4, Fig 8 |
| 7 | `run17_robustness.py` | Boundary-placement ±0.5 R_E / BS ±10% / foreshock-exclusion robustness | `run17_robustness/` | App C |
| 8 | `run18_synthetic_injection.py` | **Forward test**: the D_p ordering is degenerate with a D_p-correlated boundary error ≤0.6 R_E | `run18_synthetic/` | §8.4, App D |
| 9 | `run19_contrast_validation.py` (+`run19b`), `run20_regression.py`, `run21_contrast_stats.py`, `run22_field_specificity.py` | Pressure budget + n–\|B\| coupling; the Table-5 regression; contrast statistics; field-genericity tests | `run19_contrast_checks/` … `run22_field_generic/` | §6.4, §9.1, Table 5, App E |
| 10 | `run23_paired_provenance.py` | **The headline**: the 672-row paired fixed-axis-vs-radial table | `run23_paired/paired_1d_vs_radial.csv` | §4, Fig 2 |
| 11 | `run24_sza_reconciliation.py` | Paired shift vs solar-zenith angle (the artefact is a broad-bin effect, largest near subsolar) | `run24_sza_reconciliation/` | §4 |
| 12 | `run25_threshold_sensitivity.py`, `run26_logratio.py` | Membership-threshold sweep; log-ratio distribution | `run25_…/`, `run26_…/` | §3.4, App F |
| 13 | `fetch_omni_contributing.py` | 1-min OMNI re-fetch for the contributing windows | `data/omni_cache/` (local) | feeds 14/16/18 |
| 14 | `run27_boundary_motion_qc.py` | Within-window boundary-motion QC (median ΔMP 1.45 R_E; region-swap exclusion sweep) | `run27_boundary_motion_qc/` | §9.2, App C |
| 15 | `run28_candidate_conditions.py`, `make_candidate_atlas.py` | Candidate-vs-baseline characterisation (clean/circular axis discipline) + top-6 atlas | `run28_candidate_conditions/` | §7.1–7.4, App G |
| 16 | `run29_themis_omni_ratio.py` + `run29b_ratio_avnn.py` | THEMIS-to-OMNI ratio check: background = textbook shocked sheath (3.77 / 0.246); the ratio classifier inherits moment blindness; A-vs-baseline addendum +0.83 [+0.40, +1.18] | `run29_themis_omni_ratio/` | §7.2, §7.5, App G |
| 17 | `run30_selection_null.py` | **The denominator-selection permutation null** (5,000 permutations × 2 strata variants): the dense-background gap sits inside the null; the cone gap exceeds it (P = 0.007–0.024) | `run30_selection_null/` | §7.3–7.4, App G |
| 18 | `run31ac_repr_thetabn.py`, `run31b_cone_stability.py` (+`run31b_analysis.py`) | Representativeness (5 probes / 14 years, leave-one-out), OMNI-stable cone re-check (+11.2° at IQR ≤ 15°), analytic θ_Bn≈cone bound | `run31_candidate_context/` | §7.1, §7.4, App G |
| F | `make_figs_1_6.py`, `make_fig7.py`, `make_fig8_deepdive.py`, `make_fig9_environment.py` | All nine dissertation figures from committed outputs (+ demo/full data where noted in the tutorial) | `figures/` | Figs 1–9 |
| V | `verification/*.py` | Independent re-implementations and adversarial checks of the coordinate-artefact result | — | §4 verification |

**Counting units:** every count is an *encounter* (one 6-h window per probe per day; multiple back-and-forth magnetopause crossings inside a window = one encounter). 672 = evaluable under both coordinate constructions (§4 paired test); 661 = contributing (adds the background + membership-shell requirements); 667 = threshold-sweep evaluable set (App F note); 332 = northward validation subset; 107/60 = moment-classified / spectrally supported candidates.
