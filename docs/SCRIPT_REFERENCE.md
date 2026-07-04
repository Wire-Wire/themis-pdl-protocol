# Script reference — every script, what it does, reads and writes

Paths are repo-relative. "outputs/…" means the git-ignored working tree (seed it with `python tools/seed_outputs.py`); the committed reference copies live in `committed_outputs/`. Dissertation anchors use the **final** section/figure numbering (see the figure map in the README for repo-filename ↔ dissertation-figure).

## Root

| Script | What it does | Reads | Writes |
|---|---|---|---|
| `verify_key_numbers.py` | Re-derives all 22 headline anchors (paired test, funnel, population contrast, tiers, environment nulls, stability) and exits non-zero on any mismatch — the CI gate | `committed_outputs/` | — |
| `quickstart_fig2.py` | Draws the paired coordinate comparison (dissertation Fig 2b–c) in seconds | `committed_outputs/run23_paired/` | `outputs/quickstart_fig2.png` |

## examples/

| Script | What it does | Dissertation anchor |
|---|---|---|
| `01_one_encounter.py` | Full per-encounter protocol (coordinate → screen → contrast) on a shipped real encounter; prints D_n = 0.554, E_B = 2.279 | §3.1–3.4 |
| `02_classify_spectrum.py` | Spectral metrics + classification of the shipped deep-dive spectrogram (1.00 / 0.96 / 0.39 → sheath-like) | §3.5 |
| `03_bring_your_own_data.py` | The same loop on a synthetic encounter with a known injected depletion (recovers 0.5 / 2.0) — the template for your own data | docs/DATA.md schema |
| `04_find_candidates.py` | The one-command finder on the demo encounters (equivalent CLI: `pdl-find`) | §6.2 criteria |
| `05_contrast_profile.py` | Shell profile + population contrast + paired coordinate check (equivalent CLI: `pdl-contrast --paired-check`) | §5, §9.1, §4 |

## tools/

| Script | What it does |
|---|---|
| `get_data.py` | Data acquisition: `--check` inventory · `--event EID` (1-min OMNI + ESA spectrogram for one encounter) · `--omni` (1-min OMNI for every substrate encounter; feeds runs 27/29/31) · `--full-substrate` (full-rebuild guidance) |
| `seed_outputs.py` | Copies `committed_outputs/` into `outputs/` so any pipeline stage can run standalone |

## src/pdl_protocol/ (the installable library)

| Module | Contents |
|---|---|
| `core.py` | `load_encounter` (NPZ → radial s + background state), `member_mask` (the §3.4 screen), `shell_contrast` (D_n, E_B) — the exact pipeline operations, docs in [API.md](API.md) |
| `coords.py` | Shue-1998 magnetopause (`shue_r0/alpha/r`), Jelínek-2012 bow shock (`jelinek_r`), the radial normalised coordinate (`compute_s`) |
| `radial_models.py` | The Jelínek-2012 bow-shock surface used by the whole chain (constants + `jelinek_bs_r`, `compute_s_radial_jbs`) |
| `spectral.py` | `spectral_metrics` (peak-energy ratio, log-shape correlation, flux ratio) + `classify_spectrum` (§3.5 thresholds) |
| `psub.py` | `pmap` — threaded map over every substrate NPZ (np.load + numpy release the GIL); `list_files` |
| `config.py` | Repo path layout, `PDL_ROOT` override, and `P()` — maps the original analysis-machine paths kept in the reproduction scripts onto this layout |
| `finder.py` | **`pdl-find`** / `find_candidates()`: the whole protocol over a directory of encounters -> candidate table with spectral tiers (criteria of dissertation Sec. 6.2; scope note included) |
| `contrast.py` | **`pdl-contrast`** / `shell_profile()`, `population_contrast()`, `paired_coordinate_check()`: population profile, near-shell statistics, and the Sec.-4 fixed-axis-vs-radial artefact check |

## scripts/reproduce/ — the dissertation chain, stage order

| # | Script | What it does | Committed output | Anchor |
|---|---|---|---|---|
| 0 | `build_substrate.py` | Builds the 6,248 per-encounter NPZ time series from CDAWeb (6-h window per probe per dayside day, centred on the day's minimum solar-zenith-angle sample; SZA<30°, r=8–25 R_E) | `data/substrate/` (local) | §2.3 |
| 0b | `fetch_omni_contributing.py` | 1-min OMNI re-fetch for the contributing encounters | `data/omni_cache/` (local) | runs 27/29/31 |
| 1 | `psub.py`, `radial_models.py`, `run_hardening_checks.py`, `config.py` | Shared loaders (per-directory copies so each stage folder is self-contained) | — | §3 |
| 2 | `run2_profile_cube.py` | Radial fine-shell profile of D_n/E_B/β vs s | `run2_profile_cube/` | §5, Fig 3 |
| 3 | `run10_selection_funnel.py` | The observability funnel 6,248 → … → 661 contributing → 332 northward | `run10_selection/` | §6, Table 2, Fig 9 |
| 4 | `run7_membership_screen.py` | The temperature–flow–density membership screen + per-shell pass rates | `run7_membership/` | §3.4, §6.1 |
| 5 | `run8a_select_candidates.py` → `run9a_classify_candidates.py` | Conservative absolute moment screen on northward encounters → the 107 candidates (legacy "confirmed" labels = "spectrally supported" in the dissertation) | `run8_atlas/`, `run9_candidates/` | §6.2 |
| 6 | `run12_spectral_check_all.py` → `run13_validation_status.py` | Per-event ESA ion-spectral validation of all 107 → 60/28/4/15 tiers | `run12_spectral/`, `run13_validation/` | §3.5, §6.3, Fig 5, App B |
| 7 | `run23_paired_provenance.py` | **The headline**: the 672-row paired fixed-axis-vs-radial table | `run23_paired/` | §4, Fig 2 |
| 8 | `run_selection_function.py` | Which encounters reach which s-shells (observability accounting) | `run2_profile_cube/SELECTION_FUNCTION.txt` | §9.4 |

## scripts/reproduce/statistics/ — robustness & significance

| Script | What it does | Committed output | Anchor |
|---|---|---|---|
| `run2b_robustness.py`, `run2c_random_mc.py` | Boundary-standoff systematic sweep + per-encounter random model-error Monte Carlo | `run2_profile_cube/RUN2b/2c` | App C |
| `run3_bin_multiverse.py` | Bin/background/contamination multiverse (guards against cherry-picking) | `run3_bin_multiverse/` | App A |
| `run4_favourable_test.py` | Pre-registered favourable-regime test; its per-shell membership pass-rates are the source of the 90–96% figures | `run4_favourable/` | §6.1 |
| `run11_bz_confounder.py` | IMF-B_z split on the full membership-screened set | `run11_bz_confounder/` | §8.3 |
| `run14_mechanism_probe.py` | Pre-registered mechanism probe (T, β structure of the near shell) | `run14_mechanism/` | §6.4 context |
| `run15_subset_profile.py` | Does the contrast survive in the spectrally supported subset alone | `run15_subset_profile/` | §6.4 |
| `run16_dp_ordering.py` | The dynamic-pressure tertile ordering | `run16_dp_ordering/` | §8.1, Table 4, Fig 8 |
| `run17_robustness.py` | Boundary ±0.5 R_E / bow-shock ±10% / foreshock-exclusion stress test | `run17_robustness/` | App C |
| `run18_synthetic_injection.py` | **Forward test**: Dp-independent synthetic layers produce no ordering; a Dp-correlated boundary error ≤0.6 R_E reproduces all of it | `run18_synthetic/` | §8.4, App D |
| `run19_contrast_validation.py`, `run19b_pressure_supp.py` | Pressure-budget + n–\|B\|-coupling checks behind Appendix E | `run19_contrast_checks/` | §6.4, App E |
| `run20_regression.py` | The Table-5 log-OLS confounder regression, committed for provenance | `run20_regression/` | §8.3, Table 5 |
| `run21_contrast_stats.py` | Population-contrast supporting statistics (Wilcoxon-vs-1, bootstrap CI, r_nB group test) | `run21_contrast_stats/` | §9.1 |
| `run22_field_specificity.py` | Field enhancement is generic: clean-vs-rejected MWU, enhancement in non-depleted encounters | `run22_field_generic/` | §9.1 |
| `run24_sza_reconciliation.py` | Paired shift vs solar-zenith angle (refutes the off-subsolar reading of the artefact) | `run24_sza_reconciliation/` | §4 |
| `run25_threshold_sensitivity.py`, `run26_logratio.py` | Membership-threshold sweep; log-ratio presentation | `run25_…/`, `run26_…/` | §3.4, App F |
| `run27_boundary_motion_qc.py` | Within-window model-boundary motion QC (median ΔMP 1.45 R_E; region-swap exclusion) | `run27_boundary_motion_qc/` | §9.2, App C |
| `run28_candidate_conditions.py` | Candidate-vs-baseline characterisation with the clean/circular axis discipline | `run28_candidate_conditions/` | §7.1–7.4, App G |
| `run29_themis_omni_ratio.py`, `run29b_ratio_avnn.py` | THEMIS-to-OMNI ratio check (background validation; the classifier inherits moment blindness) + the A-vs-baseline addendum (+0.83 [+0.40, +1.18]) | `run29_themis_omni_ratio/` | §7.2, §7.5, App G |
| `run30_selection_null.py` | **The denominator-selection permutation null** (5,000 permutations, two strata variants): what preference can the density gate alone induce? Dense-background gap sits inside the null; the cone gap exceeds it | `run30_selection_null/` | §7.3–7.4, App G |
| `run31ac_repr_thetabn.py` | Representativeness of the 60 (probes/years/Herfindahl), leave-one-out stability of the headline axes, analytic θ_Bn≈cone bound | `run31_candidate_context/` | §7.1, §7.4, App G |
| `run31b_cone_stability.py`, `run31b_analysis.py` | Re-fetched within-window IMF: cone steadiness A-vs-baseline, and the preference on upstream-steady subsets (+11.2° at IQR ≤ 15°) | `run31_candidate_context/` | §7.4, App G |
| `run_hardening_checks.py` | Membership-threshold sweep + mirror-mode aliasing check (detrended residual test) | `run7_membership/HARDENING_CHECKS.txt` | §3.4, App E |

## scripts/figures/

| Script | Draws | Needs |
|---|---|---|
| `make_figs_1_6.py` | Dissertation Figures 1, 2, 3, 5, 8, 9 (repo files fig1–fig6) | seeded outputs |
| `make_fig7.py` | Dissertation Figure 4 — three event archetypes (repo fig7) | full substrate (three specific encounters) |
| `make_fig8_deepdive.py` | Dissertation Figure 6 — the 2020-12-12 deep dive with spectrogram + 1-min OMNI (repo fig8) | shipped demo data |
| `make_fig9_environment.py` | Dissertation Figure 7 — environment forest plot with the run30 null bands; self-checks against committed outputs (repo fig9) | seeded outputs |
| `make_candidate_atlas.py` | The 6 cleanest spectrally supported events on one sheet | full substrate |

## scripts/verification/ — independent checks of the coordinate result

| Script | What it does |
|---|---|
| `xcheck_independent.py` | From-scratch independent re-implementation of the fixed-axis and radial coordinates; reproduced both to the digit |
| `adversarial.py` | Adversarial geometry variant (blunt-magnetopause normal distances) — does the artefact survive hostile constructions |
| `robustness.py` | Interpolation/gap-handling robustness of the recompute |
| `bound.py` | Analytic bounds on the coordinate difference |
| `validate_C.py` | Validation checks for the dynamic-pressure non-identifiability result |

## archive/legacy_repair/ — provenance only

The five scripts of the original fixed-axis→radial recompute lineage that produced `data/derived/archive_radial_catalogue_v2.csv` (the input of `run23_paired_provenance.py`), plus their local loaders. Kept unmodified for provenance; superseded variants were removed in v1.2.0. Details: `archive/legacy_repair/README.md`.
