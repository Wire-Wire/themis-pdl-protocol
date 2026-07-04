# Changelog

## v1.2.0 — 2026-07-03

The two-door release: the repository now separates its two audiences explicitly.

- **New: one-command research tools (Door B).** `pdl-find` (= `pdl_protocol.find_candidates`): a directory of encounter NPZ files -> a candidate table (radial coordinate -> membership screen -> contrast -> the dissertation's Sec.-6.2 conservative criteria -> per-event spectral tier where a spectrum is available), with the selection-limited scope note built into the output. `pdl-contrast` (= `shell_profile` / `population_contrast` / `paired_coordinate_check`): population shell profile, near-shell statistics, and the Sec.-4 fixed-axis-vs-radial artefact check runnable on any dataset. New `examples/04_find_candidates.py`, `examples/05_contrast_profile.py`; console scripts registered via pyproject.
- **Two-door documentation.** README now routes readers at the top (Door A: verify + `docs/REPRODUCING.md` + PIPELINE_MAP; Door B: install + `pdl-find` + TUTORIAL); the layout table tags every path by door and documents the `outputs/` working tree vs `committed_outputs/` reference-tree design. TUTORIAL rewritten as the pure researcher path (real captured outputs at every level); new `docs/REPRODUCING.md` carries the audit path (verify -> seed -> re-run/diff -> figures -> full rebuild).
- **Removals** (serving neither audience): `run5_hq` (script + committed output — compared against a predecessor report's subset the dissertation never mentions), the superseded pre-fig9 forest PNGs, the redundant library copy of `run_hardening_checks.py` (API lives in `core.py`; the runnable copies remain in `scripts/`); `logratio_histogram.png` relocated to its run's committed directory.
- `archive/legacy_repair/` slimmed to the catalogue lineage (5 scripts + loaders + README); six superseded variants (early paired comparison, previews, subset sensitivity runs, a loader duplicate) removed.
- run4's reference entry corrected: its per-shell membership pass-rates are the source of the dissertation's §6.1 90–96% figures.

## v1.1.0 — 2026-07-03

Synchronised with the final dissertation text and made properly installable.

- **Backfilled the post-snapshot science** (the v1.0 snapshot predated these by hours to days):
  - `run30_selection_null` — the denominator-selection permutation null (5,000 permutations, two strata variants) that absorbs the dense-background preference and that the quasi-perpendicular cone preference exceeds (dissertation §7.3–7.4, Appendix G);
  - `run31_candidate_context` — representativeness of the 60 supported candidates (5 probes, 14 years, leave-one-out), the OMNI-stable cone re-check (+11.2° at within-window cone IQR ≤ 15°, p = 0.020) and the analytic θ_Bn≈cone bound (§7.1, §7.4, Appendix G);
  - the `run29` **A-vs-NN addendum** (background compression-ratio excess +0.83 [+0.40, +1.18], p = 3.9×10⁻⁵) appended to the committed output;
  - `fig9_environment` (dissertation **Figure 7**) — figure + generation script.
- **Packaging**: `pyproject.toml` (`pip install -e .`, extras `plot` / `fetch`); library API extracted into `pdl_protocol.core` (no path shims); examples/tests/tools fall back to the source tree only when the package is not installed.
- **Verification**: `verify_key_numbers.py` extended 15 → 22 anchored checks (runs 27, 29-addendum, 30, 31); the full suite (verifier, tests, all examples, quickstart, figure scripts 1–6/8/9) was executed in a clean venv before release.
- **Bug fixes found by that run**: `config.ROOT` pointed at `src/` instead of the repository root, so `data/substrate` never resolved on a clean checkout (v1.0's examples/tests could not actually run as shipped); all figure scripts now write into `outputs/figures/` rather than beside themselves; `examples/02` now states honestly that its simplified shell windows give flux ratio 0.46 vs the pipeline's exact 0.39 (identical classification).
- **Docs**: README rewritten around the dissertation's two principal results; new `docs/TUTORIAL.md`, `docs/SCRIPT_REFERENCE.md` (every script, one line each), `docs/API.md`; METHODS/DATA/PIPELINE_MAP updated to the final dissertation section/figure numbering (with the repo-filename ↔ dissertation-figure map).
- **Metadata**: CITATION.cff points at the final dissertation title; CI tests Python 3.10/3.12, installs the package, and runs all examples + the quickstart figure.
- Committed TXT files re-encoded to UTF-8 (16 were GBK and rendered as mojibake off the original machine); internal review-stage tags removed from script docstrings and TXT narrative headers (numbers untouched); the run19 summary file renamed `CONTRAST_VALIDATION.txt` and its writer aligned; supplementary figure copies now go to `outputs/supplementary/`.
- Terminology aligned with the dissertation's published vocabulary (e.g. "spectrally supported", "hot magnetospheric / boundary-layer contaminant"); committed outputs keep their legacy internal labels unmodified, as documented in the README terminology note.

## v1.0.0 — 2026-06-11

Initial public snapshot: importable protocol core, three worked examples, reproduction chain
(stages 0–16), committed outputs for runs 2–29, 3-encounter demo substrate, 15-check verifier, CI.
