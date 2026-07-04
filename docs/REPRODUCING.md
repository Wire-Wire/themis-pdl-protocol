# Reproducing the dissertation (the reviewer path)

This is **Door A**: you are checking that the dissertation's numbers, figures and methods are real. You will not need to install anything beyond numpy/scipy (and matplotlib for figures); if you want to *use* the protocol on data instead, go to [TUTORIAL.md](TUTORIAL.md).

## 1. Verify every headline number (30 seconds, no downloads)

```bash
pip install -e .           # or just: pip install numpy scipy
python verify_key_numbers.py
```

Expected: 22 lines of `PASS` ending in `22/22 checks passed`. The script re-derives the paired coordinate test, the selection funnel, the population contrast, the spectral tiers, the environment nulls and the stability checks **from the committed CSV/TXT outputs alone**. If this fails, the repository copy is damaged; re-clone.

## 2. The two output trees (how reproduction is checked here)

| Tree | Tracked by git | Role |
|---|---|---|
| `committed_outputs/` | yes | the **frozen reference**: exactly what the dissertation used; read-only by convention |
| `outputs/` | no (git-ignored, created locally) | the **working tree**: every re-run writes here |

This split is the reproduction mechanism: re-run a stage, then diff `outputs/runN/` against `committed_outputs/runN/`. Scripts can never overwrite the reference, so "what the dissertation used" and "what you just produced" stay distinguishable.

```bash
python tools/seed_outputs.py     # copy the committed reference into outputs/ so any stage can run standalone
```

## 3. Re-run a stage and compare

Stages that only read earlier outputs run immediately after seeding, e.g. the Table-5 regression or the selection null:

```bash
python scripts/reproduce/statistics/run20_regression.py
python scripts/reproduce/statistics/run30_selection_null.py
diff outputs/run30_selection_null/RUN30_SELECTION_NULL.txt committed_outputs/run30_selection_null/RUN30_SELECTION_NULL.txt
```

Stages that re-read the full 6,248-encounter substrate (`run2`, `run7`, `run17`, `run25`, …) need the full rebuild (section 6). [PIPELINE_MAP.md](PIPELINE_MAP.md) lists every stage in chain order with its committed output and its dissertation anchor; [SCRIPT_REFERENCE.md](SCRIPT_REFERENCE.md) describes every script in one line.

## 4. Reproduce the figures

```bash
pip install -e ".[plot]"
python quickstart_fig2.py                        # the headline comparison, from the committed CSV
python scripts/figures/make_figs_1_6.py          # dissertation Figures 1, 2, 3, 5, 8, 9
python scripts/figures/make_fig8_deepdive.py     # dissertation Figure 6 (shipped demo data)
python scripts/figures/make_fig9_environment.py  # dissertation Figure 7 - self-checks its numbers
                                                 #   against committed run28/run29/run30 before writing
```

All are written into `outputs/figures/`; compare with the committed `figures/`. (`make_fig7.py`, dissertation Figure 4, needs three encounters that are not in the demo subset — full rebuild required. Repo figure *filenames* use the pre-renumbering scheme; the map is in the README.)

## 5. What the committed evidence covers

- `committed_outputs/run2…run31` — every dissertation number's source (the README's key-number table gives the headline→directory map).
- `data/derived/archive_radial_catalogue_v2.csv` — the dual-geometry catalogue behind the §4 paired headline; `scripts/reproduce/run23_paired_provenance.py` re-derives 0.689→1.114 from it.
- `scripts/verification/` — independent re-implementations of the coordinate result (a from-scratch cross-check reproduced both constructions to the digit) and adversarial geometry variants.
- `archive/legacy_repair/` — the original recompute scripts that produced the catalogue (kept unmodified for provenance).

## 6. Full rebuild from the public archive (hours, ~2 GB)

```bash
pip install -e ".[fetch]"
python scripts/reproduce/build_substrate.py   # THEMIS FGM/ESA + OMNI 2007-2025 via CDAWeb
python tools/get_data.py --omni               # 1-min OMNI per encounter (runs 27/29/31)
```

then the numbered chain in [PIPELINE_MAP.md](PIPELINE_MAP.md) (stages 2 → 18). Honest expectations: the builder fetches from the live CDAWeb archive; if the archive re-calibrates data between now and the original build (2026-05-29), per-encounter values can differ at the margin. The committed outputs are the frozen reference of record; that is why they are committed.

## 7. Methods and thresholds

[METHODS.md](METHODS.md): what is computed and why, every threshold with its sensitivity run, the selection-discipline rules (clean vs circular axes, the permutation null), and the 14 DOI'd references behind the model choices.
