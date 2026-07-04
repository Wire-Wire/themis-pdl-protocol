# Legacy / repair-era scripts (provenance only)

The original scripts of the fixed-axis -> radial recomputation that produced
`data/derived/archive_radial_catalogue_v2.csv` — the dual-geometry catalogue from which
`scripts/reproduce/run23_paired_provenance.py` re-derives the paired headline
(D_n 0.689 -> 1.114). Kept unmodified for provenance; they use the repair era's original
working vocabulary and path strings (mapped by the local `config.py`).

| Script | Role in the lineage |
|---|---|
| `recompute_subset_jbs.py` | Validated the Jelínek-2012 bow-shock radial coordinate on a cached subset before the archive run |
| `recompute_archive_radial.py` | The full-archive dual-geometry recompute (fixed-axis AND radial per encounter) |
| `build_radial_catalogues.py` | Converted the recompute output into the catalogue schema |
| `aggregate_radial.py` | Aggregated the catalogue into the headline-level numbers |
| `paired_v2.py` | The definitive paired comparison from the v2 catalogue |
| `config.py`, `psub.py`, `radial_models.py` | Local loaders so the above run self-contained |

The maintained, documented equivalents live in `scripts/reproduce/` (see `docs/SCRIPT_REFERENCE.md`).
