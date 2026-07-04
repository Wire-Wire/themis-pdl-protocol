# API reference — `pdl_protocol`

```python
pip install -e .            # from the repo root
import pdl_protocol         # pdl_protocol.__version__ -> "1.2.0"
```

The library is deliberately small: the exact operations of the dissertation pipeline, nothing speculative. All array arguments are numpy arrays; units follow the NPZ schema in [DATA.md](DATA.md) (R_E, cm⁻³, nT, km/s, nPa, eV).

## Encounter-level protocol (`pdl_protocol.core`)

### `load_encounter(d) -> dict | None`

Turn one substrate NPZ (`d = np.load(path, allow_pickle=True)`) into an analysis-ready dict: computes the radial normalised coordinate `s` from the encounter's Shue parameters (`mp0`, `alpha`) and dynamic pressure (`dp`), keeps samples geometrically inside the model magnetosheath, and attaches the encounter's own background medians. Returns `None` when the encounter has no valid boundary solution or too few samples (the pipeline's silent-skip rule).

Returned keys: `s, n, b, T, v, beta` (arrays, geometric-sheath samples) and `bz, n_bg, b_bg, T_bg, v_bg` (scalars).

### `member_mask(e, TBAND=3.0, VFRAC=0.2, NFLOOR=0.3) -> bool array`

The §3.4 magnetosheath-membership screen on a loaded encounter: dense (`n > NFLOOR` cm⁻³), sheath-temperature (within a factor `TBAND` of the background temperature), flowing (`v > VFRAC × v_bg`). Defaults are the dissertation baseline; sensitivity in committed `run25` (Appendix F).

### `shell_contrast(e, mask, lo, hi) -> (D_n, E_B) | (None, None)`

Median near-shell/background ratios over `s ∈ [lo, hi)` for the samples passing `mask`. The dissertation's analysis shell is `(0.05, 0.20)`; background is fixed at `s ∈ [0.6, 1.0]` inside `load_encounter`.

```python
e = load_encounter(np.load("data/substrate/2020-12-12_tha.npz", allow_pickle=True))
dn, eb = shell_contrast(e, member_mask(e), 0.05, 0.20)   # -> 0.554, 2.279
```

## Boundary models & coordinate (`pdl_protocol.coords`)

| Function | Meaning |
|---|---|
| `shue_r0(dp, bz)` | Shue-1998 subsolar magnetopause standoff (R_E) from D_p (nPa) and IMF B_z (nT) |
| `shue_alpha(dp, bz)` | Shue-1998 flaring parameter |
| `shue_r(cos_theta, r0, alpha)` | Magnetopause radius at angle θ from the Sun–Earth line |
| `jelinek_r(cos_theta, dp)` | Jelínek-2012 bow-shock radius (R0 = 15.02 R_E, ε = 6.55, λ = 1.17) |
| `compute_s(x, y, z, mp0, alpha, dp)` | The radial normalised coordinate s = d_MP/(d_MP+d_BS); 0 at the magnetopause, 1 at the bow shock; NaN outside |

Positions are GSE in R_E. `compute_s` is the construction whose fixed-axis counterpart the dissertation shows to manufacture an artefact (§4); the fixed-axis variant exists only inside `paired_coordinate_check` as the explicitly labelled comparator.

## Spectral validation (`pdl_protocol.spectral`)

### `spectral_metrics(near_spec, bg_spec, energies) -> (peak_ratio, shape_corr, flux_ratio)`

Compare the time-median near-shell and background ion energy spectra (differential energy flux vs the 32 ESA channels): peak-energy ratio, log-shape Pearson correlation, integrated-flux ratio.

### `classify_spectrum(peak_ratio, shape_corr) -> "sheath-like" | "boundary" | "borderline"`

The §3.5 tiers: sheath-like when the peak stays (ratio ∈ [0.4, 2.5]) and the shape holds (corr ≥ 0.8); boundary when the peak jumps toward keV or the shape decorrelates. Thresholds are non-critical — the populations separate by an order of magnitude (peak ratio ≈ 1 vs ≈ 24).

## One-command tools (`pdl_protocol.finder`, `pdl_protocol.contrast`)

### `find_candidates(substrate_dir=None, out_csv=None) -> list[dict]`

The whole protocol over a directory of encounter NPZ files (default `data/substrate/`): per encounter, radial coordinate -> membership screen -> near-shell contrast -> the conservative candidate criteria of dissertation Sec. 6.2 -> spectral classification wherever `data/events/esa_<date>.npz` exists. Each dict: `eid, D_n, E_B, n_near, T_near_eV, beta_near, n_samples, bz_nT, candidate, status`. CLI: **`pdl-find [dir] [--csv out]`**. Scope: candidates are selection-limited indications, never confirmed layers.

### `shell_profile(substrate_dir=None, shells=...) -> list[dict]`

Population radial profile: per shell, the across-encounter median D_n / E_B / beta and contributing N (the Sec.-5 profile for your data).

### `population_contrast(substrate_dir=None, shell=(0.05, 0.20)) -> dict`

Near-shell population statistics: D_n / E_B medians with IQR and the depleted fraction (Sec. 9.1 form).

### `paired_coordinate_check(substrate_dir=None, near=(0.2, 0.4), bg=(0.6, 1.0)) -> dict`

The Sec.-4 paired design on your data: broad-bin D_n under the radial vs the fixed Sun-Earth-axis construction on the encounters evaluable under both, with the median paired shift and a Wilcoxon test. A large shift means the coordinate choice, not the plasma, produces the apparent depletion. CLI: **`pdl-contrast [dir] [--paired-check] [--csv out]`**.

## Bulk processing (`pdl_protocol.psub`)

| Function | Meaning |
|---|---|
| `pmap(fn, files=None, workers=…, with_name=False)` | Apply `fn(npz)` (or `fn(basename, npz)`) to every NPZ in `data/substrate/` in threads; returns non-None results. Order not guaranteed — aggregate order-independently. |
| `list_files()` | The substrate NPZ paths |

```python
from pdl_protocol import pmap, load_encounter, member_mask, shell_contrast

def one(d):
    e = load_encounter(d)
    if e is None: return None
    dn, eb = shell_contrast(e, member_mask(e), 0.05, 0.20)
    return (dn, eb) if dn is not None else None

results = pmap(one)      # the population loop of the dissertation, in four lines
```

## Paths (`pdl_protocol.config`)

`ROOT` (repo root, override with the `PDL_ROOT` env var), `DATA`, `SUBSTRATE`, `OMNI_CACHE`, `OUTPUTS`, `COMMITTED`, `EVENTS` — and `P(original_path)`, which maps the original analysis-machine path strings preserved in the reproduction scripts onto this layout. User code never needs `P()`; it exists so the frozen provenance scripts run unmodified.

## Constants

`KB = 1.602e-4` (so `T[eV] = p_th[nPa] / (n[cm⁻³] · KB)`); library thresholds `TBAND/VFRAC/NFLOOR` and shell constants live in `core.py`; spectral thresholds in `spectral.py`.
