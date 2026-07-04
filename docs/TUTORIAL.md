# Tutorial: using the protocol on data (the researcher path)

This is **Door B**: you want to find depletion-layer candidates or measure near-boundary contrasts in magnetosheath data. You will not need the dissertation's reproduction chain; if you are here to audit the dissertation instead, go to [REPRODUCING.md](REPRODUCING.md).

Five levels, each self-contained, each with the real expected output. Levels 0–4 need **no downloads** and finish in minutes.

## Level 0: install (once, ~1 minute)

```bash
git clone https://github.com/Wire-Wire/themis-pdl-protocol
cd themis-pdl-protocol
pip install -e .            # library core: numpy + scipy only; registers pdl-find / pdl-contrast
```

Optional extras when you need them:

| Command | Adds | Needed for |
|---|---|---|
| `pip install -e ".[plot]"` | matplotlib | figure scripts |
| `pip install -e ".[fetch]"` | cdasws, xarray | `tools/get_data.py` (CDAWeb downloads) |

## Level 1: one encounter through the protocol (10 seconds)

```bash
python examples/01_one_encounter.py
```

```text
THEMIS-A 2020-12-12: D_n = 0.554, E_B = 2.279
committed values: 0.554 / 2.279 -> a depleted-density, enhanced-field candidate interval
```

What the ~10 lines of code did: loaded a shipped 6-hour encounter (`data/substrate/2020-12-12_tha.npz`), computed the radial normalised coordinate s from the encounter's own Shue/Jelínek boundary parameters, applied the magnetosheath-membership screen, and took near-shell (s ∈ [0.05, 0.20)) to background (s ∈ [0.6, 1.0]) median ratios. That is the entire per-encounter measurement.

## Level 2: classify the ion spectrum (10 seconds)

```bash
python examples/02_classify_spectrum.py
```

```text
peak ratio 1.00, shape corr 0.96, flux ratio 0.46 -> sheath-like
```

Reading: the near-boundary ions keep the background-sheath peak energy and spectral shape while the flux falls: the signature of genuinely depleted *magnetosheath* plasma. A magnetospheric / boundary-layer contaminant jumps to keV (peak ratio ≈ 24 in the dissertation's rejected population). The example recomputes the shell windows with a simple interpolation, so its flux ratio (0.46) differs slightly from the pipeline's exact value (0.39); peak and shape reproduce exactly, and the classification is identical.

## Level 3: your own data (5 minutes)

```bash
python examples/03_bring_your_own_data.py
```

```text
synthetic encounter: D_n = 0.50, E_B = 2.02  (injected: 0.5 / 2.0)
```

The example builds a synthetic encounter with a known injected depletion and recovers it through the full screen+contrast loop. To use real data: create one NPZ per encounter matching the schema in [DATA.md](DATA.md) (positions, ion moments, |B|, per-encounter OMNI drivers and Shue/Jelínek parameters), name it `YYYY-MM-DD_prb.npz`, and drop it into `data/substrate/` (or keep your own directory and pass it to the tools below).

## Level 4: the one-command tools

**Find candidates** (the dissertation's §6 chain distilled into one command):

```bash
pdl-find                          # scans data/substrate/
pdl-find /your/encounters --csv candidates.csv
```

```text
encounter             D_n    E_B   n_near    T_near      N  status
------------------------------------------------------------------
2011-12-18_the      0.526  1.998    52.97     101.8    850  candidate (moments only; no spectrum available)
2020-12-12_tha      0.554  2.279    13.48     172.6    723  spectrally supported candidate
2024-04-18_tha      0.607  2.845    19.08     123.4    781  candidate (moments only; no spectrum available)

3 encounters analysed; 3 candidates, of which 1 spectrally supported.
```

Per encounter: radial coordinate → membership screen → near-shell contrast → the conservative candidate criteria (density ≥ 2 cm⁻³, β ∈ [0.1, 2], T ∈ [50, 800] eV, E_B > 1.2, D_n ∈ [0.40, 0.90], ≥ 15 screened samples) → spectral classification wherever `data/events/esa_<date>.npz` exists (fetch one with `python tools/get_data.py --event <eid>`). The three shipped demo encounters were chosen as candidate events, so 3/3 here is expected, not typical.

**Analyse the contrast** (profile, population statistic, and the coordinate-artefact check):

```bash
pdl-contrast --paired-check
```

```text
radial shell profile (across-encounter medians):
s shell           N     D_n     E_B    beta
[0.00, 0.05)     3   0.350   2.981   0.149
[0.05, 0.10)     3   0.373   2.820   0.215
...
near-shell population statistic, s in [0.05, 0.20):
  N = 3 encounters | D_n = 0.554 IQR [0.54, 0.581] | E_B = 2.279 IQR [2.139, 2.562] | depleted fraction 100%
```

`--paired-check` re-runs the dissertation's §4 design on *your* data: the same encounters scored under the radial and the fixed Sun–Earth-axis constructions. A large shift means the coordinate choice, not the plasma, is producing the apparent depletion. (It needs ≥ 3 encounters evaluable under both constructions; on the 3-encounter demo it reports "too few".)

**Scope note that travels with both tools:** they report *candidates* — depleted, field-enhanced, spectrally sheath-like intervals near the *model* magnetopause. They do not confirm a depletion layer: model placement alone cannot certify any interval as the layer, and candidate counts are selection-limited, not occurrence rates.

## Level 5: hunt across the full THEMIS archive (hours, ~2 GB)

```bash
pip install -e ".[fetch]"
python tools/get_data.py --full-substrate     # prints the build command + expectations
python scripts/reproduce/build_substrate.py   # THEMIS FGM/ESA + OMNI, 2007-2025 -> data/substrate/*.npz
pdl-find --csv all_candidates.csv             # then hunt the whole mission
```

The builder fetches from the live CDAWeb archive; if the archive re-calibrates data, per-encounter values can differ at the margin from the dissertation's frozen reference. For interpreting population-level results (what is generic compression vs depletion-specific, why candidate counts are not occurrence rates), read [METHODS.md](METHODS.md) sections 4–6.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ModuleNotFoundError: pdl_protocol` | run `pip install -e .` from the repo root, or run examples from the repo root so the `src/` fallback works |
| `pdl-find: command not found` | the console scripts are registered by `pip install -e .`; without installing, use `python -m pdl_protocol.finder` |
| Figures fail with a display/backend error | headless environments: `MPLBACKEND=Agg python …` |
| `cdasws` missing | only needed for downloads: `pip install -e ".[fetch]"` |
| You want the data elsewhere | pass the directory to `pdl-find`/`pdl-contrast`, or set `PDL_ROOT=/path` (see `src/pdl_protocol/config.py`) |
