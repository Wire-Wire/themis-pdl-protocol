# Data — acquisition, demo subset, and bring-your-own

All spacecraft data are public (NASA CDAWeb): THEMIS FGM + ESA [Ang08, Aus08, McF08] and OMNI [KP05].

## What ships with the repo

- `committed_outputs/` — every dissertation number's source (TXT/CSV). `verify_key_numbers.py` re-derives the headlines from these alone (no downloads).
- `data/substrate/` — a **3-encounter demo subset** (`2020-12-12_tha`, `2011-12-18_the`, `2024-04-18_tha`: the deep-dive event + two atlas events) so the event-level pipeline runs immediately.
- `data/events/esa_2020-12-12.npz`, `omni_2020-12-12.npz` — the re-fetched ESA spectrogram + 1-min OMNI for the deep-dive figure (`make_fig8_deepdive.py`); re-fetchable with `python tools/get_data.py --event 2020-12-12_tha`.

## Auto-download (tools/get_data.py)

```bash
python tools/get_data.py --check            # what is present
python tools/get_data.py --event 2020-12-12_tha   # OMNI + ESA spectrogram for one encounter
python tools/get_data.py --omni             # 1-min OMNI for every encounter in data/substrate (runs 27/29/31)
python tools/get_data.py --full-substrate   # guidance for the full 6,248-encounter rebuild (hours, ~2 GB)
```

The full rebuild (`scripts/reproduce/build_substrate.py`) fetches THEMIS moments + OMNI drivers for 2007–2025 and writes one NPZ per encounter (6-h dayside window per probe per day, centred on the day's minimum solar-zenith-angle sample within the orbital screen SZA < 30°, 8 < r < 25 R_E; multiple back-and-forth crossings within a window = ONE encounter).

## Bring your own data (substrate NPZ schema)

Drop files named `YYYY-MM-DD_prb.npz` into `data/substrate/`; the whole chain runs on whatever is present. Required keys:

| Key | Shape | Meaning |
|---|---|---|
| `t` | (N,) | Unix seconds |
| `x_re, y_re, z_re` | (N,) | GSE position, R_E |
| `n` | (N,) | ion density, cm⁻³ |
| `bmag` | (N,) | \|B\|, nT |
| `beta, p_th, p_b, p_tot` | (N,) | plasma β and pressures (nPa) |
| `vmag` | (N,) | ion bulk speed, km/s |
| `dp, bz, ma, bx, cone_deg, clock_deg, sza_deg` | scalar | encounter-averaged OMNI drivers + solar-zenith angle |
| `mp0, alpha` | scalar | Shue r0 (R_E) and flaring α for the encounter |
| `bs0` | scalar | Jelínek subsolar standoff (R_E) |
| `r_re, x_mean` | scalar | mean radial distance / X (R_E) |
| `probe, date` | scalar str | e.g. `tha`, `2020-12-12` |

Set `PDL_ROOT=/elsewhere` to relocate `data/` and `outputs/` without touching the code.
