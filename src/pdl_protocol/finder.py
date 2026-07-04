"""One-command candidate finder: run the full measurement protocol over a directory
of encounters and emit a table of depletion-layer candidates.

    pdl-find                      # scan data/substrate/ (after pip install -e .)
    pdl-find /path/to/encounters --csv candidates.csv

or from Python:

    from pdl_protocol import find_candidates
    rows = find_candidates()                  # list of per-encounter dicts

Per encounter: radial normalised coordinate -> magnetosheath-membership screen ->
near-shell contrast -> the conservative candidate criteria of dissertation Sec. 6.2
(density >= 2 cm^-3, beta in [0.1, 2], temperature 50-800 eV, E_B > 1.2,
D_n in [0.40, 0.90], >= 15 screened samples) -> where an ESA ion spectrum is
available (data/events/esa_<date>.npz), the per-event spectral classification of
Sec. 3.5.

Scope note (inherited from the dissertation): this finds *candidates* -- intervals
of depleted, field-enhanced, spectrally sheath-like plasma near the model
magnetopause. It does not confirm a depletion layer: placement by empirical
boundary models alone cannot certify any interval as the layer, and candidate
counts are selection-limited, not occurrence rates. Confirmation requires an
independently anchored boundary (crossing-calibrated positions, multi-spacecraft
timing, or global soft-X-ray imaging).
"""
import argparse
import csv
import os

import numpy as np

from .config import EVENTS
from .core import load_encounter, member_mask, shell_contrast
from .coords import compute_s
from .psub import pmap, list_files
from .spectral import spectral_metrics, classify_spectrum

# the near-magnetopause analysis shell (dissertation Sec. 3.3)
SHELL = (0.05, 0.20)
# the conservative candidate criteria (dissertation Sec. 6.2)
CAND_N_MIN = 2.0        # cm^-3
CAND_BETA = (0.1, 2.0)
CAND_T_EV = (50.0, 800.0)
CAND_EB_MIN = 1.2
CAND_DN = (0.40, 0.90)
CAND_MIN_SAMPLES = 15


def _spectral_tier(d, date):
    """Classify the encounter's ion spectrum if data/events/esa_<date>.npz exists."""
    p = os.path.join(EVENTS, "esa_%s.npz" % date)
    if not os.path.exists(p):
        return "", None
    esa = np.load(p)
    s = compute_s(d["x_re"], d["y_re"], d["z_re"],
                  float(d["mp0"]), float(d["alpha"]), float(d["dp"]))
    s_at = np.interp(esa["tunix"], d["t"].astype(float), s)
    near = np.nanmedian(esa["eflux"][(s_at >= SHELL[0]) & (s_at < SHELL[1])], axis=0)
    bg = np.nanmedian(esa["eflux"][(s_at >= 0.6) & (s_at <= 1.0)], axis=0)
    peak, shape, flux = spectral_metrics(near, bg, esa["energy"])
    return classify_spectrum(peak, shape), (peak, shape, flux)


def analyse_encounter(fname, d):
    """One encounter -> result dict, or None if it cannot contribute."""
    e = load_encounter(d)
    if e is None:
        return None
    m = member_mask(e)
    dn, eb = shell_contrast(e, m, *SHELL)
    if dn is None:
        return None
    sel = m & (e["s"] >= SHELL[0]) & (e["s"] < SHELL[1])
    n_samp = int(sel.sum())
    n_near = float(np.median(e["n"][sel]))
    t_near = float(np.nanmedian(e["T"][sel]))
    beta_near = float(np.nanmedian(e["beta"][sel]))
    candidate = (n_near >= CAND_N_MIN
                 and CAND_BETA[0] <= beta_near <= CAND_BETA[1]
                 and CAND_T_EV[0] <= t_near <= CAND_T_EV[1]
                 and eb > CAND_EB_MIN
                 and CAND_DN[0] <= dn <= CAND_DN[1]
                 and n_samp >= CAND_MIN_SAMPLES)
    eid = fname[:-4] if fname.lower().endswith(".npz") else fname
    date = eid.rsplit("_", 1)[0]
    tier, _metrics = _spectral_tier(d, date) if candidate else ("", None)
    if candidate and tier == "sheath-like":
        status = "spectrally supported candidate"
    elif candidate and tier in ("boundary", "borderline"):
        status = "candidate (spectrum: %s)" % tier
    elif candidate:
        status = "candidate (moments only; no spectrum available)"
    else:
        status = "-"
    bz = float(d["bz"]) if np.isfinite(d["bz"]) else np.nan
    return dict(eid=eid, D_n=round(dn, 3), E_B=round(eb, 3),
                n_near=round(n_near, 2), T_near_eV=round(t_near, 1),
                beta_near=round(beta_near, 3), n_samples=n_samp,
                bz_nT=round(bz, 2) if np.isfinite(bz) else "",
                candidate=candidate, status=status)


def find_candidates(substrate_dir=None, out_csv=None):
    """Run the protocol over every NPZ in substrate_dir (default data/substrate/).

    Returns the list of per-encounter dicts (all encounters that could contribute,
    with candidate=True/False); optionally writes them to out_csv."""
    files = list_files(substrate_dir) if substrate_dir else list_files()
    rows = pmap(analyse_encounter, files=files, with_name=True)
    rows.sort(key=lambda r: r["eid"])
    if out_csv:
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="pdl-find",
        description="Find plasma-depletion-layer candidates in a directory of "
                    "encounter NPZ files (schema: docs/DATA.md).")
    ap.add_argument("substrate", nargs="?", default=None,
                    help="directory of encounter NPZ files (default: data/substrate/)")
    ap.add_argument("--csv", default=None, help="also write the full table to this CSV")
    a = ap.parse_args(argv)

    rows = find_candidates(a.substrate, a.csv)
    if not rows:
        print("no analysable encounters found (see docs/DATA.md for the NPZ schema)")
        return 1
    hdr = "%-18s %6s %6s %8s %9s %6s  %s" % ("encounter", "D_n", "E_B", "n_near", "T_near", "N", "status")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print("%-18s %6.3f %6.3f %8.2f %9.1f %6d  %s" %
              (r["eid"], r["D_n"], r["E_B"], r["n_near"], r["T_near_eV"], r["n_samples"], r["status"]))
    n_cand = sum(r["candidate"] for r in rows)
    n_supp = sum(r["status"].startswith("spectrally supported") for r in rows)
    print("\n%d encounters analysed; %d candidates, of which %d spectrally supported."
          % (len(rows), n_cand, n_supp))
    print("Candidates are selection-limited indications, not confirmed depletion layers "
          "(see the scope note in pdl_protocol.finder).")
    if a.csv:
        print("table written to", a.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
