'''Per-event ion-spectral validation metrics (dissertation Sec. 3.5; docs/METHODS.md Sec. 4).'''
import numpy as np

PEAK_SHEATH = (0.4, 2.5)
SHAPE_SHEATH = 0.8
SHAPE_REJECT = 0.6


def spectral_metrics(near_spec, bg_spec, energies):
    """Median spectra near-shell vs background -> (peak-energy ratio, log-shape correlation, flux ratio)."""
    near = np.asarray(near_spec, float); bg = np.asarray(bg_spec, float); en = np.asarray(energies, float)
    ok = np.isfinite(near) & np.isfinite(bg) & (near > 0) & (bg > 0)
    peak = en[ok][np.argmax(near[ok])] / en[ok][np.argmax(bg[ok])]
    shape = float(np.corrcoef(np.log(near[ok]), np.log(bg[ok]))[0, 1])
    flux = float(np.trapz(near[ok], en[ok]) / np.trapz(bg[ok], en[ok]))
    return float(peak), shape, flux


def classify_spectrum(peak_ratio, shape_corr):
    """sheath-like / boundary / borderline. Thresholds are non-critical: the populations
    separate by an order of magnitude (peak ratio ~1 vs ~24)."""
    if PEAK_SHEATH[0] <= peak_ratio <= PEAK_SHEATH[1] and shape_corr >= SHAPE_SHEATH:
        return "sheath-like"
    if peak_ratio > PEAK_SHEATH[1] or shape_corr < SHAPE_REJECT:
        return "boundary"
    return "borderline"
