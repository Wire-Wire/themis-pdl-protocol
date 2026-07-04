'''Boundary models and the radial normalised magnetosheath coordinate (single source of truth).
Equations and references: docs/METHODS.md (Shue et al. 1998; Jelinek et al. 2012).'''
import numpy as np
from .radial_models import jelinek_bs_r


def shue_r0(dp, bz):
    """Shue-1998 subsolar magnetopause standoff (R_E) from Dp (nPa) and IMF Bz (nT)."""
    return (10.22 + 1.29 * np.tanh(0.184 * (bz + 8.14))) * dp ** (-1.0 / 6.6)


def shue_alpha(dp, bz):
    """Shue-1998 flaring parameter."""
    return (0.58 - 0.007 * bz) * (1.0 + 0.024 * np.log(dp))


def shue_r(cos_theta, r0, alpha):
    """Shue magnetopause radius at angle theta from the Sun-Earth line."""
    return r0 * (2.0 / (1.0 + cos_theta)) ** alpha


def jelinek_r(cos_theta, dp):
    """Jelinek-2012 bow-shock radius at angle theta (R0=15.02, eps=6.55, lam=1.17)."""
    return jelinek_bs_r(cos_theta, dp)


def compute_s(x_re, y_re, z_re, mp0, alpha, dp):
    """Radial normalised coordinate s = d_MP/(d_MP+d_BS); 0 at the magnetopause, 1 at the bow shock.
    Identical to the construction used throughout the reproduction chain."""
    x = np.asarray(x_re, float); y = np.asarray(y_re, float); z = np.asarray(z_re, float)
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    d_mp = r - shue_r(ct, mp0, alpha)
    d_bs = jelinek_r(ct, dp) - r
    den = d_mp + d_bs
    return np.where(den > 0, d_mp / den, np.nan)
