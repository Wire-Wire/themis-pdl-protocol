"""Proper bow-shock surfaces for the definitive radial recompute.

Jelinek et al. (2012, WDS'10 Proc. / JGR 2012) bow shock — a parabolic surface fit
& validated on THEMIS dayside crossings 2007-2009 (the exact mission/era of our data),
with an INDEPENDENT shape from the magnetopause (its own standoff exponent and flaring).

BS surface (parabolic, sigma=const):  sqrt(x^2 + lambda^2 (y^2+z^2)) + x = 2 R0 Dp^(-1/eps)
Subsolar standoff R_BS = 15.02 * Dp^(-1/6.55) Re;  lambda_BS = 1.17.
Radial distance at angle theta:  r_BS(theta) = 2 R0 Dp^(-1/eps) / (cos th + sqrt(cos^2 th + lambda^2 sin^2 th)).
At theta=0 this reduces to R_BS (verified). Magnetopause kept as exact Shue-1998 (pipeline mp0, alpha).
"""
from __future__ import annotations
import numpy as np

R0_BS_JEL = 15.02
EPS_BS_JEL = 6.55
LAMBDA_BS_JEL = 1.17


def jelinek_bs_r(cos_theta: np.ndarray, dp_nPa: float) -> np.ndarray:
    """Radial distance (Re) to the Jelinek-2012 bow shock at angle theta (given cos theta), for Dp."""
    ct = np.clip(np.asarray(cos_theta, float), -1.0, 1.0)
    sin2 = np.clip(1.0 - ct ** 2, 0.0, 1.0)
    standoff_times2 = 2.0 * R0_BS_JEL * dp_nPa ** (-1.0 / EPS_BS_JEL)
    denom = ct + np.sqrt(ct ** 2 + (LAMBDA_BS_JEL ** 2) * sin2)
    return standoff_times2 / denom


def compute_s_radial_jbs(x, y, z, mp0: float, alpha: float, dp_nPa: float) -> np.ndarray:
    """Radial normalised coordinate with exact Shue MP + Jelinek-2012 BS (two independent surfaces)."""
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha          # Shue 1998 magnetopause
    r_bs = jelinek_bs_r(ct, dp_nPa)                    # Jelinek 2012 bow shock
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, 0.5)
    return np.clip(s, 0.0, 1.0)
