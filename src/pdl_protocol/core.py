"""Core protocol operations: load an encounter, apply the magnetosheath-membership
screen, compute shell contrasts.

These are the exact operations used by every reproduction script (extracted verbatim
from the pipeline's shared loader so `pip install`ed users get them without path
tricks). Dissertation anchors: coordinate §3.1–3.2, contrasts §3.3, membership §3.4.

    >>> import numpy as np
    >>> from pdl_protocol import load_encounter, member_mask, shell_contrast
    >>> d = np.load("data/substrate/2020-12-12_tha.npz", allow_pickle=True)
    >>> e = load_encounter(d)
    >>> dn, eb = shell_contrast(e, member_mask(e), 0.05, 0.20)
"""
import numpy as np

from .radial_models import jelinek_bs_r

# Baseline thresholds (dissertation §3.3–3.4; sensitivity: committed run25 / Appendix F)
BG = (0.6, 1.0)        # background shell in s (outer sheath)
MIN_SAMP = 3           # minimum samples for a shell median (pipeline stages apply their own stricter minima)
MIN_BG = 5             # minimum background samples per encounter
KB = 1.602e-4          # Boltzmann factor for T[eV] = p_th[nPa] / (n[cm^-3] * KB)
TBAND = 3.0            # membership: temperature within a factor TBAND of the background
VFRAC = 0.2            # membership: flow speed above VFRAC x background flow
NFLOOR = 0.3           # membership: density floor, cm^-3


def load_encounter(d):
    """Turn one substrate NPZ (schema: docs/DATA.md) into an analysis-ready dict.

    Computes the radial normalised coordinate s from the per-encounter Shue
    parameters (mp0, alpha) and the Jelínek bow shock (dp), keeps the samples that
    are geometrically inside the model magnetosheath, and attaches the encounter's
    own background medians (n_bg, b_bg, T_bg, v_bg over s in [0.6, 1.0]).

    Returns None if the encounter has no valid boundary solution or too few
    samples — the same silent-skip rule as the pipeline.
    """
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float); beta = d['beta'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    bz_enc = float(d['bz']) if np.isfinite(d['bz']) else np.nan
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    geo = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]; betav = beta[geo]
    if len(sv) < MIN_BG + MIN_SAMP:
        return None
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return None
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm])
    T_bg = np.nanmedian(Tv[bgm]); v_bg = np.nanmedian(vv[bgm])
    if not (n_bg > 0 and b_bg > 0 and np.isfinite(T_bg) and T_bg > 0):
        return None
    return dict(s=sv, n=nv, b=bv, T=Tv, v=vv, beta=betav, bz=bz_enc,
                n_bg=n_bg, b_bg=b_bg, T_bg=T_bg, v_bg=v_bg)


def member_mask(e, TBAND=TBAND, VFRAC=VFRAC, NFLOOR=NFLOOR):
    """Magnetosheath-membership screen (dissertation §3.4).

    A sample counts as magnetosheath if it is dense (n > NFLOOR cm^-3), at a
    sheath-like temperature (within a factor TBAND of the encounter's background
    temperature) and flowing (v > VFRAC x background flow). Thresholds are
    deliberately permissive; sensitivity in committed run25 (Appendix F).
    """
    m = (e['n'] > NFLOOR) & np.isfinite(e['T']) & (e['T'] > e['T_bg'] / TBAND) & (e['T'] < e['T_bg'] * TBAND)
    if np.isfinite(e['v_bg']) and e['v_bg'] > 0:
        m = m & np.isfinite(e['v']) & (e['v'] > VFRAC * e['v_bg'])
    return m


def shell_contrast(e, mask, lo, hi):
    """Near-shell contrast over s in [lo, hi): (D_n, E_B) relative to the encounter's
    own background medians, or (None, None) if the shell has fewer than MIN_SAMP
    surviving samples (dissertation §3.3)."""
    sel = mask & (e['s'] >= lo) & (e['s'] < hi)
    if sel.sum() < MIN_SAMP:
        return None, None
    return float(np.median(e['n'][sel]) / e['n_bg']), float(np.median(e['b'][sel]) / e['b_bg'])


# pipeline-era aliases (the reproduction scripts use these names)
load_enc = load_encounter
shell_dn = shell_contrast
