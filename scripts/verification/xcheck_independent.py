"""
Independent adversarial cross-check of the magnetosheath-coordinate recompute.

Written from first principles per the supplied recipe.  Does NOT import
compute_s_radial (or any pdl_pilot code).  The radial-s coordinate is
implemented here directly from the Shue (1998) surface geometry.

Author: independent reviewer.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import os
import numpy as np

RE_KM = 6371.2
BASE = P(r"H:\0mssl\0mssl519\data_cache\normalized")

PASSES = [
    ("P1", "usable_aug18_6h"),
    ("P2", "usable_sep03_6h"),
    ("P3", "usable_sep13_09_6h"),
    ("P4", "usable_sep20_09_6h"),
    ("P5", "usable_sep26_09_10h"),
    ("P6", "usable_sep27_09_10h"),
    ("P7", "usable_oct24_09_6h"),
]

# bin edges
NEAR_LO, NEAR_HI = 0.2, 0.4
BG_LO, BG_HI = 0.6, 1.0


def _load(d, fn):
    return np.load(os.path.join(BASE, d, fn), allow_pickle=True)


def _t_seconds(dt64):
    """Convert datetime64[ns] array to float seconds relative to its first sample."""
    t0 = dt64[0]
    return (dt64 - t0) / np.timedelta64(1, "s"), t0


def _interp_gap_masked(t_target, t_src, y_src, gap_s=60.0):
    """Linear interpolation of y_src(t_src) onto t_target.

    Points on the timeline whose nearest *bracketing* source samples are
    more than gap_s apart are set to NaN (interpolation across a data gap).
    Also NaN outside the source coverage.
    """
    t_src = np.asarray(t_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    # keep only finite source samples for the interpolation backbone
    good = np.isfinite(t_src) & np.isfinite(y_src)
    t_src = t_src[good]
    y_src = y_src[good]
    order = np.argsort(t_src)
    t_src = t_src[order]
    y_src = y_src[order]

    out = np.interp(t_target, t_src, y_src, left=np.nan, right=np.nan)

    # locate bracketing source indices for each target time
    idx_right = np.searchsorted(t_src, t_target, side="left")
    idx_right = np.clip(idx_right, 1, len(t_src) - 1)
    idx_left = idx_right - 1
    bracket_gap = t_src[idx_right] - t_src[idx_left]
    # only mask where the target actually lies between two samples (i.e. was
    # interpolated, not extrapolated); extrapolated points are already NaN.
    inside = (t_target >= t_src[0]) & (t_target <= t_src[-1])
    out[inside & (bracket_gap > gap_s)] = np.nan
    return out


def build_pass(d):
    fgm = _load(d, "thd_l2_fgm.npz")
    mom = _load(d, "thd_l2_mom.npz")
    st = _load(d, "thd_l1_state.npz")
    om = _load(d, "omni_hro_1min.npz")

    # --- source series with their own time axes (datetime64) ---
    bt = fgm["_time"]
    bvec = np.asarray(fgm["thd_fgs_gsm"], dtype=float)
    bmag_src = np.sqrt(np.sum(bvec ** 2, axis=1))

    nt = mom["_time"]
    den_src = np.asarray(mom["thd_peim_density"], dtype=float)

    pt = st["_time"]
    pos_src = np.asarray(st["thd_pos_gsm"], dtype=float) / RE_KM  # -> Re

    # --- window spanning the pass: union of B, density, position coverage ---
    t_start = min(bt[0], nt[0], pt[0])
    t_end = max(bt[-1], nt[-1], pt[-1])

    # 10-second timeline
    span_s = (t_end - t_start) / np.timedelta64(1, "s")
    n_steps = int(np.floor(span_s / 10.0)) + 1
    tl_s = np.arange(n_steps) * 10.0  # seconds from t_start
    # source times expressed in the SAME origin (t_start)
    def to_origin(arr):
        return (arr - t_start) / np.timedelta64(1, "s")

    bmag = _interp_gap_masked(tl_s, to_origin(bt), bmag_src, 60.0)
    den = _interp_gap_masked(tl_s, to_origin(nt), den_src, 60.0)
    px = _interp_gap_masked(tl_s, to_origin(pt), pos_src[:, 0], 60.0)
    py = _interp_gap_masked(tl_s, to_origin(pt), pos_src[:, 1], 60.0)
    pz = _interp_gap_masked(tl_s, to_origin(pt), pos_src[:, 2], 60.0)

    # clip density to >= 0.01
    den = np.where(np.isfinite(den), np.clip(den, 0.01, None), np.nan)

    # --- upstream solar wind medians over finite values, restricted to the
    #     spacecraft pass window (state coverage).  The reference pipeline
    #     trims the +/-30 min OMNI padding before taking the median; using the
    #     full padded array shifts mp0/bs0 enough to bias the disturbed passes.
    ot = om["_time"]
    in_win = (ot >= pt[0]) & (ot <= pt[-1])
    pres = np.asarray(om["Pressure"], dtype=float)[in_win]
    bz_arr = np.asarray(om["BZ_GSM"], dtype=float)[in_win]
    mach = np.asarray(om["Mach_num"], dtype=float)[in_win]
    dp = np.nanmedian(pres[np.isfinite(pres)])
    bz = np.nanmedian(bz_arr[np.isfinite(bz_arr)])
    ma = np.nanmedian(mach[np.isfinite(mach)])

    return dict(tl_s=tl_s, bmag=bmag, den=den, px=px, py=py, pz=pz,
                dp=dp, bz=bz, ma=ma, t_start=t_start, t_end=t_end)


def standoffs(dp, bz, ma):
    coeff = 0.013 if bz >= 0 else 0.140
    mp0 = (11.4 + coeff * bz) * dp ** (-1.0 / 6.6)
    alpha = (0.58 - 0.010 * bz) * (1.0 + 0.010 * dp)
    if ma <= 1:
        bs0 = mp0 * 1.3
    else:
        g = 5.0 / 3.0
        bs0 = mp0 * (1.0 + 1.1 * (((g - 1.0) * ma ** 2 + 2.0) / ((g + 1.0) * ma ** 2)))
    return mp0, alpha, bs0


def s_1d(px, mp0, bs0):
    """1-D coordinate: distance along X axis."""
    d_mp = np.abs(px - mp0)
    d_bs = np.abs(px - bs0)
    denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, 0.5)
    return np.clip(s, 0.0, 1.0)


def s_radial(px, py, pz, mp0, bs0, alpha):
    """My independent radial fractional coordinate.

    Shue (1998) surface: r(theta) = r0 * (2/(1+cos theta))^alpha, theta from +X.
    The magnetopause uses mp0, the bow shock is flared with the SAME alpha from
    its subsolar standoff bs0 (it has no independent shape model).  The
    spacecraft's fractional radial position between the two surfaces, at its
    actual theta, is

        s = (r - r_MP(theta)) / (r_BS(theta) - r_MP(theta)),  clipped to [0,1].
    """
    x = np.asarray(px, dtype=float)
    y = np.asarray(py, dtype=float)
    z = np.asarray(pz, dtype=float)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_t = np.where(r > 0, x / r, 1.0)
    cos_t = np.clip(cos_t, -1.0, 1.0)
    flare = (2.0 / (1.0 + cos_t)) ** alpha
    r_mp = mp0 * flare
    r_bs = bs0 * flare
    num = r - r_mp            # signed radial distance from MP
    den = r_bs - r_mp         # sheath thickness at this theta
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(den != 0, num / den, np.nan)
    return np.clip(s, 0.0, 1.0)


def s_radial_abs(px, py, pz, mp0, bs0, alpha):
    """Alternative radial form using absolute distances (mirrors the 1-D |.|)."""
    x = np.asarray(px, dtype=float); y = np.asarray(py, dtype=float); z = np.asarray(pz, dtype=float)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_t = np.where(r > 0, x / r, 1.0)
    cos_t = np.clip(cos_t, -1.0, 1.0)
    flare = (2.0 / (1.0 + cos_t)) ** alpha
    r_mp = mp0 * flare
    r_bs = bs0 * flare
    d_mp = np.abs(r - r_mp); d_bs = np.abs(r - r_bs)
    denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, 0.5)
    return np.clip(s, 0.0, 1.0)


def dn_and_occ(s, den):
    """Density ratio (near/background) and near-bin occupancy."""
    finite = np.isfinite(s)
    n_tot = int(np.sum(finite))
    near_mask = finite & (s >= NEAR_LO) & (s < NEAR_HI)
    bg_mask = finite & (s >= BG_LO) & (s < BG_HI)
    near_occ = (np.sum(near_mask) / n_tot) if n_tot > 0 else np.nan

    d_near = den[near_mask]
    d_bg = den[bg_mask]
    d_near = d_near[np.isfinite(d_near)]
    d_bg = d_bg[np.isfinite(d_bg)]
    if d_near.size == 0 or d_bg.size == 0:
        dn = np.nan
    else:
        m_bg = np.nanmedian(d_bg)
        dn = np.nanmedian(d_near) / m_bg if m_bg != 0 else np.nan
    return dn, near_occ, int(near_mask.sum()), int(bg_mask.sum()), n_tot


def fmt(v, nd=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "EMPTY/undef"
    return f"{v:.{nd}f}"


print(f"{'P':>3} {'dp':>6} {'bz':>6} {'ma':>6} {'mp0':>6} {'bs0':>6} {'alpha':>6} | "
      f"{'1D_Dn':>8} {'1D_occ':>7} | {'rad_Dn':>10} {'rad_occ':>8} {'rNear':>6} {'rBg':>6} | "
      f"{'radAbs_Dn':>10} {'radAbs_occ':>10}")
print("-" * 150)

rows = []
for name, d in PASSES:
    P = build_pass(d)
    mp0, alpha, bs0 = standoffs(P["dp"], P["bz"], P["ma"])
    s1 = s_1d(P["px"], mp0, bs0)
    sr = s_radial(P["px"], P["py"], P["pz"], mp0, bs0, alpha)
    sa = s_radial_abs(P["px"], P["py"], P["pz"], mp0, bs0, alpha)

    dn1, occ1, n1n, n1b, nt1 = dn_and_occ(s1, P["den"])
    dnr, occr, nrn, nrb, ntr = dn_and_occ(sr, P["den"])
    dna, occa, nan_, nab, nta = dn_and_occ(sa, P["den"])

    rows.append((name, P, mp0, alpha, bs0, dn1, occ1, dnr, occr, nrn, nrb, dna, occa))
    print(f"{name:>3} {P['dp']:6.3f} {P['bz']:6.2f} {P['ma']:6.2f} {mp0:6.2f} {bs0:6.2f} {alpha:6.3f} | "
          f"{fmt(dn1):>8} {occ1:7.3f} | {fmt(dnr):>10} {occr:8.4f} {nrn:6d} {nrb:6d} | "
          f"{fmt(dna):>10} {occa:10.4f}")

# extra diagnostics: theta range, r range, x range for each pass
print("\n--- geometry diagnostics (finite samples) ---")
for name, d in PASSES:
    P = build_pass(d)
    x, y, z = P["px"], P["py"], P["pz"]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    fin = np.isfinite(r)
    cos_t = np.where(r > 0, x / r, np.nan)
    th = np.degrees(np.arccos(np.clip(cos_t, -1, 1)))
    print(f"{name}: x[{np.nanmin(x):.1f},{np.nanmax(x):.1f}]Re  "
          f"r[{np.nanmin(r):.1f},{np.nanmax(r):.1f}]Re  "
          f"theta[{np.nanmin(th):.1f},{np.nanmax(th):.1f}]deg  "
          f"Nfin={int(fin.sum())}")
