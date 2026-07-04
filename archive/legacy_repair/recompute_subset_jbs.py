"""Validate the Jelinek-2012-BS radial coordinate on the cached 6-pass subset before the
definitive archive run. Compares 1D, same-alpha radial, and Jelinek-BS radial Dn.
Expectation (per Agent C's blunt-BS test): P3/P4/P7 still empty; surviving Dn near or
slightly above the same-alpha radial values; nowhere does it restore deep 1D depletion.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import sys
sys.path.insert(0, P(r"H:\0mssl\0mssl519\src"))
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
import numpy as np
from pdl_pilot.config import load_config
from pdl_pilot.data.live_provider import LiveProvider
from pdl_pilot.cli.run_pilot import _apply_fill_masking, _build_upstream
from pdl_pilot.mapping.s_mapper import compute_s_with_uncertainty, compute_s_radial
from pdl_pilot.boundaries.shue1998 import shue1998_alpha
from pdl_pilot.metrics.calculator import _select_bin
from radial_models import compute_s_radial_jbs, jelinek_bs_r

CONFIG = P(r"H:\0mssl\0mssl519\configs\pilot_live_usable.yaml")
CACHE = P(r"H:\0mssl\0mssl519\data_cache")
PASS = {"usable_aug18_6h": ("P1", 0.12), "usable_sep03_6h": ("P2", 2.31), "usable_sep13_09_6h": ("P3", 0.39),
        "usable_sep20_09_6h": ("P4", 0.97), "usable_sep26_09_10h": ("P5", 0.94), "usable_sep27_09_10h": ("P6", 1.31),
        "usable_oct24_09_6h": ("P7", 2.19)}


def dn(s, dens):
    nm = _select_bin(s, 0.2, 0.4); bm = _select_bin(s, 0.6, 1.0)
    nn = np.nanmedian(dens[nm]) if nm.sum() > 0 and np.isfinite(dens[nm]).any() else np.nan
    nb = np.nanmedian(dens[bm]) if bm.sum() > 0 and np.isfinite(dens[bm]).any() else np.nan
    return ((nn / nb) if (np.isfinite(nn) and np.isfinite(nb) and nb > 0) else np.nan), float(nm.mean())


cfg = load_config(CONFIG); cfg.live.cache_dir = CACHE
prov = LiveProvider(cfg.live)
print(f"{'pass':5}{'sza':>5}{'bs0_merka':>10}{'bs0_jel':>9} | {'Dn_1D':>7}{'Dn_radSA':>9}{'occSA':>7} | {'Dn_radJBS':>10}{'occJBS':>7}")
for spec in cfg.encounters:
    pid_exp = PASS.get(spec.encounter_id)
    if not pid_exp:
        continue
    pid, _ = pid_exp
    ed = prov.fetch(spec); _apply_fill_masking(ed, "live", "auto"); up = _build_upstream(ed, "live")
    dp = up.dp_nPa or 2.0; bz = up.bz_gsm_nT or 0.0; ma = up.mach_alfven or 8.0
    x, y, z, dens = ed.x_gsm_re, ed.y_gsm_re, ed.z_gsm_re, ed.density_cm3
    s1, _, _, mp0, bs0 = compute_s_with_uncertainty(x, dp_nPa=dp, bz_nT=bz, mach_alfven=ma, unc=cfg.uncertainty)
    a = shue1998_alpha(dp, bz)
    s_sa = compute_s_radial(x, y, z, mp0, bs0, a)
    s_jbs = compute_s_radial_jbs(x, y, z, mp0, a, dp)
    bs0_jel = float(jelinek_bs_r(np.array([1.0]), dp)[0])  # subsolar Jelinek standoff
    sza = float(np.degrees(np.arctan2(np.sqrt(np.nanmean(y)**2 + np.nanmean(z)**2), np.nanmean(x))))
    d1, _ = dn(s1, dens); dsa, osa = dn(s_sa, dens); djbs, ojbs = dn(s_jbs, dens)
    print(f"{pid:5}{sza:5.1f}{bs0:10.2f}{bs0_jel:9.2f} | {d1:7.3f}{dsa:9.3f}{osa:7.3f} | {djbs:10.3f}{ojbs:7.3f}")
