"""RUN 9a — refine the "confirmed-PDL" event set with ABSOLUTE sheath-plasma criteria.

The atlas showed the deepest-Dn candidates are often MP-crossing / boundary-layer (beta~0.02,
T~keV, n~0.5) — low-beta + hot + tenuous, not clean shocked sheath. A relative-to-background
T screen passes them when the background is also hot. So classify each event's near-MP layer
(membership samples in the PDL shell [0.05,0.20)) by ABSOLUTE plasma values:

  CONFIRMED_PDL (northward, clean sheath PDL):
     Bz>0 AND Dn_mem in [0.40,0.90] (depleted but not collapsed)
     AND n_near >= 2 cm^-3 (sheath density, not magnetosphere-tenuous)
     AND 0.10 <= beta_near <= 2.0 (depleted but sheath, not beta<<1 boundary layer)
     AND 50 <= T_near <= 800 eV (shocked sheath, not keV boundary/magnetosphere)
     AND EB_mem > 1.2 AND n_mem >= 15
  BOUNDARY_LAYER:  beta_near < 0.10 OR n_near < 1.0 OR T_near > 1200 eV
  AMBIGUOUS:       otherwise

Writes confirmed_pdl_candidates.csv (with time windows for the Run-9b spectral confirmation).
Threaded via psub.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, datetime as dt
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run9_candidates")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0); PDL = (0.05, 0.20); MIN_BG = 5; KB = 1.602e-4
N_FLOOR = 0.3; TBAND = 3.0; VFRAC = 0.2


def one(name, d):
    eid = name[:-4]; probe = eid.rsplit('_', 1)[-1]
    t = d['t'].astype(float)
    if len(t) < 10:
        return None
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float); beta = d['beta'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    bz = float(d['bz']) if np.isfinite(d['bz']) else np.nan
    cone = float(d['cone_deg']) if np.isfinite(d['cone_deg']) else np.nan
    sza = float(d['sza_deg']) if np.isfinite(d['sza_deg']) else np.nan
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha; r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    geo = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]; betav = beta[geo]
    if len(sv) < MIN_BG + 3:
        return None
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return None
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm]); T_bg = np.nanmedian(Tv[bgm]); v_bg = np.nanmedian(vv[bgm])
    if not (n_bg > 0 and b_bg > 0 and np.isfinite(T_bg) and T_bg > 0):
        return None
    mem = (nv > N_FLOOR) & np.isfinite(Tv) & (Tv > T_bg / TBAND) & (Tv < T_bg * TBAND)
    if np.isfinite(v_bg) and v_bg > 0:
        mem = mem & np.isfinite(vv) & (vv > VFRAC * v_bg)
    pm = mem & (sv >= PDL[0]) & (sv < PDL[1])
    if pm.sum() < 5:
        return None
    return dict(eid=eid, probe=probe, t0=float(t.min()), t1=float(t.max()), bz=bz, cone=cone, sza=sza, dp=dp,
                n_mem=int(pm.sum()), Dn_mem=float(np.median(nv[pm]) / n_bg), EB_mem=float(np.median(bv[pm]) / b_bg),
                n_near=float(np.median(nv[pm])), beta_near=float(np.nanmedian(betav[pm])),
                T_near=float(np.nanmedian(Tv[pm])), T_bg=float(T_bg), n_bg=float(n_bg))


def classify(r):
    if r['beta_near'] < 0.10 or r['n_near'] < 1.0 or r['T_near'] > 1200:
        return 'BOUNDARY_LAYER'
    if (np.isfinite(r['bz']) and r['bz'] > 0 and 0.40 <= r['Dn_mem'] <= 0.90 and r['n_near'] >= 2.0
            and 0.10 <= r['beta_near'] <= 2.0 and 50 <= r['T_near'] <= 800 and r['EB_mem'] > 1.2 and r['n_mem'] >= 15):
        return 'CONFIRMED_PDL'
    return 'AMBIGUOUS'


def isofmt(ep):
    return dt.datetime.utcfromtimestamp(ep).strftime("%Y-%m-%dT%H:%M:%SZ")


def main():
    res = [r for r in pmap(one, with_name=True) if r is not None]
    for r in res:
        r['cls'] = classify(r)
    nor = [r for r in res if np.isfinite(r['bz']) and r['bz'] > 0]
    from collections import Counter
    print(f"events with PDL-layer metrics: {len(res)}  (northward: {len(nor)})", flush=True)
    print("classification (northward):", dict(Counter(r['cls'] for r in nor)), flush=True)
    print("classification (all):", dict(Counter(r['cls'] for r in res)), flush=True)
    conf = [r for r in res if r['cls'] == 'CONFIRMED_PDL']
    conf.sort(key=lambda r: -r['n_mem'])           # best-sampled first
    path = os.path.join(OUT, "confirmed_pdl_candidates.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rank', 'eid', 'probe', 'date', 't_start', 't_end', 'Dn_mem', 'EB_mem', 'n_near', 'beta_near',
                    'T_near', 'T_bg', 'n_mem', 'bz', 'cone', 'sza', 'dp'])
        for i, r in enumerate(conf):
            w.writerow([i + 1, r['eid'], r['probe'], r['eid'].rsplit('_', 1)[0], isofmt(r['t0']), isofmt(r['t1']),
                        round(r['Dn_mem'], 3), round(r['EB_mem'], 3), round(r['n_near'], 2), round(r['beta_near'], 3),
                        round(r['T_near'], 0), round(r['T_bg'], 0), r['n_mem'], round(r['bz'], 1),
                        round(r['cone'], 1), round(r['sza'], 1), round(r['dp'], 2)])
    print(f"\nCONFIRMED_PDL: {len(conf)} events -> {path}", flush=True)
    print(f"\n{'eid':22s} {'Dn':>5s} {'EB':>5s} {'n_near':>6s} {'beta':>5s} {'T_near':>6s} {'n':>4s} {'bz':>5s}", flush=True)
    for r in conf[:25]:
        print(f"{r['eid']:22s} {r['Dn_mem']:5.2f} {r['EB_mem']:5.2f} {r['n_near']:6.2f} {r['beta_near']:5.2f} "
              f"{r['T_near']:6.0f} {r['n_mem']:4d} {r['bz']:5.1f}", flush=True)


if __name__ == "__main__":
    main()
