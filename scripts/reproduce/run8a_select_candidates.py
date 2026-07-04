"""RUN 8a — select ~25 candidate events for the per-event spectra atlas.

20 strongest NORTHWARD-IMF, membership-verified PDL candidates + 5 southward-IMF controls
(no PDL, same favourable geometry). Outputs candidates.csv with each event's time window
(for the Run-8b spectra re-fetch). Local/threaded via psub.

PDL strength = lowest Dn_mem in the PDL shell [0.05,0.20) among membership-confirmed sheath,
under Bz>0, cone>45, SZA<20, 1<Dp<4, with >=10 member samples and EB_mem>1.
Controls = Bz<0, same geometry, Dn_mem in [0.9,1.1], well-sampled.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, datetime as dt
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run8_atlas")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0); PDL = (0.05, 0.20); MIN_BG = 5; KB = 1.602e-4
N_FLOOR = 0.3; TBAND = 3.0; VFRAC = 0.2; MIN_MEM_PDL = 10


def one(name, d):
    eid = name[:-4]; probe = eid.rsplit('_', 1)[-1]
    t = d['t'].astype(float)
    if len(t) < 10:
        return None
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float)
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
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]
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
    rec = dict(eid=eid, probe=probe, t0=float(t.min()), t1=float(t.max()),
               bz=bz, cone=cone, sza=sza, dp=dp, n_mem=int(pm.sum()))
    if pm.sum() >= 3:
        rec['Dn_mem'] = float(np.median(nv[pm]) / n_bg)
        rec['EB_mem'] = float(np.median(bv[pm]) / b_bg)
    return rec


def isofmt(ep):
    return dt.datetime.utcfromtimestamp(ep).strftime("%Y-%m-%dT%H:%M:%SZ")


def main():
    res = [r for r in pmap(one, with_name=True) if r is not None and 'Dn_mem' in r]
    print(f"encounters with PDL-shell metrics: {len(res)}", flush=True)

    def fav(r):
        return (np.isfinite(r['bz']) and np.isfinite(r['cone']) and np.isfinite(r['sza'])
                and r['cone'] > 45 and r['sza'] < 20 and 1.0 < r['dp'] < 4.0 and r['n_mem'] >= MIN_MEM_PDL)

    pdl = [r for r in res if fav(r) and r['bz'] > 0 and r.get('EB_mem', 0) > 1.0]
    pdl.sort(key=lambda r: r['Dn_mem'])           # strongest depletion first
    pdl = pdl[:20]
    ctrl = [r for r in res if fav(r) and r['bz'] < 0 and 0.9 <= r['Dn_mem'] <= 1.1]
    ctrl.sort(key=lambda r: -r['n_mem'])
    ctrl = ctrl[:5]

    rows = [('PDL', i + 1, r) for i, r in enumerate(pdl)] + [('CONTROL', i + 1, r) for i, r in enumerate(ctrl)]
    path = os.path.join(OUT, "candidates.csv")
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['category', 'rank', 'eid', 'probe', 'date', 't_start', 't_end', 't0_epoch', 't1_epoch',
                    'Dn_mem', 'EB_mem', 'n_mem', 'bz', 'cone', 'sza', 'dp'])
        for cat, rk, r in rows:
            w.writerow([cat, rk, r['eid'], r['probe'], r['eid'].rsplit('_', 1)[0],
                        isofmt(r['t0']), isofmt(r['t1']), round(r['t0'], 1), round(r['t1'], 1),
                        round(r['Dn_mem'], 3), round(r.get('EB_mem', float('nan')), 3), r['n_mem'],
                        round(r['bz'], 1), round(r['cone'], 1), round(r['sza'], 1), round(r['dp'], 2)])
    print(f"\nselected {len(pdl)} PDL + {len(ctrl)} control -> {path}\n", flush=True)
    print(f"{'cat':8s} {'eid':22s} {'Dn_mem':>7s} {'EB_mem':>7s} {'n':>4s} {'bz':>6s} {'cone':>5s} {'sza':>5s}", flush=True)
    for cat, rk, r in rows:
        print(f"{cat:8s} {r['eid']:22s} {r['Dn_mem']:7.3f} {r.get('EB_mem',float('nan')):7.3f} {r['n_mem']:4d} "
              f"{r['bz']:6.1f} {r['cone']:5.0f} {r['sza']:5.0f}", flush=True)


if __name__ == "__main__":
    main()
