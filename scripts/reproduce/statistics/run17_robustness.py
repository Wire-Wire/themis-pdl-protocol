"""RUN 17 (robustness appendix) — boundary-placement + foreshock/quasi-radial-IMF stress test.

Answers the examiner question: the result rests on a Shue magnetopause + Jelinek bow shock + radial s.
If the model boundaries are mis-placed, or quasi-radial IMF / foreshock displaces the true magnetopause
from the model, do Dn, EB, and the dynamic-pressure ordering survive? This HARDENS the existing result;
it opens no new claim. numbers-traceable: re-derives s on the FROZEN substrate under perturbed boundaries.

Part A — boundary placement (substrate pass). Per encounter, recompute the membership-screened near-MP
contrast (s in [0.05,0.20)) under fixed boundary perturbations: magnetopause standoff offset +/-0.5 R_E
(correlated worst case), bow-shock scale +/-10%, and one realistic random per-encounter draw
(offset ~ N(0,0.5) R_E, scale ~ U[0.9,1.1]). Aggregate median Dn/EB over the contributing set per
scenario, and Dn across dynamic-pressure tertiles per scenario.

Part B — foreshock / quasi-radial IMF (CSV). Recompute the 3-population contrast and the Dp tertiles
excluding low-cone-angle encounters (cone<30 deg, cone<45 deg), the quasi-parallel-shock / foreshock-prone
regime where the boundary models are least reliable.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, collections
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap, list_files

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run17_robustness")
os.makedirs(OUT, exist_ok=True)
CONTRIB = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
SPEC = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run12_spectral\spectral_metrics.csv")
BG = (0.6, 1.0); PDL = (0.05, 0.20); MIN_BG = 5; MIN_SHELL = 5
KB = 1.602e-4; N_FLOOR = 0.3; TBAND = 3.0; VFRAC = 0.2

SCEN = {'nominal': (0.0, 1.0), 'mp-0.5': (-0.5, 1.0), 'mp+0.5': (+0.5, 1.0),
        'bs-10%': (0.0, 0.90), 'bs+10%': (0.0, 1.10), 'random': None}


def contrast(d, mp_off, bs_scale):
    """Membership-screened Dn_mem, EB_mem in [0.05,0.20) under perturbed boundaries; None if not contributing."""
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = (mp0 + mp_off) * (2.0 / (1.0 + ct)) ** alpha
    r_bs = jelinek_bs_r(ct, dp) * bs_scale
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    geo = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]
    if len(sv) < MIN_BG + MIN_SHELL:
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
    if pm.sum() < MIN_SHELL:
        return None
    return float(np.median(nv[pm]) / n_bg), float(np.median(bv[pm]) / b_bg)


def one(name, d):
    dp = float(d['dp']) if np.isfinite(d['dp']) else np.nan
    rec = {'dp': dp, 'scen': {}}
    seed = abs(hash(name)) % (2 ** 32)
    rng = np.random.default_rng(seed)
    for sc, params in SCEN.items():
        if sc == 'random':
            off = float(rng.normal(0.0, 0.5)); scale = float(rng.uniform(0.90, 1.10))
        else:
            off, scale = params
        r = contrast(d, off, scale)
        if r is not None:
            rec['scen'][sc] = r
    return rec if 'nominal' in rec['scen'] else None


def med(a):
    a = np.array([v for v in a if np.isfinite(v)])
    return float(np.median(a)) if len(a) else np.nan


def tertile_dn(rows, sc):
    sub = [(r['dp'], r['scen'][sc][0]) for r in rows if sc in r['scen'] and np.isfinite(r['dp'])]
    if len(sub) < 30:
        return None
    dps = np.array([a for a, _ in sub]); q1, q2 = np.percentile(dps, [100 / 3, 200 / 3])
    lo = [dn for dp, dn in sub if dp <= q1]; mid = [dn for dp, dn in sub if q1 < dp <= q2]; hi = [dn for dp, dn in sub if dp > q2]
    return med(lo), med(mid), med(hi)


def part_a(o):
    rows = pmap(one, with_name=True)
    o.append("PART A — boundary-placement perturbation (membership Dn/EB in s[0.05,0.20), contributing set)")
    o.append(f"{'scenario':10s} {'N_contrib':>9s} {'median Dn':>10s} {'median EB':>10s}   Dn by Dp tertile (low/mid/high)")
    for sc in SCEN:
        dn = [r['scen'][sc][0] for r in rows if sc in r['scen']]
        eb = [r['scen'][sc][1] for r in rows if sc in r['scen']]
        t = tertile_dn(rows, sc)
        tstr = f"{t[0]:.3f} / {t[1]:.3f} / {t[2]:.3f}" if t else "n/a"
        o.append(f"{sc:10s} {len(dn):9d} {med(dn):10.3f} {med(eb):10.3f}   {tstr}")
    o.append("(nominal here recomputes the frozen result; perturbed rows show the effect of mis-placed boundaries.)")
    return rows


def fnum(s):
    try:
        v = float(s); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def part_b(o):
    rows = list(csv.DictReader(open(CONTRIB)))
    for r in rows:
        for k in ('Dn_mem', 'EB_mem', 'cone', 'dp'):
            r[k] = fnum(r[k])
    sheath = {r['eid'] for r in csv.DictReader(open(SPEC)) if r['status'] == 'SHEATH_CONSISTENT'}
    pops = {'all contributing': lambda r: True,
            'moment-classified (107)': lambda r: r['cls'] == 'CONFIRMED_PDL',
            'spectrally supported (60)': lambda r: r['eid'] in sheath}
    o.append("")
    o.append("PART B — foreshock / quasi-radial-IMF exclusion (low cone angle removed)")
    for cut_name, cut in (('cone>=0 (all)', 0.0), ('cone>=30', 30.0), ('cone>=45', 45.0)):
        o.append(f"  [{cut_name}]")
        for pname, sel in pops.items():
            sub = [r for r in rows if sel(r) and np.isfinite(r['cone']) and r['cone'] >= cut]
            dn = med([r['Dn_mem'] for r in sub]); eb = med([r['EB_mem'] for r in sub])
            o.append(f"    {pname:28s} N={len(sub):4d}  Dn={dn:.3f}  EB={eb:.3f}")
        # Dp tertiles on all-contributing under this cone cut
        sub = [r for r in rows if np.isfinite(r['cone']) and r['cone'] >= cut and np.isfinite(r['dp']) and np.isfinite(r['Dn_mem'])]
        dps = np.array([r['dp'] for r in sub]); q1, q2 = np.percentile(dps, [100 / 3, 200 / 3])
        lo = med([r['Dn_mem'] for r in sub if r['dp'] <= q1]); mid = med([r['Dn_mem'] for r in sub if q1 < r['dp'] <= q2]); hi = med([r['Dn_mem'] for r in sub if r['dp'] > q2])
        o.append(f"    Dp tertile Dn (all-contributing): {lo:.3f} / {mid:.3f} / {hi:.3f}  (N={len(sub)})")


def main():
    o = ["RUN 17 — boundary-placement + foreshock robustness (hardening test; no new claim)"]
    o.append("")
    part_a(o)
    part_b(o)
    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "ROBUSTNESS.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
