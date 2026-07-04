"""Definitive paired comparison from the v2 catalogue (1D vs same-alpha radial vs
proper Jelinek-2012-BS radial), on encounters evaluable under both geometries.
Run after the v2 archive run finishes.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import csv
import numpy as np

CAT = P(r"H:\0mssl\review\repair\option3\derived\archive_radial_catalogue_v2.csv")
CONE = ["quasi-radial", "low-cone", "intermediate", "perpendicular"]


def F(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


rows = list(csv.DictReader(open(CAT)))
print(f"catalogue rows: {len(rows)}")

# full-archive evaluable counts + Dn medians per geometry
print("\n=== FULL-ARCHIVE (per-geometry evaluable set) ===")
for g in ("1d", "rad", "jbs"):
    ev = [F(r[f"Dn_{g}"]) for r in rows if str(r.get(f"evaluable_{g}")) == "True"]
    ev = [x for x in ev if x and x > 0]
    n = sum(1 for r in rows if str(r.get(f"evaluable_{g}")) == "True")
    print(f"  {g:4s}: N_evaluable={n:4d}  Dn median={np.median(ev):.3f}  (frac Dn<1 = {np.mean(np.array(ev)<1):.1%})")


def paired(g):
    P = []
    for r in rows:
        if str(r.get("evaluable_1d")) == "True" and str(r.get(f"evaluable_{g}")) == "True":
            d1 = F(r["Dn_1d"]); dg = F(r[f"Dn_{g}"]); e1 = F(r["EB_1d"]); eg = F(r[f"EB_{g}"]); sza = F(r["sza_deg"])
            if d1 and dg and d1 > 0 and dg > 0:
                P.append((r["cone_bin"], d1, dg, e1, eg, sza))
    n = len(P)
    d1 = np.array([p[1] for p in P]); dg = np.array([p[2] for p in P])
    e1 = np.array([p[3] for p in P]); eg = np.array([p[4] for p in P])
    print(f"\n=== PAIRED 1D vs {g.upper()} (evaluable in BOTH): N={n} ===")
    print(f"  Dn median: 1D={np.median(d1):.3f} -> {g}={np.median(dg):.3f}   (median per-enc ratio {np.median(dg/d1):.3f})")
    print(f"  EB median: 1D={np.median(e1):.3f} -> {g}={np.median(eg):.3f}")
    print(f"  frac Dn<1: 1D={np.mean(d1<1):.1%} -> {g}={np.mean(dg<1):.1%};  flip <1->>=1: {np.mean((d1<1)&(dg>=1)):.1%}")
    low = d1 < 0.5
    print(f"  Dn_1d<0.5 'strong depletion' (N={low.sum()}): {g} median={np.median(dg[low]):.3f}, frac>=0.8={np.mean(dg[low]>=0.8):.1%}")
    for cb in CONE:
        m = np.array([p[0] == cb for p in P])
        if m.sum():
            print(f"     {cb:14s} N={m.sum():4d}  1D={np.median(d1[m]):.3f} -> {g}={np.median(dg[m]):.3f}")


for g in ("rad", "jbs"):
    paired(g)
