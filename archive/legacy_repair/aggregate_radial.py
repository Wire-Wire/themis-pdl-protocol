"""Aggregate the dual-geometry archive catalogue into the manuscript-level numbers
(full cone medians, Dp-matched Table 5, Tier-A Table 6, (Dn,EB) IQR box) for BOTH the
1D geometry (cross-check vs frozen V8.6 values) and the radial geometry (the repair).
Replicates the documented recipe (Dp band [1.0,2.5), tier rubric, bootstrap seed 42 B=1000,
IQR 25/75). Run after recompute_archive_radial.py finishes.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import csv, os
import numpy as np

CAT = P(r"H:\0mssl\review\repair\option3\derived\archive_radial_catalogue_v2.csv")
CONE = ["quasi-radial", "low-cone", "intermediate", "perpendicular"]
# Frozen V8.6 values (from derived_*.csv, per comprehension spec)
FROZEN_FULL = {"quasi-radial": 0.891, "low-cone": 0.795, "intermediate": 0.785, "perpendicular": 0.735}
FROZEN_DPM = {"quasi-radial": 0.674, "low-cone": 0.824, "intermediate": 0.871, "perpendicular": 0.798}
FROZEN_TIERA = {"quasi-radial": 0.177, "low-cone": 0.179, "intermediate": 0.161, "perpendicular": 0.229}
FROZEN_TIERA_N = {"quasi-radial": 1, "low-cone": 14, "intermediate": 62, "perpendicular": 116}
FROZEN_BOX = {"Dn": (0.0607, 0.5854), "EB": (2.1276, 4.1846)}


def F(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def boot_ci(vals, B=1000, seed=42):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], float)
    if len(v) == 0:
        return (np.nan, np.nan, np.nan)
    med = float(np.median(v))
    if len(v) < 3:
        return (med, np.nan, np.nan)
    rng = np.random.default_rng(seed)
    b = np.array([np.median(rng.choice(v, len(v), replace=True)) for _ in range(B)])
    lo, hi = np.percentile(b, [2.5, 97.5])
    return (med, float(lo), float(hi))


def build(rows, geom):
    out = []
    for r in rows:
        if str(r.get(f"evaluable_{geom}")) != "True":
            continue
        dn = F(r.get(f"Dn_{geom}")); eb = F(r.get(f"EB_{geom}"))
        no = F(r.get(f"near_occ_{geom}")); bo = F(r.get(f"bg_occ_{geom}"))
        dp = F(r.get("dp_nPa")); cb = r.get("cone_bin")
        if dn is None or eb is None or cb not in CONE:
            continue
        qtc = "clean" if (dn < 1 and eb > 1) else ("mixed" if (dn < 1 or eb > 1) else "unclear")
        qd = "undisturbed" if (no and bo and no > 0.1 and bo > 0.05) else "uncertain"
        qbm = "stable" if (dp and 2 <= dp <= 6) else "uncertain"
        if qtc == "clean" and qd == "undisturbed" and qbm == "stable":
            tier = "A_best"
        elif qtc == "unclear" or qd == "uncertain":
            tier = "C_lowest"
        else:
            tier = "B_middle"
        out.append({"cb": cb, "Dn": dn, "EB": eb, "dp": dp, "tier": tier})
    return out


def med_pos(data, cb, key="Dn"):
    v = [d[key] for d in data if d["cb"] == cb and d[key] > 0]
    return (float(np.median(v)) if v else np.nan), len(v)


def main():
    rows = list(csv.DictReader(open(CAT)))
    print(f"catalogue rows: {len(rows)}")
    for geom in ("1d", "rad", "jbs"):
        data = build(rows, geom)
        print(f"\n{'='*64}\nGEOM = {geom}   evaluable N = {len(data)}\n{'='*64}")
        cc = {cb: sum(1 for d in data if d["cb"] == cb) for cb in CONE}
        print("cone counts:", cc)
        print("\nFull-archive Dn medians:")
        for cb in CONE:
            m, n = med_pos(data, cb)
            extra = f"  (frozen {FROZEN_FULL[cb]})" if geom == "1d" else ""
            print(f"   {cb:14s} N={n:4d}  Dn_med={m:.3f}{extra}")
        print("\nDp-matched [1.0,2.5) Dn medians (Table 5):")
        dpm = [d for d in data if d["dp"] and 1.0 <= d["dp"] < 2.5]
        for cb in CONE:
            v = [d["Dn"] for d in dpm if d["cb"] == cb and d["Dn"] > 0]
            m, lo, hi = boot_ci(v)
            extra = f"  (frozen {FROZEN_DPM[cb]})" if geom == "1d" else ""
            print(f"   {cb:14s} N={len(v):4d}  med={m:.3f} CI=[{lo:.3f},{hi:.3f}]{extra}")
        tA = [d for d in data if d["tier"] == "A_best"]
        tiers = {t: sum(1 for d in data if d["tier"] == t) for t in ("A_best", "B_middle", "C_lowest")}
        print(f"\nTier rubric: {tiers}  (high-quality = A_best N={len(tA)})")
        print("Tier-A Dn medians (Table 6):")
        for cb in CONE:
            m, n = med_pos(tA, cb)
            extra = f"  (frozen {FROZEN_TIERA[cb]}, N={FROZEN_TIERA_N[cb]})" if geom == "1d" else ""
            print(f"   {cb:14s} N={n:4d}  Dn_med={m:.3f}{extra}")
        if tA:
            ldn = np.log10([max(d["Dn"], 1e-3) for d in tA if d["Dn"] > 0])
            leb = np.log10([max(d["EB"], 1e-3) for d in tA if d["EB"] > 0])
            box = (10**np.percentile(ldn, 25), 10**np.percentile(ldn, 75),
                   10**np.percentile(leb, 25), 10**np.percentile(leb, 75))
            print(f"(Dn,EB) IQR box: Dn[{box[0]:.3f},{box[1]:.3f}] EB[{box[2]:.3f},{box[3]:.3f}]"
                  + (f"   (frozen Dn{FROZEN_BOX['Dn']} EB{FROZEN_BOX['EB']})" if geom == "1d" else ""))
        if geom == "1d":
            print("\n-- 1D self-validation vs frozen full medians --")
            ok = all(abs(med_pos(data, cb)[0] - FROZEN_FULL[cb]) < 0.03 for cb in CONE)
            print("   PASS" if ok else "   CHECK — differences >0.03 (likely screening 6903 vs 6118 diff)")


if __name__ == "__main__":
    main()
