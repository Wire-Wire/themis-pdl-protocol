"""One-command contrast analysis: population profile, near-shell contrast statistics,
and the paired coordinate-construction check, over a directory of encounters.

    pdl-contrast                       # scan data/substrate/
    pdl-contrast /path --paired-check  # add the fixed-axis vs radial comparison
    pdl-contrast --csv profile.csv

or from Python:

    from pdl_protocol import shell_profile, population_contrast, paired_coordinate_check

`shell_profile` is the population-level radial profile (dissertation Sec. 5, Fig 3);
`population_contrast` is the near-shell population statistic (D_n, E_B medians + IQR,
Sec. 9.1); `paired_coordinate_check` reproduces the Sec. 4 paired design on YOUR data:
the same encounters scored under the radial construction and under the fixed
Sun-Earth-axis construction (identical boundary models, only the measurement
direction differs), to show how much of an apparent near-boundary depletion is
manufactured by the coordinate choice alone.

Scope note: the fixed-axis construction exists here only as the comparator in the
paired check; it is the geometrically inconsistent choice the dissertation shows to
manufacture an artefact, and is deliberately not part of the analysis API.
"""
import argparse
import csv

import numpy as np

from .core import load_encounter, member_mask
from .psub import pmap, list_files

DEFAULT_SHELLS = [(0.00, 0.05), (0.05, 0.10), (0.10, 0.15), (0.15, 0.20),
                  (0.20, 0.30), (0.30, 0.40), (0.40, 0.60), (0.60, 0.80), (0.80, 1.00)]
ANALYSIS_SHELL = (0.05, 0.20)
MIN_SAMP = 3


def _per_encounter_shells(d, shells):
    e = load_encounter(d)
    if e is None:
        return None
    m = member_mask(e)
    out = []
    for lo, hi in shells:
        sel = m & (e["s"] >= lo) & (e["s"] < hi)
        if sel.sum() >= MIN_SAMP:
            out.append((lo, hi,
                        float(np.median(e["n"][sel]) / e["n_bg"]),
                        float(np.median(e["b"][sel]) / e["b_bg"]),
                        float(np.nanmedian(e["beta"][sel]))))
    return out or None


def shell_profile(substrate_dir=None, shells=DEFAULT_SHELLS):
    """Population radial profile: per shell, the across-encounter median D_n, E_B,
    beta and the number of contributing encounters."""
    files = list_files(substrate_dir) if substrate_dir else list_files()
    per_enc = pmap(lambda d: _per_encounter_shells(d, shells), files=files)
    prof = []
    for lo, hi in shells:
        dn = [v[2] for enc in per_enc for v in enc if v[0] == lo]
        eb = [v[3] for enc in per_enc for v in enc if v[0] == lo]
        bt = [v[4] for enc in per_enc for v in enc if v[0] == lo]
        if dn:
            prof.append(dict(s_lo=lo, s_hi=hi, N=len(dn),
                             D_n=round(float(np.median(dn)), 3),
                             E_B=round(float(np.median(eb)), 3),
                             beta=round(float(np.nanmedian(bt)), 3)))
    return prof


def population_contrast(substrate_dir=None, shell=ANALYSIS_SHELL):
    """Near-shell population statistics over the analysis shell (default [0.05, 0.20))."""
    files = list_files(substrate_dir) if substrate_dir else list_files()

    def one(d):
        e = load_encounter(d)
        if e is None:
            return None
        sel = member_mask(e) & (e["s"] >= shell[0]) & (e["s"] < shell[1])
        if sel.sum() < MIN_SAMP:
            return None
        return (float(np.median(e["n"][sel]) / e["n_bg"]),
                float(np.median(e["b"][sel]) / e["b_bg"]))

    pairs = pmap(one, files=files)
    if not pairs:
        return None
    dn = np.array([p[0] for p in pairs]); eb = np.array([p[1] for p in pairs])
    q = lambda x: [round(float(v), 3) for v in np.percentile(x, [25, 50, 75])]
    return dict(N=len(pairs),
                D_n_median=q(dn)[1], D_n_IQR=[q(dn)[0], q(dn)[2]],
                E_B_median=q(eb)[1], E_B_IQR=[q(eb)[0], q(eb)[2]],
                frac_depleted=round(float(np.mean(dn < 1.0)), 3))


def paired_coordinate_check(substrate_dir=None, near=(0.20, 0.40), bg=(0.60, 1.00)):
    """The Sec.-4 paired design on your data: broad-bin D_n under the radial vs the
    fixed Sun-Earth-axis construction, on the encounters evaluable under BOTH.

    Fixed-axis: the spacecraft is projected onto the Sun-Earth line and the boundary
    intercepts taken there, s_1d = (x - r_MP0) / (r_BS0 - r_MP0). No membership
    screen is applied (the minimal Sec.-4 design isolates the coordinate effect)."""
    files = list_files(substrate_dir) if substrate_dir else list_files()

    def one(d):
        e = load_encounter(d)          # radial geometric-sheath samples
        if e is None:
            return None
        mp0 = float(d["mp0"]); bs0 = float(d["bs0"])
        if not (np.isfinite(bs0) and bs0 > mp0):
            return None
        x = d["x_re"].astype(float); n = d["n"].astype(float)
        s1 = (x - mp0) / (bs0 - mp0)
        ok = np.isfinite(n) & (n > 0) & np.isfinite(s1) & (s1 > 0) & (s1 < 1)

        def ratio(s_arr, n_arr, valid):
            nb = valid & (s_arr >= near[0]) & (s_arr < near[1])
            bb = valid & (s_arr >= bg[0]) & (s_arr <= bg[1])
            if nb.sum() < MIN_SAMP or bb.sum() < MIN_SAMP:
                return None
            return float(np.median(n_arr[nb]) / np.median(n_arr[bb]))

        d1 = ratio(s1, n, ok)
        dr = ratio(e["s"], e["n"], np.ones(len(e["s"]), bool))
        if d1 is None or dr is None:
            return None
        return (d1, dr)

    pairs = pmap(one, files=files)
    if len(pairs) < 3:
        return dict(N=len(pairs), note="too few encounters evaluable under both constructions")
    d1 = np.array([p[0] for p in pairs]); dr = np.array([p[1] for p in pairs])
    res = dict(N=len(pairs),
               D_n_fixed_axis=round(float(np.median(d1)), 3),
               D_n_radial=round(float(np.median(dr)), 3),
               median_paired_shift=round(float(np.median(dr - d1)), 3),
               frac_shift_up=round(float(np.mean(dr > d1)), 3))
    try:
        from scipy import stats
        w = stats.wilcoxon(dr - d1, alternative="two-sided")
        res["wilcoxon_p"] = float(w.pvalue)
    except Exception:
        pass
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="pdl-contrast",
        description="Population contrast analysis over a directory of encounter NPZ "
                    "files: radial shell profile, near-shell statistics, and "
                    "(optionally) the paired fixed-axis vs radial coordinate check.")
    ap.add_argument("substrate", nargs="?", default=None,
                    help="directory of encounter NPZ files (default: data/substrate/)")
    ap.add_argument("--paired-check", action="store_true",
                    help="add the Sec.-4 fixed-axis vs radial paired comparison")
    ap.add_argument("--csv", default=None, help="write the shell profile to this CSV")
    a = ap.parse_args(argv)

    prof = shell_profile(a.substrate)
    if not prof:
        print("no analysable encounters found (see docs/DATA.md for the NPZ schema)")
        return 1
    print("radial shell profile (across-encounter medians):")
    print("%-14s %4s %7s %7s %7s" % ("s shell", "N", "D_n", "E_B", "beta"))
    for r in prof:
        print("[%.2f, %.2f)  %4d %7.3f %7.3f %7.3f" % (r["s_lo"], r["s_hi"], r["N"], r["D_n"], r["E_B"], r["beta"]))
    pop = population_contrast(a.substrate)
    if pop:
        print("\nnear-shell population statistic, s in [%.2f, %.2f):" % ANALYSIS_SHELL)
        print("  N = %d encounters | D_n = %.3f IQR %s | E_B = %.3f IQR %s | depleted fraction %.0f%%"
              % (pop["N"], pop["D_n_median"], pop["D_n_IQR"], pop["E_B_median"], pop["E_B_IQR"],
                 100 * pop["frac_depleted"]))
    if a.paired_check:
        pc = paired_coordinate_check(a.substrate)
        print("\npaired coordinate-construction check (broad bin %s vs background):" % "[0.20, 0.40)")
        for k, v in pc.items():
            print("  %s = %s" % (k, v))
        print("  (a large fixed-axis-to-radial shift means the coordinate choice, not the plasma,")
        print("   is producing the apparent near-boundary depletion in your data)")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(prof[0].keys()))
            w.writeheader(); w.writerows(prof)
        print("\nprofile written to", a.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
