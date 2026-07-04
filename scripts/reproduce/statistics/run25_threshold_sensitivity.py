"""run25 — membership density-threshold sensitivity (a reviewer threshold question).

The reviewer asked: the membership density floor (0.3 cm^-3) may be too low; how does the
near-magnetopause contrast change if it is raised? Reuses the EXACT committed membership +
contrast functions from run_hardening_checks (run7), sweeps the density floor over
{0.1, 0.3(baseline), 0.5, 1.0, 2.0} cm^-3 (and a few combined density+temperature settings),
and reports Dn_mem (EB_mem) for the headline near-MP shell s in [0.05,0.20). Does NOT touch run7.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os
import numpy as np
from run_hardening_checks import load_enc, member_mask, shell_dn, KEYSHELLS, MIN_NENC
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run25_threshold_sensitivity")
os.makedirs(OUT, exist_ok=True)


def sweep(data, settings):
    rows = []
    for nm, tb, vf, nf in settings:
        cells = []
        for sh, lo, hi in KEYSHELLS:
            dns, ebs = [], []
            for e in data:
                dn, eb = shell_dn(e, member_mask(e, tb, vf, nf), lo, hi)
                if dn is not None:
                    dns.append(dn); ebs.append(eb)
            cells.append(f"{np.median(dns):.3f}({np.median(ebs):.2f})N{len(dns)}"
                         if len(dns) >= MIN_NENC else f"(N{len(dns)}<30)")
        rows.append(f"{nm:24s} " + " ".join(f"{c:>20s}" for c in cells))
    return rows


def main():
    data = [e for e in pmap(load_enc) if e is not None]
    out = [f"loaded {len(data)} encounters", "",
           "=== MEMBERSHIP DENSITY-THRESHOLD SENSITIVITY  — Dn_mem (EB_mem) N per shell ===",
           "baseline screen = TBAND 3 (T within x3 of background), VFRAC 0.2 (v>0.2 v_bg), n>floor",
           f"{'setting':24s} " + " ".join(f"{nm:>20s}" for nm, _, _ in KEYSHELLS)]
    # (A) density floor sweep
    out += sweep(data, [('n>0.1', 3.0, 0.2, 0.1), ('n>0.3 (baseline)', 3.0, 0.2, 0.3),
                        ('n>0.5', 3.0, 0.2, 0.5), ('n>1.0', 3.0, 0.2, 1.0), ('n>2.0', 3.0, 0.2, 2.0)])
    # (B) combined density + tighter temperature band
    out += ["", "--- combined: raise density floor AND tighten temperature band ---"]
    out += sweep(data, [('n>1.0 & TBAND=2', 2.0, 0.2, 1.0), ('n>2.0 & TBAND=2', 2.0, 0.2, 2.0)])
    out += ["",
            "Reading: the headline shell is the first column [0.05,0.20). If Dn_mem stays ~0.88-0.94",
            "as the density floor rises 0.1 -> 2.0, the near-MP contrast is NOT an artefact of an",
            "over-permissive density threshold."]
    txt = "\n".join(out)
    print(txt, flush=True)
    with open(os.path.join(OUT, "THRESHOLD_SENSITIVITY.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved -> run25_threshold_sensitivity/THRESHOLD_SENSITIVITY.txt", flush=True)


if __name__ == "__main__":
    main()
