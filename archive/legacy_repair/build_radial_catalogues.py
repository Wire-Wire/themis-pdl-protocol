"""Convert the dual-geometry archive output into the schema the existing aggregation +
harmonisation pipelines expect, for BOTH geometries:
  - <geom>_encounter_catalogue_extra.csv  (24-col, same as encounter_catalogue_extra.csv)
  - <geom>_crossprobe_overlap_groups_extra.csv  (same-day multi-probe overlap groups)
Also prints OOM disagreement counts (Table 4 / Appendix C input) per geometry.
Run after recompute_archive_radial.py. Outputs to option3/derived/.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import csv, os
from collections import defaultdict
import numpy as np

CAT = P(r"H:\0mssl\review\repair\option3\derived\archive_radial_catalogue.csv")
OUTDIR = P(r"H:\0mssl\review\repair\option3\derived")
CONE = ["quasi-radial", "low-cone", "intermediate", "perpendicular"]
EXTRA_COLS = ["encounter_id", "date", "year", "month", "probe", "sza_deg", "mlt", "cone_deg",
              "clock_deg", "dp_nPa", "bz_nT", "near_occ", "bg_occ", "Dn", "EB", "cone_bin_calc",
              "retained", "qc_transition_cleanliness", "qc_disturbance", "qc_boundary_motion",
              "omni_context_quality_note", "boundary_uncertainty_note", "r_re", "x_gsm_re"]
OVERLAP_COLS = ["date", "probe", "cone_deg", "dp_nPa", "Dn", "EB", "qc_transition", "boundary_note", "cone_bin"]


def F(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def derive(geom, r):
    dn = F(r.get(f"Dn_{geom}")); eb = F(r.get(f"EB_{geom}"))
    no = F(r.get(f"near_occ_{geom}")); bo = F(r.get(f"bg_occ_{geom}")); dp = F(r.get("dp_nPa"))
    qtc = "clean" if (dn is not None and eb is not None and dn < 1 and eb > 1) else \
          ("mixed" if (dn is not None and eb is not None and (dn < 1 or eb > 1)) else "unclear")
    qd = "undisturbed" if (no and bo and no > 0.1 and bo > 0.05) else "uncertain"
    qbm = "stable" if (dp and 2 <= dp <= 6) else "uncertain"
    bun = "plausible" if (dp and 2 <= dp <= 6) else "uncertain"
    return dn, eb, no, bo, qtc, qd, qbm, bun


def main():
    rows = list(csv.DictReader(open(CAT)))
    for geom in ("1d", "rad"):
        recs = []
        for r in rows:
            if str(r.get(f"evaluable_{geom}")) != "True":
                continue
            dn, eb, no, bo, qtc, qd, qbm, bun = derive(geom, r)
            recs.append({"encounter_id": r["encounter_id"], "date": r["date"], "year": r.get("year"),
                         "month": r.get("month"), "probe": r["probe"], "sza_deg": r.get("sza_deg"),
                         "mlt": "", "cone_deg": r.get("cone_deg"), "clock_deg": r.get("clock_deg"),
                         "dp_nPa": r.get("dp_nPa"), "bz_nT": r.get("bz_nT"), "near_occ": no, "bg_occ": bo,
                         "Dn": dn, "EB": eb, "cone_bin_calc": r.get("cone_bin"), "retained": True,
                         "qc_transition_cleanliness": qtc, "qc_disturbance": qd, "qc_boundary_motion": qbm,
                         "omni_context_quality_note": "good", "boundary_uncertainty_note": bun,
                         "r_re": r.get("r_re"), "x_gsm_re": r.get("x_gsm_re")})
        cat_path = os.path.join(OUTDIR, f"{geom}_encounter_catalogue_extra.csv")
        with open(cat_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=EXTRA_COLS, extrasaction="ignore"); w.writeheader()
            for e in recs:
                w.writerow(e)
        # overlap groups: same-day, >=2 probes
        byd = defaultdict(list)
        for e in recs:
            byd[e["date"]].append(e)
        ov_path = os.path.join(OUTDIR, f"{geom}_crossprobe_overlap_groups_extra.csv")
        n_groups = 0
        with open(ov_path, "w", newline="") as f:
            w = csv.writer(f); w.writerow(OVERLAP_COLS)
            for date, es in sorted(byd.items()):
                if len({e["probe"] for e in es}) < 2:
                    continue
                n_groups += 1
                for e in es:
                    w.writerow([date, e["probe"], e["cone_deg"], e["dp_nPa"], e["Dn"], e["EB"],
                                e["qc_transition_cleanliness"], e["boundary_uncertainty_note"], e["cone_bin_calc"]])
        # OOM count (LC/QR groups with >1 dex Dn spread)
        oom = 0; lcqr = 0
        for date, es in byd.items():
            if len({e["probe"] for e in es}) < 2:
                continue
            dns = [e["Dn"] for e in es if e["Dn"] and e["Dn"] > 0]
            if any(e["cone_bin_calc"] in ("low-cone", "quasi-radial") for e in es):
                lcqr += 1
                if len(dns) >= 2 and (max(np.log10(dns)) - min(np.log10(dns))) > 1.0:
                    oom += 1
        print(f"[{geom}] retained={len(recs)}  overlap_groups={n_groups}  LC/QR_groups={lcqr}  OOM_LCQR={oom}")
        print(f"        wrote {cat_path}")
        print(f"        wrote {ov_path}")


if __name__ == "__main__":
    main()
