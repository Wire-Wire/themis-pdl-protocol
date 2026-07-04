"""RUN 19b — supplementary Appendix-E stats (NO substrate re-pass; reads the run19 per-encounter CSV).
Reports median/IQR/fractions of the pressure-budget ratios per population (for the revision: report
spread, not just medians), and checks whether the substrate carries a velocity VECTOR (for the v_normal
sensitivity) or only the speed |v|. Pure reporting on frozen data — not a new analysis branch.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, glob
import numpy as np

CSV = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run19_contrast_checks\contrast_per_encounter.csv")
SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
OUT = os.path.dirname(CSV)


def f(x):
    try:
        v = float(x); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def stats(vals):
    a = np.array([v for v in vals if np.isfinite(v)])
    if not len(a):
        return None
    return dict(n=len(a), med=np.median(a), q1=np.percentile(a, 25), q3=np.percentile(a, 75),
                f_gt125=float(np.mean(a > 1.25)), f_in=float(np.mean((a >= 0.8) & (a <= 1.25))),
                f_lt08=float(np.mean(a < 0.8)))


def main():
    rows = list(csv.DictReader(open(CSV)))
    pops = {
        'all-contributing (661)': rows,
        '60 spectrally-supported': [r for r in rows if r['cls'] == 'CONFIRMED_PDL' and r['spec'] == 'SHEATH_CONSISTENT'],
        '32 spectral false-pos': [r for r in rows if r['cls'] == 'CONFIRMED_PDL' and r['spec'] in ('HOT_BOUNDARY_FLAG', 'SHAPE_FLAG')],
    }
    out = ["RUN 19b — Appendix-E supplementary stats (median/IQR/fractions; from run19 CSV, no re-pass)", ""]
    for name, rs in pops.items():
        out.append(name)
        for key in ('R_stat', 'R_Ptot', 'R_Pdyn'):
            s = stats([f(r[key]) for r in rs])
            if s:
                out.append(f"  {key:7s} N={s['n']:3d}  median={s['med']:.3f}  IQR[{s['q1']:.3f},{s['q3']:.3f}]  "
                           f"frac>1.25={s['f_gt125']:.2f}  frac_in[0.8,1.25]={s['f_in']:.2f}  frac<0.8={s['f_lt08']:.2f}")
        out.append("")
    # velocity-vector availability (for the v_normal sensitivity)
    files = glob.glob(os.path.join(SUB, "*.npz"))
    if files:
        keys = sorted(np.load(files[0], allow_pickle=True).files)
        out.append("substrate npz keys (sample file): " + ", ".join(keys))
        has_vvec = any(k.lower() in ('vx', 'vy', 'vz', 'vx_gse', 'vy_gse', 'vz_gse', 'v_gse', 'vel', 'vgse') for k in keys)
        out.append(f"velocity VECTOR present? {has_vvec}   (if False: only speed |v|='vmag' -> v_normal NOT directly available; "
                   f"P_dyn=rho*|v|^2 is an UPPER estimate of the normal dynamic-pressure compensation)")
    txt = "\n".join(out)
    print(txt)
    with open(os.path.join(OUT, "APPENDIX_E_SUPP.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\nsaved ->", os.path.join(OUT, "APPENDIX_E_SUPP.txt"))


if __name__ == "__main__":
    main()
