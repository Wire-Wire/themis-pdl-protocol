"""RUN 15 — does the near-MP density-and-field contrast survive in the
SPECTRALLY-SUPPORTED subset? 3-population comparison from the per-encounter PDL-shell metrics.

Populations (all use the SAME per-encounter near-MP-shell s in [0.05,0.20) membership metrics):
  P1 all contributing membership-screened encounters
  P2 moment-classified clean candidates
  P3 spectrally-supported sheath-like candidates (P2 with SHEATH_CONSISTENT ion spectrum)
Per population: N; median+IQR of Dn, EB, beta_near, Dp; median SZA; median cone; probe distribution;
N with usable spectra. Question: is the Dn<1 / EB>1 contrast carried by moment-only candidates, or
does it persist in the spectrally-supported subset? Reads run10 + run12 CSVs (no substrate re-pass).
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, collections
import numpy as np

CONTRIB = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
SPEC = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run12_spectral\spectral_metrics.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run15_subset_profile")
os.makedirs(OUT, exist_ok=True)


def fnum(s):
    try:
        v = float(s); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def mi(vals):
    a = np.array([v for v in vals if np.isfinite(v)])
    if not len(a):
        return "n/a"
    return f"{np.median(a):.3f} [{np.percentile(a,25):.3f},{np.percentile(a,75):.3f}]"


def med(vals):
    a = np.array([v for v in vals if np.isfinite(v)])
    return f"{np.median(a):.1f}" if len(a) else "n/a"


def main():
    rows = list(csv.DictReader(open(CONTRIB)))
    for r in rows:
        for k in ('Dn_mem', 'EB_mem', 'beta_near', 'dp', 'sza', 'cone'):
            r[k] = fnum(r[k])
    spec = {r['eid']: r['status'] for r in csv.DictReader(open(SPEC))}
    sheath = {eid for eid, st in spec.items() if st == 'SHEATH_CONSISTENT'}
    usable = {eid for eid, st in spec.items() if st != 'NO_DATA'}

    p1 = rows
    p2 = [r for r in rows if r['cls'] == 'CONFIRMED_PDL']
    p3 = [r for r in p2 if r['eid'] in sheath]
    pops = [("all contributing (membership-screened)", p1),
            ("moment-classified clean candidates", p2),
            ("spectrally-supported sheath-like candidates", p3)]

    lines = ["RUN 15 — spectrally-supported subset: does the near-MP contrast survive?",
             "(near-MP-shell s in [0.05,0.20) membership metrics, per encounter; Dn<1 = depletion, EB>1 = field pile-up)",
             ""]
    hdr = f"{'population':46s}{'N':>5s}  {'Dn (med[IQR])':>22s}  {'EB (med[IQR])':>22s}  {'beta_near':>20s}  {'Dp':>20s}  {'SZA':>5s}  {'cone':>5s}  {'N_spec':>6s}"
    lines.append(hdr)
    md = ["# Spectrally-supported subset profile — 3-population comparison", "",
          "Near-magnetopause-shell (s∈[0.05,0.20)) density-and-field contrast per encounter; Dn<1 = depletion, EB>1 = field pile-up.",
          "",
          "| population | N | Dn median [IQR] | EB median [IQR] | β_near median [IQR] | Dp median [IQR] | SZA med | cone med | probe distribution | N w/ usable spectra |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for name, pop in pops:
        probe = collections.Counter(r['probe'] for r in pop)
        probe_s = " ".join(f"{k}:{probe[k]}" for k in sorted(probe))
        nspec = sum(1 for r in pop if r['eid'] in usable)
        lines.append(f"{name:46s}{len(pop):5d}  {mi([r['Dn_mem'] for r in pop]):>22s}  {mi([r['EB_mem'] for r in pop]):>22s}  "
                     f"{mi([r['beta_near'] for r in pop]):>20s}  {mi([r['dp'] for r in pop]):>20s}  {med([r['sza'] for r in pop]):>5s}  "
                     f"{med([r['cone'] for r in pop]):>5s}  {nspec:>6d}")
        md.append(f"| {name} | {len(pop)} | {mi([r['Dn_mem'] for r in pop])} | {mi([r['EB_mem'] for r in pop])} | "
                  f"{mi([r['beta_near'] for r in pop])} | {mi([r['dp'] for r in pop])} | {med([r['sza'] for r in pop])} | "
                  f"{med([r['cone'] for r in pop])} | {probe_s} | {nspec} |")

    # quick verdict
    d3 = np.median([r['Dn_mem'] for r in p3]); e3 = np.median([r['EB_mem'] for r in p3])
    d2 = np.median([r['Dn_mem'] for r in p2]); e2 = np.median([r['EB_mem'] for r in p2])
    verdict = (f"\nVERDICT: spectrally-supported subset Dn={d3:.3f} (<1 => depletion {'RETAINED' if d3<1 else 'LOST'}), "
               f"EB={e3:.3f} (>1 => pile-up {'RETAINED' if e3>1 else 'LOST'}). "
               f"Moment-classified Dn={d2:.3f}, EB={e2:.3f}. "
               f"{'Contrast persists in the spectrally-supported subset (not carried by moment-only candidates).' if (d3<1 and e3>1) else 'Contrast weakens in the spectrally-supported subset.'}")
    note = ("\nNote: spectra were fetched only for the candidate subset; for P1 'N w/ usable spectra' reflects that overlap, "
            "not a spectral screen of all 661. P3 is within P2 by construction (spectral check applied to the moment-classified candidates).")
    lines.append(verdict); lines.append(note)
    md.append(verdict.replace("\n", "\n\n")); md.append(note)
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "SUBSET_PROFILE.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(OUT, "SUBSET_PROFILE.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\nsaved ->", OUT)


if __name__ == "__main__":
    main()
