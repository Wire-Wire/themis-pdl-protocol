"""RUN 13 — VALIDATION-STATUS TABLE for the 107 candidates.

User refinement: NOT "splitting 107 into weak subsets" — a transparent per-event validation-status table.
Joins the moment classifier (run9a, 107 rows) with the complete ESA spectral check (run12, 107 rows) so
every event carries BOTH its moment metrics AND its spectral-validation status. Derives one status label:

  A_SPECTRAL_CONFIRMED   : moment-clean PDL AND spectrum stays shocked-sheath   (SHEATH_CONSISTENT)
  B_MOMENT_CLEAN_SPEC_AMBIG : moment-clean PDL, spectrum borderline             (AMBIGUOUS_SPEC)
  C_SPECTRAL_FALSE_POS   : moment-clean BUT spectrum reveals keV boundary/atypical (HOT_BOUNDARY_FLAG/SHAPE_FLAG)
  D_NO_SPECTRA           : no usable spectra                                     (NO_DATA)

This is the honest re-tier the review asked for: the headline becomes "N_A spectrally-confirmed clean of
107 moment-classified", with the boundary false-positives shown explicitly.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, collections

CONF = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run9_candidates\confirmed_pdl_candidates.csv")
SPEC = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run12_spectral\spectral_metrics.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run13_validation")
os.makedirs(OUT, exist_ok=True)

TIER = {'SHEATH_CONSISTENT': 'A_SPECTRAL_CONFIRMED', 'AMBIGUOUS_SPEC': 'B_MOMENT_CLEAN_SPEC_AMBIG',
        'HOT_BOUNDARY_FLAG': 'C_SPECTRAL_FALSE_POS', 'SHAPE_FLAG': 'C_SPECTRAL_FALSE_POS',
        'NO_DATA': 'D_NO_SPECTRA'}
TIER_ORDER = ['A_SPECTRAL_CONFIRMED', 'B_MOMENT_CLEAN_SPEC_AMBIG', 'C_SPECTRAL_FALSE_POS', 'D_NO_SPECTRA']


def main():
    spec = {r['eid']: r for r in csv.DictReader(open(SPEC))}
    conf = list(csv.DictReader(open(CONF)))
    rows = []
    for c in conf:
        s = spec.get(c['eid'], {})
        st = s.get('status', 'NO_DATA')
        rows.append(dict(rank=int(c['rank']), eid=c['eid'], probe=c['probe'], date=c['date'],
                         Dn_mem=c['Dn_mem'], EB_mem=c['EB_mem'], n_near=c['n_near'], beta_near=c['beta_near'],
                         T_near=c['T_near'], bz=c['bz'], cone=c['cone'], sza=c['sza'], dp=c['dp'],
                         spec_status=st, peak_ratio=s.get('peak_ratio', ''), shape_corr=s.get('shape_corr', ''),
                         flux_ratio=s.get('flux_ratio', ''), tier=TIER.get(st, 'D_NO_SPECTRA')))
    cnt = collections.Counter(r['tier'] for r in rows)
    scnt = collections.Counter(r['spec_status'] for r in rows)

    cols = ['rank', 'eid', 'probe', 'date', 'tier', 'spec_status', 'Dn_mem', 'flux_ratio', 'peak_ratio',
            'shape_corr', 'n_near', 'beta_near', 'T_near', 'bz', 'cone', 'sza', 'dp', 'EB_mem']
    with open(os.path.join(OUT, "validation_status.csv"), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
        for r in sorted(rows, key=lambda r: (TIER_ORDER.index(r['tier']), r['rank'])):
            w.writerow(r)

    md = ["# Validation-status table — 107 moment-classified PDL candidates",
          "",
          "Join of the moment classifier (Run 9a) with the COMPLETE ESA ion-spectral check (Run 12, all 107).",
          "Each event carries its moment metrics AND its spectral-validation status. This is the honest",
          "re-tier the external review required — NOT a split into weak subsets.",
          "",
          "## Headline",
          f"- **{cnt['A_SPECTRAL_CONFIRMED']} of 107 are spectrally-confirmed clean** (moment-clean PDL AND the ion",
          "  spectrum stays shocked-sheath into the near-MP shell: peak energy unchanged, shape preserved, flux depleted).",
          f"- **{cnt['C_SPECTRAL_FALSE_POS']} are spectral false-positives** — moment-clean but the near-MP ion spectrum jumps to keV",
          "  (boundary-layer / magnetosphere) or changes shape; the moment temperature screen could not see this.",
          f"- **{cnt['B_MOMENT_CLEAN_SPEC_AMBIG']} are moment-clean but spectrally borderline**; {cnt['D_NO_SPECTRA']} have no usable spectra.",
          "",
          "→ Correct phrasing: **107 moment-classified clean candidates, of which %d are spectrally confirmed.**" % cnt['A_SPECTRAL_CONFIRMED'],
          "   The earlier '107 spectrally confirmed' was wrong (only the atlas subset had spectra) — now corrected for all 107.",
          "",
          "## Tier counts",
          "| tier | meaning | spectral status | N |",
          "|---|---|---|---|"]
    meaning = {'A_SPECTRAL_CONFIRMED': 'spectrally-confirmed clean sheath PDL', 'B_MOMENT_CLEAN_SPEC_AMBIG': 'moment-clean, spectrum borderline',
               'C_SPECTRAL_FALSE_POS': 'spectral false-positive (keV boundary / atypical)', 'D_NO_SPECTRA': 'no usable spectra'}
    statusmap = {'A_SPECTRAL_CONFIRMED': 'SHEATH_CONSISTENT', 'B_MOMENT_CLEAN_SPEC_AMBIG': 'AMBIGUOUS_SPEC',
                 'C_SPECTRAL_FALSE_POS': 'HOT_BOUNDARY_FLAG / SHAPE_FLAG', 'D_NO_SPECTRA': 'NO_DATA'}
    for t in TIER_ORDER:
        md.append(f"| {t} | {meaning[t]} | {statusmap[t]} | {cnt[t]} |")
    md.append(f"| **total** | | | **{len(rows)}** |")
    md.append("")
    md.append(f"raw spectral-status counts: {dict(scnt)}")
    md.append("")
    md.append("## Full table (sorted by tier, then moment rank)")
    md.append("| rank | eid | tier | spec | Dn_mem | flux_ratio | peak_ratio | shape_corr | n_near | beta | T_near | Bz | cone | Dp |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: (TIER_ORDER.index(r['tier']), r['rank'])):
        md.append(f"| {r['rank']} | {r['eid']} | {r['tier'][0]} | {r['spec_status']} | {r['Dn_mem']} | {r['flux_ratio']} | "
                  f"{r['peak_ratio']} | {r['shape_corr']} | {r['n_near']} | {r['beta_near']} | {r['T_near']} | {r['bz']} | {r['cone']} | {r['dp']} |")
    txt = "\n".join(md)
    with open(os.path.join(OUT, "VALIDATION_STATUS.md"), 'w', encoding='utf-8') as f:
        f.write(txt + "\n")
    print(f"tiers: {dict(cnt)}")
    print(f"spectral status: {dict(scnt)}")
    print(f"\nsaved -> {OUT}\\validation_status.csv , VALIDATION_STATUS.md")
    print(f"\nHEADLINE: {cnt['A_SPECTRAL_CONFIRMED']} of 107 spectrally-confirmed clean; "
          f"{cnt['C_SPECTRAL_FALSE_POS']} spectral false-positives; {cnt['B_MOMENT_CLEAN_SPEC_AMBIG']} borderline.")


if __name__ == "__main__":
    main()
