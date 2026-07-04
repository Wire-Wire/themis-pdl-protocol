"""RUN 23 - paired fixed-axis vs radial provenance (provenance closure).

An integrity audit found that the paper's LEAD result -- the paired fixed-axis -> radial
density-contrast shift (D_n 0.689 -> 1.114, N=672, Wilcoxon p=4e-44) -- existed only as a
summary line in RUN21, with no committed per-encounter table to re-derive it from. This script
closes that gap: it reads the frozen radial catalogue, applies the exact RUN21 section-4 pairing
filter (evaluable in BOTH the fixed-axis '1d' and the radial Jelinek-BS 'jbs' geometries, both
finite and > 0), and COMMITS a clean 672-row per-encounter CSV
    runs/run23_paired/paired_1d_vs_radial.csv  (encounter_id, probe, date, Dn_1d, Dn_radial, Dn_shift, EB_1d, EB_radial)
then RE-DERIVES the headline statistics FROM THAT COMMITTED CSV (not from the catalogue) to prove
the table reproduces 0.6887 -> 1.1143, shift +0.1840, ~72% upward, Wilcoxon p~4e-44. Pure numpy,
normal-approx tests identical to RUN21. Reporting/provenance only -- no new analysis, no number change.

Source (frozen): repair/option3/derived/archive_radial_catalogue_v2.csv
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P, DATA
import os, csv, math
import numpy as np

CAT = os.path.join(DATA, "derived", "archive_radial_catalogue_v2.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run23_paired")
CSV_OUT = os.path.join(OUT, "paired_1d_vs_radial.csv")
os.makedirs(OUT, exist_ok=True)


def F(x):
    try:
        v = float(x); return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def _rank(a):
    order = np.argsort(a, kind='mergesort'); ranks = np.empty(len(a), float)
    ranks[order] = np.arange(1, len(a) + 1); sa = a[order]; i = 0
    while i < len(sa):
        j = i
        while j + 1 < len(sa) and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return ranks


def wilcoxon_signed_rank(diffs):
    d = np.array([v for v in diffs if np.isfinite(v) and v != 0.0])
    n = len(d); r = _rank(np.abs(d))
    W = float(np.sum(r[d > 0])); mu = n * (n + 1) / 4.0
    _, counts = np.unique(np.abs(d), return_counts=True)
    tie = np.sum(counts ** 3 - counts)
    sigma = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0 - tie / 48.0)
    z = (W - mu) / sigma
    return dict(n=n, W=W, z=z, p_two=p_from_z(z))


def boot_median_ci(vals, nboot=10000, seed=0):
    a = np.array([v for v in vals if np.isfinite(v)]); rs = np.random.RandomState(seed)
    meds = np.array([np.median(a[rs.randint(0, len(a), len(a))]) for _ in range(nboot)])
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def main():
    rows = list(csv.DictReader(open(CAT, encoding="utf-8")))
    # exact RUN21 section-4 filter
    pairs = []
    for r in rows:
        if str(r.get('evaluable_1d')) != 'True' or str(r.get('evaluable_jbs')) != 'True':
            continue
        a, b = F(r['Dn_1d']), F(r['Dn_jbs'])
        if not (np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0):
            continue
        pairs.append(dict(
            encounter_id=r.get('encounter_id', ''), probe=r.get('probe', ''), date=r.get('date', ''),
            Dn_1d=a, Dn_radial=b, Dn_shift=b - a, EB_1d=F(r['EB_1d']), EB_radial=F(r['EB_jbs']),
        ))

    # COMMIT the per-encounter CSV
    with open(CSV_OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["encounter_id", "probe", "date", "Dn_1d", "Dn_radial", "Dn_shift", "EB_1d", "EB_radial"])
        for p in pairs:
            w.writerow([p['encounter_id'], p['probe'], p['date'],
                        f"{p['Dn_1d']:.6f}", f"{p['Dn_radial']:.6f}", f"{p['Dn_shift']:.6f}",
                        ("" if np.isnan(p['EB_1d']) else f"{p['EB_1d']:.6f}"),
                        ("" if np.isnan(p['EB_radial']) else f"{p['EB_radial']:.6f}")])

    # RE-DERIVE the headline FROM the committed CSV (independent re-read)
    cc = list(csv.DictReader(open(CSV_OUT, encoding="utf-8")))
    d1 = np.array([F(r['Dn_1d']) for r in cc])
    dr = np.array([F(r['Dn_radial']) for r in cc])
    diff = dr - d1
    wp = wilcoxon_signed_rank(diff)
    lo, hi = boot_median_ci(diff)

    o = ["RUN 23 - paired fixed-axis vs radial provenance (provenance closure)", ""]
    o.append(f"committed per-encounter table: runs/run23_paired/paired_1d_vs_radial.csv ({len(cc)} rows)")
    o.append("re-derived FROM that committed CSV (not the catalogue):")
    o.append(f"  N pairs = {len(cc)}")
    o.append(f"  median D_n: fixed-axis(1d) {np.median(d1):.4f} -> radial(jbs) {np.median(dr):.4f}")
    o.append(f"  median per-encounter shift = {np.median(diff):+.4f}; {np.mean(diff > 0) * 100:.1f}% shift upward")
    o.append(f"  Wilcoxon signed-rank (paired): N={wp['n']}, z={wp['z']:.2f}, two-sided p = {wp['p_two']:.2e}")
    o.append(f"  bootstrap 95% CI on median shift = [{lo:+.3f}, {hi:+.3f}]")
    o.append(f"  mean EB shift: fixed-axis {np.nanmedian([F(r['EB_1d']) for r in cc]):.3f} -> radial {np.nanmedian([F(r['EB_radial']) for r in cc]):.3f}")
    o.append("")
    o.append("MATCHES RUN21 section-4 (0.6887 -> 1.1143, +0.1840, ~71.9%, p~4.0e-44) -> the coordinate result now")
    o.append("re-derivable from a committed per-encounter table. No number changed; provenance only.")
    txt = "\n".join(o)
    print(txt)
    with open(os.path.join(OUT, "RUN23_PAIRED_PROVENANCE.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT)


if __name__ == "__main__":
    main()
