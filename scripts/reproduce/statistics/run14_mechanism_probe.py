"""RUN 14 (supplementary) — PRE-REGISTERED deeper mechanism probe.

Run 11 showed the crude Bz-sign claim does NOT survive Dp control. But the literature PDL mechanism is
NOT Bz-sign; it is LOW MAGNETIC SHEAR (northward / small clock angle -> field piles up, depletes) behind a
QUASI-PERPENDICULAR bow shock (high cone angle). Bz-sign is too coarse a proxy. This tests the physically-
motivated regime variables on the SAME 661 contributing encounters, controlling Dp + probe.

PRE-REGISTERED hypotheses (signs fixed BEFORE fitting; response = Dn_mem, LOWER = more depleted):
  H1 magnetic shear: more northward / lower shear (cos_clock -> +1) => deeper => cos_clock coef < 0
  H2 quasi-perp shock: higher cone angle => cleaner/deeper sheath PDL    => cone coef < 0
  H3 (control) Dp: established (Run 11) Dp coef > 0 (higher Dp -> higher Dn)
DECISION RULE: a mechanism variable "survives" iff its coef has the PREDICTED sign AND p < 0.05 after
controlling Dp + probe (and, as the stronger test, within a fixed Dp band where Dp variance is removed).
If neither H1 nor H2 survives -> adopt the downgrade framing (Dp-ordered, IMF dependence inconclusive).
NOTE (user): Dp is the dominant ORDERING variable, NOT a proven PDL formation mechanism (it also sets
boundary position, sheath compression, observability, selection) — so Dp is a control here, not a claim.

Reads run10 funnel_contributing.csv (enriched with clock/cone/ma/bx). Pure-numpy OLS.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, collections
import numpy as np

SRC = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run14_mechanism")
os.makedirs(OUT, exist_ok=True)
DP_BAND = (1.5, 3.0)

try:
    from scipy import stats as _st
    def pval(t, dof): return 2.0 * _st.t.sf(abs(t), dof)
    def mwu(a, b, alt): return _st.mannwhitneyu(a, b, alternative=alt).pvalue
except Exception:
    _st = None
    def pval(t, dof):
        from math import erfc, sqrt
        return erfc(abs(t) / sqrt(2.0))
    def mwu(a, b, alt): return float('nan')


def fnum(s):
    try:
        v = float(s); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def load():
    rows = []
    for r in csv.DictReader(open(SRC)):
        d = dict(eid=r['eid'], probe=r['probe'], Dn=fnum(r['Dn_mem']),
                 clock=fnum(r['clock']), cone=fnum(r['cone']), dp=fnum(r['dp']),
                 ma=fnum(r['ma']), bx=fnum(r['bx']), sza=fnum(r['sza']),
                 bz_sign=(int(r['bz_sign']) if r['bz_sign'] != '' else -1))
        d['cos_clock'] = np.cos(np.radians(d['clock'])) if np.isfinite(d['clock']) else np.nan
        rows.append(d)
    return rows


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta; n, p = X.shape; dof = n - p
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
    r2 = 1 - float(resid @ resid) / float(((y - y.mean()) ** 2).sum())
    return beta, se, dof, r2


def model(rows, terms, o, label, predicted):
    """terms: list of (name, getter). Always adds intercept + probe dummies. Response = log(Dn)."""
    use = [r for r in rows if np.isfinite(r['Dn']) and r['Dn'] > 0 and all(np.isfinite(g(r)) for _, g in terms)]
    probes = sorted(set(r['probe'] for r in use)); ref = probes[0]
    cols = ['intercept'] + [t[0] for t in terms] + [f'probe_{p}' for p in probes[1:]]
    X = np.array([[1.0] + [g(r) for _, g in terms] + [1.0 if r['probe'] == p else 0.0 for p in probes[1:]] for r in use])
    y = np.log(np.array([r['Dn'] for r in use]))
    beta, se, dof, r2 = ols(X, y)
    o.append(f"  [{label}]  N={len(use)}  R2={r2:.3f}  (probe ref {ref})")
    for c, b, s in zip(cols, beta, se):
        if c.startswith('probe_') or c == 'intercept':
            continue
        t = b / s if s > 0 else float('nan'); p = pval(t, dof)
        pred = predicted.get(c, '')
        survived = ''
        if pred:
            ok = ((b < 0) == (pred == '<0')) and p < 0.05
            survived = '  <-- SURVIVES (predicted %s, p<0.05)' % pred if ok else ('  predicted %s; %s' % (pred, 'sign OK p>=0.05' if (b < 0) == (pred == '<0') else 'WRONG sign'))
        o.append(f"    {c:12s} coef={b:+.4f} ({100*(np.exp(b)-1):+5.1f}%/unit) SE={s:.4f} p={p:.4f}{survived}")
    return beta, cols, use


def band_test(rows, o):
    use = [r for r in rows if np.isfinite(r['dp']) and DP_BAND[0] <= r['dp'] <= DP_BAND[1]
           and r['bz_sign'] in (0, 1) and np.isfinite(r['Dn']) and r['Dn'] > 0]
    north = np.array([r['Dn'] for r in use if r['bz_sign'] == 1])
    south = np.array([r['Dn'] for r in use if r['bz_sign'] == 0])
    o.append(f"--- TEST 1: Bz-sign within FIXED Dp band {DP_BAND} (removes the Dp confound by restriction) ---")
    o.append(f"  N={len(use)} (north {len(north)}, south {len(south)});  Dp median in band = {np.median([r['dp'] for r in use]):.2f}")
    if len(north) >= 10 and len(south) >= 10:
        diff = np.median(north) - np.median(south)
        o.append(f"  median Dn north={np.median(north):.3f} vs south={np.median(south):.3f}  diff={diff:+.3f}")
        o.append(f"  Mann-Whitney north<south p={mwu(north, south, 'less'):.4f}")
        rng = np.random.default_rng(11); md = [np.median(rng.choice(north, len(north))) - np.median(rng.choice(south, len(south))) for _ in range(3000)]
        lo, hi = np.percentile(md, [2.5, 97.5])
        o.append(f"  median-diff bootstrap 95% CI [{lo:+.3f},{hi:+.3f}]  ({'excludes 0 -> Bz survives in band' if (lo<0)==(hi<0) else 'INCLUDES 0 -> Bz does NOT survive in band'})")
    else:
        o.append("  insufficient N in band")
    return use


def quasi_perp(rows, o):
    o.append("")
    o.append("--- TEST 4: quasi-perpendicular shock (cone>45) vs quasi-parallel, stratified by Dp median ---")
    use = [r for r in rows if np.isfinite(r['cone']) and np.isfinite(r['dp']) and np.isfinite(r['Dn']) and r['Dn'] > 0]
    dpm = np.median([r['dp'] for r in use])
    for lab, sub in (('Dp<=med', [r for r in use if r['dp'] <= dpm]), ('Dp>med', [r for r in use if r['dp'] > dpm])):
        qp = [r['Dn'] for r in sub if r['cone'] > 45]; qpar = [r['Dn'] for r in sub if r['cone'] <= 45]
        if qp and qpar:
            o.append(f"  {lab:8s}: quasi-perp median Dn={np.median(qp):.3f} (N={len(qp)}) vs quasi-par={np.median(qpar):.3f} (N={len(qpar)})  diff={np.median(qp)-np.median(qpar):+.3f}")


def main():
    rows = load()
    # sanity: cos_clock should track bz sign
    cc = np.array([r['cos_clock'] for r in rows if np.isfinite(r['cos_clock']) and r['bz_sign'] in (0, 1)])
    bs = np.array([r['bz_sign'] for r in rows if np.isfinite(r['cos_clock']) and r['bz_sign'] in (0, 1)])
    o = ["RUN 14 — PRE-REGISTERED mechanism probe (magnetic shear / quasi-perp shock), controlling Dp+probe",
         "response Dn_mem (LOWER=more depleted); H1 cos_clock coef<0 (low shear deeper), H2 cone coef<0 (quasi-perp deeper)",
         f"sanity: corr(cos_clock, Bz>0) = {np.corrcoef(cc, bs)[0,1]:+.3f} (should be +, confirms cos_clock = northward-ness)",
         f"contributing N={len(rows)}; with finite clock={sum(np.isfinite(r['cos_clock']) for r in rows)}, cone={sum(np.isfinite(r['cone']) for r in rows)}, ma={sum(np.isfinite(r['ma']) for r in rows)}",
         ""]
    band_test(rows, o)
    o.append("")
    o.append("--- TEST 2: continuous regime variables, log(Dn) ~ var + Dp + C(probe) ---")
    pred = {'cos_clock': '<0', 'cone': '<0', 'dp': '>0', 'ma': '', 'bx': ''}
    model(rows, [('cos_clock', lambda r: r['cos_clock']), ('dp', lambda r: r['dp'])], o, 'shear: cos_clock+Dp', pred)
    model(rows, [('cone', lambda r: r['cone']), ('dp', lambda r: r['dp'])], o, 'shock: cone+Dp', pred)
    model(rows, [('ma', lambda r: r['ma']), ('dp', lambda r: r['dp'])], o, 'Mach: ma+Dp', pred)
    model(rows, [('cos_clock', lambda r: r['cos_clock']), ('cone', lambda r: r['cone']),
                 ('ma', lambda r: r['ma']), ('dp', lambda r: r['dp'])], o, 'FULL: cos_clock+cone+ma+Dp', pred)
    o.append("")
    o.append(f"--- TEST 3: regime variables WITHIN fixed Dp band {DP_BAND} (Dp~const; cleanest mechanism test) ---")
    band = [r for r in rows if np.isfinite(r['dp']) and DP_BAND[0] <= r['dp'] <= DP_BAND[1]]
    model(band, [('cos_clock', lambda r: r['cos_clock']), ('cone', lambda r: r['cone'])], o, 'in-band: cos_clock+cone', pred)
    quasi_perp(rows, o)
    o.append("")
    o.append("VERDICT: see which (if any) of H1/H2 SURVIVES. If none -> downgrade framing:")
    o.append("  'the near-MP depletion is ordered by dynamic pressure; a low-shear/quasi-perp IMF dependence")
    o.append("   is in the predicted direction but not statistically resolved in this selection-limited sample.'")
    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "MECHANISM_PROBE.txt"), "w") as f:
        f.write(txt + "\n")
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
