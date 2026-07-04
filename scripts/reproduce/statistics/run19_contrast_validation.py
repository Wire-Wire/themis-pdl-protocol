"""RUN 19 — CONTRAST VALIDATION (review-prompted). Two substrate checks.

Review raised a critical objection to the contrast result (the near-MP contrast
D_n~0.68 / E_B~2.4): (i) the headline 0.68 is conditional/circular (candidates are pre-selected to
0.40<=D_n<=0.90), and (ii) E_B~2.4 may be a generic background-referenced inward field gradient, a
pressure-budget mismatch, or an ESA near-boundary density underestimate, rather than genuine pile-up.

This script answers both, on the IDENTICAL geometry/membership/classify as run10 (NO new selection):

CHECK 1 — per-encounter pressure budget. For all-contributing(661) / moment-candidates(107) /
  spectrally-supported(60) / spectral-false-pos(32), the per-encounter NEAR/BACKGROUND median ratios of
    P_B   = B^2/2mu0                 (magnetic)
    P_th  = p_th                     (ION thermal only -> the static budget below is an
                                      ION-THERMAL + MAGNETIC budget, NOT full total pressure)
    P_stat= P_B + P_th               (ion-thermal + magnetic static budget)
    P_dyn = rho v^2 (ion)            (dynamic; reported so a static excess can be tested against
                                      flow stagnation, which is the pile-up mechanism)
    P_tot = P_B + P_th + P_dyn       (ion total; omits electron thermal ~10-20%)
  plus D_n, E_B, v near/bg, beta_near, T_near. A genuine pile-up either (a) conserves the static
  budget (P_stat ~ 1) or (b) shows a static excess balanced by a dynamic-pressure drop (P_tot ~ 1,
  P_dyn < 1) = stagnation. A static AND total over-pressure with no dynamic compensation points to an
  ESA density-underestimate / background-referencing artefact.

CHECK 2 — event-wise density-field coupling. The defining PDL signature is the n-|B| ANTICORRELATION.
  (a) event-level Spearman(D_n, E_B) across each population (negative => deeper depletion goes with
      stronger field, as a PDL requires); computed on the UNRESTRICTED 661 (where D_n is not truncated)
      and on the 60. (b) per-event Spearman(n, |B|) over the inner membership samples s in [0.05,0.40),
      with the distribution (median r, frac<0) compared across 60 vs 32 vs all-contributing.

Read-out rule (user): if the 60 show reasonable pressure compensation AND the D_n-E_B anticorrelation
holds -> keep "selected candidates have D_n~0.68/E_B~2.4" as a CONDITIONAL physical contrast while the
main text also reports the unconditioned all-contributing D_n (~0.92). Otherwise -> downgrade the contrast result to
the unconditioned screened profile + a field enhancement that is a robust but not uniquely-coupled
inward gradient, with 0.68 as a selected-candidate descriptive statistic only.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, collections
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap, list_files

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run19_contrast_checks")
SPEC = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run12_spectral\spectral_metrics.csv")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0); PDL = (0.05, 0.20); INNER = (0.05, 0.40); MIN_BG = 5
KB = 1.602e-4          # nPa*cm3/eV   : T[eV] = p_th/(n*KB)
KMAG = 3.97887e-4      # nPa/nT^2     : P_B = KMAG*B^2
KDYN = 1.6726e-6       # nPa/(cm-3*(km/s)^2) : P_dyn = KDYN*n*v^2
N_FLOOR = 0.3; TBAND = 3.0; VFRAC = 0.2; MIN_SHELL = 5


def classify(r):
    """IDENTICAL to run10/run9a.classify()."""
    if r['beta_near'] < 0.10 or r['n_near'] < 1.0 or r['T_near'] > 1200:
        return 'BOUNDARY_LAYER'
    if (np.isfinite(r['bz']) and r['bz'] > 0 and 0.40 <= r['Dn_mem'] <= 0.90 and r['n_near'] >= 2.0
            and 0.10 <= r['beta_near'] <= 2.0 and 50 <= r['T_near'] <= 800 and r['EB_mem'] > 1.2 and r['n_mem'] >= 15):
        return 'CONFIRMED_PDL'
    return 'AMBIGUOUS'


def _rank(a):
    a = np.asarray(a, float)
    order = np.argsort(a, kind='mergesort')
    ranks = np.empty(len(a), float); ranks[order] = np.arange(1, len(a) + 1)
    sa = a[order]; i = 0
    while i < len(sa):
        j = i
        while j + 1 < len(sa) and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return ranks


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, int(m.sum())
    rx = _rank(x[m]); ry = _rank(y[m])
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, int(m.sum())
    return float(np.corrcoef(rx, ry)[0, 1]), int(m.sum())


def perm_p(x, y, obs, niter=10000, seed=0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y); x = x[m]; y = y[m]
    if len(x) < 4 or not np.isfinite(obs):
        return np.nan
    rx = _rank(x); ry = _rank(y); rs = np.random.RandomState(seed); cnt = 0
    for _ in range(niter):
        r = np.corrcoef(rx, rs.permutation(ry))[0, 1]
        if abs(r) >= abs(obs) - 1e-12:
            cnt += 1
    return (cnt + 1) / (niter + 1)


def funnel_one(name, d):
    eid = name[:-4]; probe = eid.rsplit('_', 1)[-1]
    base = dict(eid=eid, probe=probe, stage='x')
    t = d['t'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if len(t) < 10 or not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return {**base, 'stage': 'no_params'}
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float); beta = d['beta'].astype(float)
    pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    bz = float(d['bz']) if np.isfinite(d['bz']) else np.nan
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha; r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    geo = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[geo]; nv = n[geo]; bv = b[geo]; Tv = T[geo]; vv = v[geo]; betav = beta[geo]; pthv = pth[geo]
    if len(sv) < MIN_BG + 3:
        return {**base, 'stage': 'no_sheath', 'bz': bz}
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return {**base, 'stage': 'no_bg', 'bz': bz}
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm]); T_bg = np.nanmedian(Tv[bgm]); v_bg = np.nanmedian(vv[bgm])
    if not (n_bg > 0 and b_bg > 0 and np.isfinite(T_bg) and T_bg > 0):
        return {**base, 'stage': 'no_bg', 'bz': bz}
    mem = (nv > N_FLOOR) & np.isfinite(Tv) & (Tv > T_bg / TBAND) & (Tv < T_bg * TBAND)
    if np.isfinite(v_bg) and v_bg > 0:
        mem = mem & np.isfinite(vv) & (vv > VFRAC * v_bg)
    pm = mem & (sv >= PDL[0]) & (sv < PDL[1])
    if pm.sum() < MIN_SHELL:
        return {**base, 'stage': 'no_shell', 'bz': bz}

    # --- pressure terms per sample (nPa) ---
    PB = KMAG * bv * bv
    Pth = pthv
    Pdyn = KDYN * nv * vv * vv
    Pstat = PB + Pth
    Ptot = Pstat + Pdyn

    def ratio(arr):
        an = np.nanmedian(arr[pm]); ab = np.nanmedian(arr[bgm])
        return float(an / ab) if (np.isfinite(an) and np.isfinite(ab) and ab > 0) else np.nan

    # beta cross-check: KMAG*B^2 should equal p_th/beta if substrate beta = P_th/P_B
    with np.errstate(invalid='ignore', divide='ignore'):
        pb_from_beta = np.where(betav > 0, pthv / betav, np.nan)
        xc = np.nanmedian((PB / pb_from_beta)[bgm])

    inner = mem & (sv >= INNER[0]) & (sv < INNER[1])
    r_nB, n_in = (spearman(nv[inner], bv[inner]) if inner.sum() >= 5 else (np.nan, int(inner.sum())))

    rec = dict(eid=eid, probe=probe, stage='contrib', bz=bz,
               Dn_mem=float(np.median(nv[pm]) / n_bg), EB_mem=float(np.median(bv[pm]) / b_bg),
               n_near=float(np.median(nv[pm])), beta_near=float(np.nanmedian(betav[pm])),
               T_near=float(np.nanmedian(Tv[pm])), n_mem=int(pm.sum()),
               R_PB=ratio(PB), R_Pth=ratio(Pth), R_stat=ratio(Pstat), R_Pdyn=ratio(Pdyn),
               R_Ptot=ratio(Ptot), R_v=ratio(vv),
               xcheck=float(xc) if np.isfinite(xc) else np.nan, r_nB=r_nB, n_in=n_in)
    rec['cls'] = classify(rec)
    return rec


def mi(vals, fmt="{:.3f}"):
    a = np.array([v for v in vals if np.isfinite(v)])
    if not len(a):
        return "n/a"
    return f"{fmt.format(np.median(a))} [{fmt.format(np.percentile(a,25))},{fmt.format(np.percentile(a,75))}]"


def main():
    files = list_files()
    rows = pmap(funnel_one, files=files, with_name=True)
    contrib = [r for r in rows if r.get('stage') == 'contrib']
    spec = {r['eid']: r['status'] for r in csv.DictReader(open(SPEC))}
    for r in contrib:
        r['spec'] = spec.get(r['eid'], '')

    P_all = contrib
    P_cand = [r for r in contrib if r['cls'] == 'CONFIRMED_PDL']
    P_60 = [r for r in P_cand if r['spec'] == 'SHEATH_CONSISTENT']
    P_32 = [r for r in P_cand if r['spec'] in ('HOT_BOUNDARY_FLAG', 'SHAPE_FLAG')]
    pops = [("all contributing", P_all), ("moment candidates (107)", P_cand),
            ("spectrally-supported (60)", P_60), ("spectral false-pos (32)", P_32)]

    o = []
    o.append("RUN 19 — CONTRAST VALIDATION (pressure budget + event-wise density-field coupling)")
    o.append(f"populations: all-contributing={len(P_all)}  moment-cand={len(P_cand)}  "
             f"spectrally-supported={len(P_60)}  spectral-false-pos={len(P_32)}")
    xc_all = np.nanmedian([r['xcheck'] for r in P_all if np.isfinite(r['xcheck'])])
    o.append(f"unit cross-check (median over backgrounds of KMAG*B^2 / (p_th/beta)) = {xc_all:.3f}  "
             f"(should be ~1 if substrate beta = P_th/P_B; ratios are near/bg so a constant offset cancels)")
    o.append("")
    o.append("=== CHECK 1: per-encounter NEAR/BACKGROUND pressure ratios (median [IQR]) ===")
    o.append("P_th = ION thermal only; 'static' = ion-thermal+magnetic (NOT full total pressure);")
    o.append("P_dyn, P_tot are ION terms (omit electron thermal ~10-20%). A genuine pile-up: P_stat~1, OR")
    o.append("P_stat>1 compensated by P_dyn<1 so P_tot~1 (stagnation). P_stat AND P_tot >>1 => artefact/under-count.")
    o.append("")
    cols = [('Dn_mem', 'D_n', "{:.3f}"), ('EB_mem', 'E_B', "{:.3f}"), ('R_PB', 'P_B n/bg', "{:.2f}"),
            ('R_Pth', 'P_th n/bg', "{:.2f}"), ('R_stat', 'P_stat(th+mag)', "{:.2f}"),
            ('R_Pdyn', 'P_dyn n/bg', "{:.2f}"), ('R_Ptot', 'P_tot(+dyn)', "{:.2f}"),
            ('R_v', 'v n/bg', "{:.2f}"), ('beta_near', 'beta_near', "{:.2f}"), ('T_near', 'T_near eV', "{:.0f}")]
    for key, lab, fmt in cols:
        line = f"  {lab:16s}"
        for name, pop in pops:
            line += f" | {name.split(' (')[0][:18]:>18s}: {mi([r[key] for r in pop], fmt):>24s}"
        o.append(line)
    o.append("")
    # compensation diagnostics on the 60
    s60 = np.array([r['R_stat'] for r in P_60 if np.isfinite(r['R_stat'])])
    t60 = np.array([r['R_Ptot'] for r in P_60 if np.isfinite(r['R_Ptot'])])
    d60 = np.array([r['R_Pdyn'] for r in P_60 if np.isfinite(r['R_Pdyn'])])
    o.append(f"  [60] static budget P_stat near/bg: median {np.median(s60):.2f}; frac>1.25 = {np.mean(s60>1.25):.2f}")
    o.append(f"  [60] ion total  P_tot  near/bg: median {np.median(t60):.2f}; frac in [0.8,1.25] = {np.mean((t60>=0.8)&(t60<=1.25)):.2f}")
    o.append(f"  [60] dynamic    P_dyn  near/bg: median {np.median(d60):.2f}  (flow stagnation if <1)")
    o.append("")

    o.append("=== CHECK 2a: event-level Spearman(D_n, E_B) across each population ===")
    o.append("(PDL => NEGATIVE: deeper depletion (low D_n) goes with stronger field (high E_B))")
    for name, pop in pops:
        dn = [r['Dn_mem'] for r in pop]; eb = [r['EB_mem'] for r in pop]
        rho, nn = spearman(dn, eb); p = perm_p(dn, eb, rho)
        o.append(f"  {name:26s} N={nn:4d}  Spearman(D_n,E_B) = {rho:+.3f}  (perm p = {p:.4f})")
    o.append("")
    o.append("=== CHECK 2b: per-event Spearman(n,|B|) over inner membership samples s in [0.05,0.40) ===")
    o.append("(PDL => most events NEGATIVE; compare clean 60 vs false-pos 32 vs all-contributing)")
    for name, pop in pops:
        rs = np.array([r['r_nB'] for r in pop if np.isfinite(r['r_nB'])])
        if len(rs):
            o.append(f"  {name:26s} N={len(rs):4d}  median r = {np.median(rs):+.3f}  "
                     f"frac<0 = {np.mean(rs<0):.2f}  frac<-0.3 = {np.mean(rs<-0.3):.2f}")
        else:
            o.append(f"  {name:26s} n/a")
    o.append("")

    # ---- read-out per user's rule ----
    rho60, _ = spearman([r['Dn_mem'] for r in P_60], [r['EB_mem'] for r in P_60])
    p60 = perm_p([r['Dn_mem'] for r in P_60], [r['EB_mem'] for r in P_60], rho60)
    rhoall, _ = spearman([r['Dn_mem'] for r in P_all], [r['EB_mem'] for r in P_all])
    pall = perm_p([r['Dn_mem'] for r in P_all], [r['EB_mem'] for r in P_all], rhoall)
    rnb60 = np.array([r['r_nB'] for r in P_60 if np.isfinite(r['r_nB'])])
    comp_ok = (np.median(t60) >= 0.8 and np.median(t60) <= 1.30)         # ion total ~ conserved
    anti_ok = ((np.isfinite(rho60) and rho60 < 0 and p60 < 0.05) or (np.isfinite(rhoall) and rhoall < 0 and pall < 0.05)) \
        and (len(rnb60) and np.mean(rnb60 < 0) >= 0.6)
    dn_all = np.median([r['Dn_mem'] for r in P_all]); eb_all = np.median([r['EB_mem'] for r in P_all])
    dn60 = np.median([r['Dn_mem'] for r in P_60]); eb60 = np.median([r['EB_mem'] for r in P_60])
    o.append("=== READ-OUT (user's decision rule) ===")
    o.append(f"  unconditioned all-contributing: D_n = {dn_all:.3f}, E_B = {eb_all:.3f} (N={len(P_all)})")
    o.append(f"  conditional spectrally-supported(60): D_n = {dn60:.3f}, E_B = {eb60:.3f}")
    o.append(f"  pressure compensation (ion total near/bg median = {np.median(t60):.2f}): {'REASONABLE' if comp_ok else 'NOT SATISFIED'}")
    o.append(f"  density-field anticorrelation (event-level rho60={rho60:+.3f} p={p60:.4f}; "
             f"all rho={rhoall:+.3f} p={pall:.4f}; per-event n-|B| frac<0(60)={np.mean(rnb60<0):.2f}): "
             f"{'HOLDS' if anti_ok else 'NOT ESTABLISHED'}")
    if comp_ok and anti_ok:
        o.append("  => VERDICT: KEEP 0.68/2.4 as a CONDITIONAL physical contrast for the selected sheath-like")
        o.append("     candidates, while the main text reports the unconditioned all-contributing D_n (~0.92).")
    else:
        o.append("  => VERDICT: DOWNGRADE the contrast result to the unconditioned screened profile (D_n~0.92) + a field")
        o.append("     enhancement that is a robust but not uniquely-coupled inward gradient; 0.68 is a")
        o.append("     selected-candidate descriptive statistic only.")
    o.append("")
    o.append("NOTE: this ion budget omits electron thermal pressure (~10-20% of ion in the sheath); P_dyn uses")
    o.append("ESA bulk speed. The check is diagnostic, not a closed total-pressure proof.")

    txt = "\n".join(o)
    print(txt, flush=True)
    with open(os.path.join(OUT, "CONTRAST_VALIDATION.txt"), "w") as f:
        f.write(txt + "\n")

    with open(os.path.join(OUT, "contrast_per_encounter.csv"), "w", newline="") as f:
        c = ['eid', 'probe', 'cls', 'spec', 'Dn_mem', 'EB_mem', 'beta_near', 'T_near', 'n_mem',
             'R_PB', 'R_Pth', 'R_stat', 'R_Pdyn', 'R_Ptot', 'R_v', 'r_nB', 'n_in', 'xcheck']
        w = csv.writer(f); w.writerow(c)
        rr = lambda v, k=4: (round(v, k) if (isinstance(v, (int, float)) and np.isfinite(v)) else '')
        for r in contrib:
            w.writerow([r['eid'], r['probe'], r['cls'], r['spec'], rr(r['Dn_mem']), rr(r['EB_mem']),
                        rr(r['beta_near']), rr(r['T_near'], 1), r['n_mem'], rr(r['R_PB'], 3), rr(r['R_Pth'], 3),
                        rr(r['R_stat'], 3), rr(r['R_Pdyn'], 3), rr(r['R_Ptot'], 3), rr(r['R_v'], 3),
                        rr(r['r_nB'], 3), r['n_in'], rr(r['xcheck'], 3)])
    print("\nsaved ->", OUT, flush=True)


if __name__ == "__main__":
    main()
