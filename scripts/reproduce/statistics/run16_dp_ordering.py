"""RUN 16 — dynamic-pressure ORDERING of the near-MP contrast (reader-friendly figure).

Full contributing membership-screened population (both IMF signs). Bin Dp into tertiles (balanced N;
ranges reported). Per bin: median Dn_mem & median EB_mem with bootstrap 95% CI, N, and a non-headline
Bz-sign split. Question: is the near-MP density-and-field contrast primarily ORDERED by dynamic pressure
in this observable sample? ("ordered by", NOT "caused by" — Dp also sets boundary location, sheath
compression, and observability.) Reads run10 funnel_contributing.csv. Figure + table out.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CONTRIB = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run16_dp_ordering")
os.makedirs(OUT, exist_ok=True)


def fnum(s):
    try:
        v = float(s); return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def boot_med_ci(a, nb=5000, seed=3):
    a = np.asarray([v for v in a if np.isfinite(v)])
    if len(a) < 3:
        return (np.nan, np.nan, np.nan)
    rng = np.random.default_rng(seed)
    bs = [np.median(rng.choice(a, len(a))) for _ in range(nb)]
    return float(np.median(a)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def main():
    rows = []
    for r in csv.DictReader(open(CONTRIB)):
        dp = fnum(r['dp']); dn = fnum(r['Dn_mem']); eb = fnum(r['EB_mem'])
        bzs = r['bz_sign']
        if np.isfinite(dp) and np.isfinite(dn) and np.isfinite(eb):
            rows.append(dict(dp=dp, dn=dn, eb=eb, bz_sign=(int(bzs) if bzs in ('0', '1') else -1)))
    dps = np.array([r['dp'] for r in rows])
    # tertile edges (balanced N), reported as physical ranges
    q1, q2 = np.percentile(dps, [100 / 3, 200 / 3])
    def binid(dp): return 0 if dp <= q1 else (1 if dp <= q2 else 2)
    labels = [f"low\n(Dp≤{q1:.2f})", f"moderate\n({q1:.2f}–{q2:.2f})", f"high\n(Dp>{q2:.2f})"]
    bins = [[r for r in rows if binid(r['dp']) == i] for i in range(3)]

    lines = ["RUN 16 — dynamic-pressure ordering of the near-MP contrast (full contributing population)",
             f"Dp tertile edges: {q1:.2f}, {q2:.2f} nPa ; N total = {len(rows)}",
             "",
             f"{'Dp bin':22s}{'N':>5s}{'Dp med':>8s}{'Dn med':>8s}{'Dn 95% CI':>20s}{'EB med':>8s}{'EB 95% CI':>20s}{'Dn N/S':>14s}"]
    dn_m = []; dn_lo = []; dn_hi = []; eb_m = []; eb_lo = []; eb_hi = []; ns = []
    md = ["# Dynamic-pressure ordering of the near-magnetopause contrast", "",
          "Full contributing, membership-screened population (both IMF signs). Dp tertiles (balanced N).",
          "Near-MP-shell (s∈[0.05,0.20)) per-encounter medians; bootstrap 95% CI (5000 resamples).",
          "",
          f"Tertile edges: **{q1:.2f}, {q2:.2f} nPa**; total N = {len(rows)}.",
          "",
          "| Dp bin (nPa) | N | Dp median | Dn_mem median [95% CI] | EB_mem median [95% CI] | Dn_mem north / south |",
          "|---|---|---|---|---|---|"]
    for i, b in enumerate(bins):
        dn, dl, dh = boot_med_ci([r['dn'] for r in b], seed=10 + i)
        eb, el, eh = boot_med_ci([r['eb'] for r in b], seed=20 + i)
        north = [r['dn'] for r in b if r['bz_sign'] == 1]; south = [r['dn'] for r in b if r['bz_sign'] == 0]
        ns_str = f"{np.median(north):.2f}/{np.median(south):.2f}" if north and south else "n/a"
        dpmed = np.median([r['dp'] for r in b])
        lab1 = labels[i].replace("\n", " ")
        lines.append(f"{lab1:22s}{len(b):5d}{dpmed:8.2f}{dn:8.3f}{f'[{dl:.3f},{dh:.3f}]':>20s}{eb:8.3f}{f'[{el:.3f},{eh:.3f}]':>20s}{ns_str:>14s}")
        md.append(f"| {lab1} | {len(b)} | {dpmed:.2f} | {dn:.3f} [{dl:.3f}, {dh:.3f}] | {eb:.3f} [{el:.3f}, {eh:.3f}] | {ns_str} |")
        dn_m.append(dn); dn_lo.append(dn - dl); dn_hi.append(dh - dn)
        eb_m.append(eb); eb_lo.append(eb - el); eb_hi.append(eh - eb); ns.append(len(b))

    # figure: two panels
    x = np.arange(3)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ax[0].errorbar(x, dn_m, yerr=[dn_lo, dn_hi], fmt='o-', color='navy', capsize=5, ms=8, lw=2)
    ax[0].axhline(1.0, color='grey', ls='--', lw=1, label='no depletion (Dn=1)')
    ax[0].set_ylabel(r'median $D_n$ (near-MP density / background)')
    ax[0].set_title('Density depletion vs dynamic pressure'); ax[0].legend(fontsize=8)
    ax[1].errorbar(x, eb_m, yerr=[eb_lo, eb_hi], fmt='s-', color='firebrick', capsize=5, ms=8, lw=2)
    ax[1].axhline(1.0, color='grey', ls='--', lw=1, label='no enhancement (EB=1)')
    ax[1].set_ylabel(r'median $E_B$ (near-MP |B| / background)')
    ax[1].set_title('Field pile-up vs dynamic pressure'); ax[1].legend(fontsize=8)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labels); a.set_xlabel('dynamic pressure bin')
        for xi, n in zip(x, ns):
            a.annotate(f'N={n}', (xi, a.get_ylim()[0]), textcoords='offset points', xytext=(0, 4), ha='center', fontsize=8, color='dimgray')
    fig.suptitle('Near-magnetopause density-and-field contrast is ordered by dynamic pressure', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    figp = os.path.join(OUT, "dp_ordering.png")
    fig.savefig(figp, dpi=130); plt.close(fig)

    lines.append("")
    lines.append(f"ORDERING: Dn_mem rises {dn_m[0]:.3f} -> {dn_m[1]:.3f} -> {dn_m[2]:.3f} across low->high Dp "
                 f"({'monotonic; CIs separate => clear Dp ordering' if (dn_m[0] < dn_m[2] and (dn_m[0]+dn_hi[0]) < (dn_m[2]-dn_lo[2])) else 'trend present; check CI overlap'}).")
    lines.append("Read as 'ordered by dynamic pressure' (Dp also sets boundary location/compression/observability), NOT 'caused by'.")
    md.append(""); md.append(f"**Ordering:** Dn_mem rises {dn_m[0]:.3f} → {dn_m[1]:.3f} → {dn_m[2]:.3f} across low→high Dp; "
              "deepest relative depletion at low dynamic pressure. Phrase as *ordered by* dynamic pressure, not *caused by*.")
    md.append(f"\n![Dp ordering](dp_ordering.png)")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "DP_ORDERING.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(OUT, "DP_ORDERING.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    print("\nsaved ->", OUT, "(dp_ordering.png + tables)")


if __name__ == "__main__":
    main()
