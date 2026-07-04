"""run26 — log-ratio presentation of the near-MP density contrast (log-scale presentation).

Motivation: plain ratios are asymmetric (a x2 enhancement looks farther from 1 than a x2
depletion); use log(n_near) - log(n_bg), symmetric about 0; and check the small near-zero group.
Reuses the committed membership pipeline (run7 functions) to dump per-encounter Dn_mem for the
headline shell s in [0.05,0.20), then reports the log10 distribution + renders the histogram.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from run_hardening_checks import load_enc, member_mask, shell_dn
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run26_logratio")
os.makedirs(OUT, exist_ok=True)


def per_enc_dn(d):
    e = load_enc(d)
    if e is None:
        return None
    dn, _ = shell_dn(e, member_mask(e, 3.0, 0.2, 0.3), 0.05, 0.20)
    return dn


def main():
    dns = np.array([d for d in pmap(per_enc_dn) if d is not None and d > 0])
    lg = np.log10(dns)
    N = len(dns)
    med, q1, q3 = np.median(dns), np.percentile(dns, 25), np.percentile(dns, 75)
    lmed, lq1, lq3 = np.median(lg), np.percentile(lg, 25), np.percentile(lg, 75)
    out = [
        f"N encounters (membership-screened, shell [0.05,0.20)) = {N}",
        "",
        "=== RAW RATIO D_n = n_near / n_bg ===",
        f"  median = {med:.3f}   IQR = [{q1:.3f}, {q3:.3f}]",
        f"  fraction D_n < 1 (depleted)       = {np.mean(dns < 1):.1%}",
        f"  fraction D_n < 0.5 (deep/near-zero) = {np.mean(dns < 0.5):.1%}   (N={int(np.sum(dns < 0.5))})",
        f"  fraction D_n < 0.3                = {np.mean(dns < 0.3):.1%}   (N={int(np.sum(dns < 0.3))})",
        "",
        "=== LOG-RATIO log10(D_n) = log10(n_near) - log10(n_bg)  [symmetric about 0] ===",
        f"  median = {lmed:+.3f}   IQR = [{lq1:+.3f}, {lq3:+.3f}]   (10^median = {10**lmed:.3f})",
        f"  fraction < 0 (depletion) = {np.mean(lg < 0):.1%}    fraction > 0 (enhancement) = {np.mean(lg > 0):.1%}",
        "",
        "Note: the raw-ratio median (0.92) and the log median agree (10^median ~ 0.92); the log view",
        "only makes depletion vs enhancement visually symmetric. The deep D_n<0.5 group is the",
        "contamination-prone tail; per-event ion-spectral validation (run12/run19) shows ~1/3 of the",
        "deepest moment-selected candidates are boundary-layer/magnetospheric false positives.",
    ]
    txt = "\n".join(out)
    print(txt, flush=True)
    with open(os.path.join(OUT, "LOGRATIO.txt"), "w") as f:
        f.write(txt + "\n")

    # --- histogram figure ---
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 3.6))
    axL.hist(dns, bins=np.linspace(0, 3, 46), color="#4C72B0", edgecolor="white", linewidth=0.4)
    axL.axvline(1.0, color="0.3", ls="--", lw=1); axL.axvline(med, color="#C44E52", lw=1.5)
    axL.set_xlabel(r"$D_n = n_{\rm near}/n_{\rm bg}$ (raw ratio)"); axL.set_ylabel("encounters")
    axL.set_title("Raw ratio (asymmetric)", fontsize=10)
    axL.text(med, axL.get_ylim()[1]*0.9, f" median {med:.2f}", color="#C44E52", fontsize=8)
    axR.hist(lg, bins=np.linspace(-1.2, 1.2, 46), color="#55A868", edgecolor="white", linewidth=0.4)
    axR.axvline(0.0, color="0.3", ls="--", lw=1); axR.axvline(lmed, color="#C44E52", lw=1.5)
    axR.set_xlabel(r"$\log_{10}(D_n)$ (symmetric)"); axR.set_ylabel("encounters")
    axR.set_title("Log-ratio (symmetric about 0)", fontsize=10)
    axR.text(lmed, axR.get_ylim()[1]*0.9, f" median {lmed:+.2f}", color="#C44E52", fontsize=8)
    fig.suptitle(f"Near-magnetopause density contrast, membership-screened (N={N})", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "logratio_histogram.png"), dpi=140)
    print("\nsaved -> run26_logratio/LOGRATIO.txt + logratio_histogram.png", flush=True)


if __name__ == "__main__":
    main()
