"""Candidate mini-atlas: the 6 cleanest spectrally supported events on one sheet.

Answers the reviewer's "have you looked at anything with regards to those events?" with a
contact sheet: per event, stacked n / |B| / T / s(t) panels with the magnetosphere (purple),
near-MP shell (orange) and background (green) bands, footed with the committed run13 metrics.
Framing: spectrally supported CANDIDATES, not confirmed PDLs.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r

SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
VAL = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run13_validation\validation_status.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions")
MEET = os.path.join(P(r"H:\0mssl\review\01_CURRENT__rebuild\runs"), "supplementary")
KB = 1.602e-4

with open(VAL, newline="") as f:
    rows = [r for r in csv.DictReader(f) if r["spec_status"] == "SHEATH_CONSISTENT"]
rows.sort(key=lambda r: int(r["rank"]))
top = rows[:6]

plt.rcParams.update({"font.family": "serif", "font.size": 7})
fig = plt.figure(figsize=(12.5, 8.2))
outer = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.24, left=0.05, right=0.985, top=0.90, bottom=0.05)

for k, v in enumerate(top):
    eid = v["eid"]
    d = np.load(os.path.join(SUB, eid + ".npz"), allow_pickle=True)
    t = d["t"].astype(float)
    x, y, z = (d[c].astype(float) for c in ("x_re", "y_re", "z_re"))
    n = d["n"].astype(float); b = d["bmag"].astype(float); pth = d["p_th"].astype(float)
    mp0, alpha, dp = float(d["mp0"]), float(d["alpha"]), float(d["dp"])
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.where(r > 0, x / r, 1.0), -1.0, 1.0)
    dmp = r - mp0 * (2 / (1 + ct)) ** alpha
    s = dmp / (dmp + (jelinek_bs_r(ct, dp) - r))
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    hr = (t - t[0]) / 3600.0
    near = np.isfinite(s) & (s >= 0.05) & (s < 0.20)
    bg = np.isfinite(s) & (s >= 0.6) & (s <= 1.0)
    msp = np.isfinite(s) & (s < 0)

    inner = outer[k].subgridspec(4, 1, hspace=0.10)
    axs = [fig.add_subplot(inner[i]) for i in range(4)]

    def shade(ax):
        lo, hi = ax.get_ylim()
        ax.fill_between(hr, lo, hi, where=msp, color="#7E57C2", alpha=0.15, step="mid")
        ax.fill_between(hr, lo, hi, where=near, color="#F4A261", alpha=0.30, step="mid")
        ax.fill_between(hr, lo, hi, where=bg, color="#7DB37D", alpha=0.15, step="mid")
        ax.set_ylim(lo, hi)

    axs[0].plot(hr, n, lw=0.5, color="#1f3a6e"); axs[0].set_yscale("log")
    axs[0].set_ylabel("n\n(cm$^{-3}$)", fontsize=6)
    axs[1].plot(hr, b, lw=0.5, color="#8c2d2d"); axs[1].set_ylabel("|B|\n(nT)", fontsize=6)
    axs[2].plot(hr, T, lw=0.5, color="#3a6e3a"); axs[2].set_yscale("log")
    axs[2].axhline(1000, color="0.4", ls=":", lw=0.7); axs[2].set_ylabel("T\n(eV)", fontsize=6)
    axs[3].plot(hr, s, lw=0.6, color="0.25"); axs[3].axhline(0, color="#7E57C2", lw=0.7)
    axs[3].set_ylabel("s", fontsize=6); axs[3].set_ylim(-0.25, 1.05)
    axs[3].set_xlabel("hours from window start", fontsize=6)
    for ax in axs:
        shade(ax)
        ax.tick_params(labelsize=5.5)
        if ax is not axs[3]:
            ax.set_xticklabels([])
    axs[0].set_title(
        f"{eid}  (rank {v['rank']})\n"
        f"$D_n$={float(v['Dn_mem']):.2f}  $E_B$={float(v['EB_mem']):.2f}  |  "
        f"pk {float(v['peak_ratio']):.2f}  shp {float(v['shape_corr']):.2f}  flx {float(v['flux_ratio']):.2f}",
        fontsize=6.5, pad=2)

fig.suptitle(
    "The six cleanest spectrally supported near-magnetopause depletion CANDIDATES (of 60; ranked by run13 spectral cleanliness)\n"
    "orange = near-MP shell s∈[0.05,0.20) · green = background s∈[0.6,1.0] · purple = model magnetosphere s<0 · dotted = 1 keV",
    fontsize=9.5)
fig.savefig(os.path.join(OUT, "candidate_atlas_top6.png"), dpi=165)
import shutil
shutil.copy(os.path.join(OUT, "candidate_atlas_top6.png"), os.path.join(MEET, "candidate_atlas_top6.png"))
print("saved -> run28_candidate_conditions/candidate_atlas_top6.png (+ supplementary copy)")
print("events:", ", ".join(v["eid"] for v in top))
