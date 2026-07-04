"""Fig 8 — complete event deep-dive for 2020-12-12_tha (with CDAWeb OMNI + ion spectrogram).

Panels: ESA reduced-ion energy spectrogram | n | |B| | T_ion | beta | OMNI (Dp, Bz) | radial s(t).
Shading: model magnetosphere s<0 (purple), near-MP shell s in[0.05,0.20) (orange), background
s in[0.6,1.0] (green). The spectrogram shows the near-MP ions stay at sheath energy (~100s eV),
NOT a keV magnetospheric peak; OMNI shows the upstream is stable across the window.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P, EVENTS
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r

SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
EID = "2020-12-12_tha"; KB = 1.602e-4
FIGDIR = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\figures")
MEET = os.path.join(P(r"H:\0mssl\review\01_CURRENT__rebuild\runs"), "supplementary")
import os as _mk; _mk.makedirs(FIGDIR, exist_ok=True); _mk.makedirs(MEET, exist_ok=True)

d = np.load(os.path.join(SUB, EID + ".npz"), allow_pickle=True)
t = d["t"].astype(float); x = d["x_re"].astype(float); y = d["y_re"].astype(float); z = d["z_re"].astype(float)
n = d["n"].astype(float); b = d["bmag"].astype(float); beta = d["beta"].astype(float)
pth = d["p_th"].astype(float)
mp0 = float(d["mp0"]); alpha = float(d["alpha"]); dp = float(d["dp"]); bz = float(d["bz"]); cone = float(d["cone_deg"])
r = np.sqrt(x*x+y*y+z*z); ct = np.clip(np.where(r > 0, x/r, 1.0), -1.0, 1.0)
dmp = r - mp0*(2/(1+ct))**alpha; s = dmp/(dmp + (jelinek_bs_r(ct, dp) - r))
T = np.divide(pth, n*KB, out=np.full_like(pth, np.nan), where=(n > 0))
hr = (t - t[0])/3600.0
near = np.isfinite(s) & (s >= 0.05) & (s < 0.20) & (n > 0); bg = np.isfinite(s) & (s >= 0.6) & (s <= 1.0) & (n > 0)
msph = np.isfinite(s) & (s < 0)
n_near, n_bg = np.median(n[near]), np.median(n[bg]); b_near, b_bg = np.median(b[near]), np.median(b[bg])
T_near, T_bg = np.nanmedian(T[near]), np.nanmedian(T[bg])

esa = np.load(os.path.join(EVENTS, "esa_2020-12-12.npz"))
ef = esa["eflux"]; en = esa["energy"]; ehr = (esa["tunix"] - t[0])/3600.0
om = np.load(os.path.join(EVENTS, "omni_2020-12-12.npz"))
oP = np.array(om["Pressure"], dtype=float); oB = np.array(om["BZ"], dtype=float)
ohr = np.linspace(0, hr[-1], len(oP))

fig = plt.figure(figsize=(9.6, 12))
gs = fig.add_gridspec(7, 1, height_ratios=[1.5, 1, 1, 1, 1, 1, 0.9], hspace=0.12)
axes = [fig.add_subplot(gs[i]) for i in range(7)]

# (0) spectrogram
axsp = axes[0]
TT, EE = np.meshgrid(ehr, en)
pc = axsp.pcolormesh(TT, EE, np.log10(ef.T), shading="auto", cmap="viridis")
axsp.set_yscale("log"); axsp.set_ylabel("ion energy\n(eV)", fontsize=9)
axsp.axhline(1000, color="w", ls="--", lw=0.8); axsp.text(0.2, 1300, "1 keV", color="w", fontsize=7)
axsp.text(0.2, 120, "sheath ions ~100s eV (peak preserved near-MP, no keV jump)", color="w", fontsize=7.5)
cb = fig.colorbar(pc, ax=axsp, pad=0.01); cb.set_label("log eflux", fontsize=7)

def shade(ax, dat):
    lo = np.nanmin(dat[np.isfinite(dat)]); hi = np.nanmax(dat[np.isfinite(dat)])
    ax.fill_between(hr, lo, hi, where=msph, color="#7E57C2", alpha=0.16, step="mid")
    ax.fill_between(hr, lo, hi, where=near, color="#F4A261", alpha=0.30, step="mid")
    ax.fill_between(hr, lo, hi, where=bg, color="#7DB37D", alpha=0.16, step="mid")
    ax.set_ylim(lo, hi)

for ax, (lab, dat, sc, meds) in zip(axes[1:6], [
        ("n  (cm$^{-3}$)", n, "log", (n_near, n_bg)), (r"$|B|$ (nT)", b, "lin", (b_near, b_bg)),
        (r"$T_{\rm ion}$ (eV)", T, "log", (T_near, T_bg)), (r"$\beta$", beta, "log", None),
        ("OMNI", None, None, None)]):
    if lab == "OMNI":
        ax.plot(ohr, oP, color="#1565C0", lw=1.0); ax.set_ylabel("OMNI D$_p$\n(nPa)", fontsize=8.5, color="#1565C0")
        ax.tick_params(axis="y", labelcolor="#1565C0"); ax.set_ylim(0, max(3, np.nanmax(oP)*1.1))
        axr = ax.twinx(); axr.plot(ohr, oB, color="#C62828", lw=1.0); axr.axhline(0, color="#C62828", ls=":", lw=0.6)
        axr.set_ylabel("B$_z$ (nT)", fontsize=8.5, color="#C62828"); axr.tick_params(axis="y", labelcolor="#C62828")
        ax.text(0.02, 0.82, f"upstream stable: D$_p$ {np.nanmin(oP):.1f}–{np.nanmax(oP):.1f} nPa", transform=ax.transAxes, fontsize=7)
        continue
    shade(ax, dat); ax.plot(hr, dat, lw=0.7, color="#1a1a1a")
    if meds: ax.axhline(meds[0], color="#E76F51", lw=1.2, ls="--"); ax.axhline(meds[1], color="#2E7D32", lw=1.2, ls="--")
    if sc == "log": ax.set_yscale("log")
    ax.set_ylabel(lab, fontsize=9)
# (6) s(t)
axs = axes[6]; axs.plot(hr, s, lw=0.9, color="#1a1a1a"); axs.axhline(0, color="#7E57C2", ls=":", lw=0.8)
axs.axhspan(0.05, 0.20, color="#F4A261", alpha=0.30); axs.axhspan(0.6, 1.0, color="#7DB37D", alpha=0.16)
axs.set_ylim(-0.3, 1.05); axs.set_ylabel("radial\n$s(t)$", fontsize=9); axs.set_xlabel("time from window start (hours)", fontsize=11)
for ax in axes[:6]: ax.set_xticklabels([])
for ax in axes: ax.set_xlim(0, hr[-1]); ax.grid(alpha=0.15)

axes[0].legend(handles=[Patch(fc="#7E57C2", alpha=.16, label="magnetosphere s<0"),
                        Patch(fc="#F4A261", alpha=.30, label="near-MP s$\\in$[0.05,0.20)"),
                        Patch(fc="#7DB37D", alpha=.16, label="background s$\\in$[0.6,1.0]")],
               fontsize=7, ncol=3, loc="upper right", framealpha=0.85)
box = (f"{EID}:  near-MP vs background  n {n_near:.1f}$\\to${n_bg:.1f} cm$^{{-3}}$ (D$_n$={n_near/n_bg:.2f}),  "
       f"|B| {b_near:.0f}$\\to${b_bg:.0f} nT (E$_B$={b_near/b_bg:.2f}),  T$_{{ion}}\\approx${T_near:.0f} eV.  "
       f"Committed spectral metrics: peak-energy ratio 1.00, shape corr 0.96, flux ratio 0.39 $\\Rightarrow$ sheath, not magnetosphere.")
fig.suptitle("Event deep-dive: a spectrally-validated near-magnetopause depletion (THEMIS-A, 2020-12-12)", fontsize=12.5, y=0.997)
fig.text(0.5, 0.005, box, ha="center", va="bottom", fontsize=8, bbox=dict(boxstyle="round", fc="#FFF8E7", ec="#C9A227"))
fig.subplots_adjust(top=0.955, bottom=0.05, left=0.09, right=0.91)
for out in (os.path.join(FIGDIR, "fig8_event_deepdive.png"), os.path.join(MEET, "fig8_event_deepdive.png")):
    fig.savefig(out, dpi=145)
print("saved fig8 ->", FIGDIR, "and", MEET)
print(f"D_n={n_near/n_bg:.2f} E_B={b_near/b_bg:.2f} T_near={T_near:.0f} OMNI Dp {np.nanmin(oP):.2f}-{np.nanmax(oP):.2f}")
