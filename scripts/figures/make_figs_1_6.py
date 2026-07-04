#!/usr/bin/env python
"""
Publication figures 1-6 for the PDL radial-coordinate manuscript.
Every plotted number traces to a committed run output (verified_only):
  Fig 1: schematic from Shue [1] MP + Jelinek [2] BS formulae at sample-median Dp~2.4 nPa
         (subsolar MP 9.97 R_E, BS 13.1 R_E per manuscript Sec 8.2).
  Fig 2: runs/run23_paired/paired_1d_vs_radial.csv  (N=672)
  Fig 3: runs/run2_profile_cube/profile_sheath.csv
  Fig 4: runs/run10_selection/SELECTION_FUNNEL.txt + run12 spectral (60 sheath-like)
  Fig 5: runs/run12_spectral/spectral_metrics.csv  (N=107)
  Fig 6: runs/run16_dp_ordering/DP_ORDERING.txt
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrow
from matplotlib.lines import Line2D

RUNS = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs")  # -> repo outputs/ via config.P
OUT  = os.path.join(RUNS, "figures")  # -> repo outputs/figures

# ---- shared publication style ----
plt.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200,
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
    "axes.linewidth": 0.9, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "legend.fontsize": 9, "legend.frameon": True, "legend.framealpha": 0.9,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})
C_DEN = "#1f4e9b"   # density  (navy)
C_FLD = "#b4202a"   # field    (firebrick)
C_BETA= "#2a7f3f"   # beta     (green)
C_GREY= "#666666"

def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)

# =====================================================================
# FIGURE 1 — dayside geometry schematic + radial coordinate s
# =====================================================================
def fig1():
    mp0, bs0, alpha, lam = 9.97, 13.1, 0.59, 1.17   # Sec 8.2 anchors + Shue/Jelinek forms
    th = np.linspace(0, np.deg2rad(135), 400)
    r_mp = mp0 * (2.0/(1.0+np.cos(th)))**alpha
    r_bs = bs0 * 2.0/(np.cos(th)+np.sqrt(np.cos(th)**2 + lam**2*np.sin(th)**2))
    def xy(r): return r*np.cos(th), r*np.sin(th)
    xmp,ymp = xy(r_mp); xbs,ybs = xy(r_bs)

    def s_contour(s):
        r = r_mp + s*(r_bs - r_mp); return r*np.cos(th), r*np.sin(th)

    fig, ax = plt.subplots(figsize=(7.2,6.2), constrained_layout=True)
    XL = -2.0   # left plot edge; close all region polygons here so no wedge/gap remains
    # solar wind = whole panel background (everything OUTSIDE the bow shock keeps this colour)
    ax.add_patch(plt.Rectangle((XL,-11.4), 23.5-XL, 22.8, color="#f2f6fc", zorder=-1, lw=0))
    # magnetosheath = inside the bow shock (close the bowl along the left edge -> fills solidly)
    xbowl=np.r_[xbs,[XL,XL],xbs[::-1]]; ybowl=np.r_[ybs,[ybs[-1],-ybs[-1]],-ybs[::-1]]
    ax.fill(xbowl, ybowl, color="#d6e4f4", zorder=1, lw=0)
    # magnetosphere = inside the magnetopause (solid closed polygon -> no white wedge at Earth)
    xmpol=np.r_[xmp,[XL,XL],xmp[::-1]]; ympol=np.r_[ymp,[ymp[-1],-ymp[-1]],-ymp[::-1]]
    ax.fill(xmpol, ympol, color="#e9e3f2", zorder=2, lw=0)
    # near-shell s in [0.05,0.20)
    xa,ya = s_contour(0.05); xb,yb = s_contour(0.20)
    for sgn in (1,-1):
        ax.fill(np.r_[xb, xa[::-1]], sgn*np.r_[yb, ya[::-1]],
                color="#f6c9a8", alpha=0.95, zorder=3, lw=0)
    # background s in [0.6,1.0]
    xc,yc = s_contour(0.60)
    for sgn in (1,-1):
        ax.fill(np.r_[xbs, xc[::-1]], sgn*np.r_[ybs, yc[::-1]],
                color="#cfe6cf", alpha=0.80, zorder=3, lw=0)
    # boundary curves
    for sgn in (1,-1):
        ax.plot(xmp, sgn*ymp, color="#3b2f6b", lw=2.0, zorder=4)
        ax.plot(xbs, sgn*ybs, color="#7a2d2d", lw=2.0, zorder=4)
    # Earth
    ax.add_patch(plt.Circle((0,0), 0.5, color="#2a4d69", zorder=6))
    ax.plot(0,0, marker=r"$\oplus$", ms=11, color="white", zorder=7)
    # Sun arrow
    ax.annotate("", xy=(20.3,0), xytext=(17.5,0),
                arrowprops=dict(arrowstyle="-|>", color="#d98a00", lw=2.2))
    ax.text(20.5,0.0,"to Sun", va="center", ha="left", color="#d98a00", fontsize=10)

    ax.text(11.0, 0.0, "MP\n$s\\,{=}\\,0$", ha="center", va="center", fontsize=9.5,
            color="#3b2f6b", fontweight="bold")
    ax.text(13.6, 4.6, "BS  $s\\,{=}\\,1$", ha="center", va="center", fontsize=9.5,
            color="#7a2d2d", fontweight="bold", rotation=-38)
    ax.text(10.6, 2.7, "near-MP shell\n$s\\in[0.05,0.20)$", ha="center", va="center",
            fontsize=9, color="#8a4a1a")
    ax.text(8.2, 7.6, "background\n$s\\in[0.6,1.0]$", ha="center", va="center",
            fontsize=9, color="#2a5c2a")
    ax.text(5.0, 0.0, "magnetosphere", ha="center", va="center", fontsize=9, color="#5b4b86")
    ax.text(9.4, -5.6, "magnetosheath", ha="center", va="center", fontsize=9.5, color="#2f4f7a")
    ax.text(17.0, 6.4, "solar wind", ha="center", va="center", fontsize=9.5, color="#5f7fb0")

    ax.set_xlabel(r"$X_{\rm GSE}$  ($R_E$)"); ax.set_ylabel(r"$\rho=\sqrt{Y^2+Z^2}$  ($R_E$)")
    ax.set_title("Dayside boundary geometry and the radial normalised coordinate $s$")
    ax.set_xlim(-2, 23.5); ax.set_ylim(-11, 11); ax.set_aspect("equal")
    ax.grid(alpha=0.15)
    save(fig, "fig1_geometry.png")

# =====================================================================
# FIGURE 2 — (a) the two coordinate constructions + (b,c) paired comparison
# =====================================================================
def fig2():
    rows = list(csv.DictReader(open(os.path.join(RUNS,"run23_paired","paired_1d_vs_radial.csv"))))
    d1=np.array([float(r["Dn_1d"]) for r in rows]); dr=np.array([float(r["Dn_radial"]) for r in rows])
    e1=np.array([float(r["EB_1d"]) for r in rows]); er=np.array([float(r["EB_radial"]) for r in rows])
    N=len(rows)
    def med(a): return float(np.median(a))
    rng=np.random.default_rng(7)

    # ---- panel (a): same Shue/Jelinek surfaces as Fig 1 (sample-median Dp anchors) ----
    mp0, bs0, alpha, lam = 9.97, 13.1, 0.59, 1.17
    th = np.linspace(0, np.deg2rad(62), 320)
    r_mp = mp0 * (2.0/(1.0+np.cos(th)))**alpha
    r_bs = bs0 * 2.0/(np.cos(th)+np.sqrt(np.cos(th)**2 + lam**2*np.sin(th)**2))
    xmp,ymp = r_mp*np.cos(th), r_mp*np.sin(th)
    xbs,ybs = r_bs*np.cos(th), r_bs*np.sin(th)
    def s_arc(s):
        r = r_mp + s*(r_bs - r_mp); return r*np.cos(th), r*np.sin(th)
    # illustrative sample: solar-zenith angle 30 deg (the orbital-screen edge), chosen inside
    # the fixed-axis broad near-boundary bin; radially the same point is outer sheath
    thS = np.deg2rad(30.0)
    xS  = 10.9; rS = xS/np.cos(thS); yS = rS*np.sin(thS)
    slab_lo, slab_hi = mp0+0.2*(bs0-mp0), mp0+0.4*(bs0-mp0)   # fixed-axis bin 0.2-0.4 along X

    fig = plt.figure(figsize=(9.8, 9.4), constrained_layout=True)
    gs  = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.08])

    def schematic(ax, mode):
        ax.add_patch(plt.Rectangle((6,0), 10.2, 9.4, color="#f2f6fc", zorder=0, lw=0))
        ax.fill(np.r_[xmp, xbs[::-1]], np.r_[ymp, ybs[::-1]], color="#d6e4f4", zorder=1, lw=0)
        ax.fill(np.r_[xmp,[6,6]], np.r_[ymp,[ymp[-1],0]], color="#e9e3f2", zorder=1, lw=0)
        ax.plot(xmp, ymp, color="#3b2f6b", lw=2.0, zorder=4)
        ax.plot(xbs, ybs, color="#7a2d2d", lw=2.0, zorder=4)
        ax.text(9.15, 0.55, "MP", color="#3b2f6b", fontsize=8.5, fontweight="bold")
        ax.text(12.55, 0.55, "BS", color="#7a2d2d", fontsize=8.5, fontweight="bold")
        ax.text(7.0, 1.5, "magneto-\nsphere", fontsize=8, color="#5b4b86", ha="center")
        ax.text(8.8, 7.15, "magneto-\nsheath", fontsize=7.5, color="#2f4f7a", ha="center")
        if mode=="fixed":
            ax.axvspan(slab_lo, slab_hi, color="#f6c9a8", alpha=0.85, zorder=2, lw=0)
            for xe in (slab_lo, slab_hi):
                ax.plot([xe,xe],[0,9.4], ls="--", color="#c07830", lw=1.0, zorder=3)
            ax.annotate("", xy=(xS, 0.25), xytext=(xS, yS),
                        arrowprops=dict(arrowstyle="->", color=C_GREY, lw=1.0, ls=":"), zorder=5)
            ax.text(0.5*(slab_lo+slab_hi), 8.75, "near-boundary bin\n(fixed-axis $s\\in[0.2,0.4)$)",
                    fontsize=7.8, color="#8a4a1a", ha="center", va="top", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#d8c0a0", lw=0.5, alpha=0.95))
            ax.annotate("sample falls inside\nthe near-boundary bin", (xS,yS), xytext=(13.15,7.6),
                        fontsize=8, ha="left", va="center", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.25", fc="#fbf7e8", ec="#caa", lw=0.6),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
            ax.set_title("distances measured along the fixed Sun–Earth axis", fontsize=9.5)
            ax.text(0.02, 0.965, "(a)", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
        else:
            xa,ya = s_arc(0.2); xb,yb = s_arc(0.4)
            ax.fill(np.r_[xa, xb[::-1]], np.r_[ya, yb[::-1]], color="#f6c9a8", alpha=0.85, zorder=2, lw=0)
            for (xe,ye) in ((xa,ya),(xb,yb)):
                ax.plot(xe, ye, ls="--", color="#c07830", lw=1.0, zorder=3)
            ax.plot([6, 16.2], [6*np.tan(thS), 16.2*np.tan(thS)], ls=":", color=C_GREY, lw=1.2, zorder=3)
            ax.text(6.5, 3.9, "Earth–spacecraft direction", fontsize=7.6,
                    color=C_GREY, ha="left", va="bottom", rotation=30)
            ax.text(10.05, 8.75, "near-boundary bin\n(radial $s\\in[0.2,0.4)$)",
                    fontsize=7.8, color="#8a4a1a", ha="center", va="top", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#d8c0a0", lw=0.5, alpha=0.95))
            ax.annotate("the same sample sits\nin the outer sheath", (xS,yS), xytext=(13.15,3.1),
                        fontsize=8, ha="left", va="center", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.25", fc="#fbf7e8", ec="#caa", lw=0.6),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
            ax.set_title("distances measured along the Earth–spacecraft direction", fontsize=9.5)
        ax.plot([xS],[yS], marker="o", ms=9, mfc="white", mec="black", mew=1.8, zorder=8)
        ax.set_xlim(6,16.2); ax.set_ylim(0,9.4); ax.set_aspect("equal")
        ax.set_xlabel(r"$X_{\rm GSE}$  ($R_E$)", fontsize=8.5)
        ax.set_ylabel(r"$\rho$  ($R_E$)", fontsize=8.5)
        ax.tick_params(labelsize=8); ax.grid(alpha=0.15)

    axA = fig.add_subplot(gs[0,0]); schematic(axA, "fixed")
    axB = fig.add_subplot(gs[0,1]); schematic(axB, "radial")

    axs = [fig.add_subplot(gs[1,0]), fig.add_subplot(gs[1,1])]
    def panel(ax, a, b, lab, color, refline, ymax):
        jx=rng.normal(0,0.05,len(a))
        ax.scatter(1+jx, a, s=6, color=color, alpha=0.18, zorder=2, edgecolors="none")
        ax.scatter(2+jx, b, s=6, color=color, alpha=0.18, zorder=2, edgecolors="none")
        bp=ax.boxplot([a,b], positions=[1,2], widths=0.28, showfliers=False,
                      patch_artist=True, zorder=3,
                      medianprops=dict(color="black",lw=1.8),
                      boxprops=dict(facecolor="white",alpha=0.85,lw=1.2),
                      whiskerprops=dict(lw=1.0), capprops=dict(lw=1.0))
        if refline is not None:
            ax.axhline(refline, ls="--", color=C_GREY, lw=1.0, zorder=0)
        ma, mb = med(a), med(b)
        ax.annotate(f"{ma:.3f}", (1,ma), xytext=(0.62,ma), fontsize=10, fontweight="bold",
                    color=color, va="center", ha="right")
        ax.annotate(f"{mb:.3f}", (2,mb), xytext=(2.38,mb), fontsize=10, fontweight="bold",
                    color=color, va="center", ha="left")
        ax.set_xlim(0.4,2.6); ax.set_ylim(0,ymax)
        ax.set_xticks([1,2]); ax.set_xticklabels(["fixed-axis\ncoordinate","radial\ncoordinate"])
        ax.set_ylabel(lab)
    panel(axs[0], d1, dr, r"density contrast $D_n$ (near-MP / background)", C_DEN, 1.0, 2.6)
    panel(axs[1], e1, er, r"field contrast $E_B$ (near-MP / background)", C_FLD, 1.0, 5.0)
    axs[0].set_title("(b) density: apparent depletion is a coordinate artefact", fontsize=10)
    axs[1].set_title("(c) field: enhancement moves oppositely", fontsize=10)
    axs[0].text(0.5,0.985,
        f"paired, same {N} encounters\nmedian paired shift $+{np.median(dr-d1):.3f}$, {100*np.mean(dr>d1):.0f}% up\nWilcoxon $p=4\\times10^{{-44}}$",
        transform=axs[0].transAxes, fontsize=8.5, va="top", ha="center",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fbf7e8", ec="#caa", lw=0.6))
    fig.suptitle("Same boundary models, different measurement direction:\n"
                 "the fixed-axis coordinate manufactures the broad-bin depletion",
                 fontsize=12, fontweight="bold")
    save(fig, "fig2_paired_artefact.png")

# =====================================================================
# FIGURE 3 — radial fine-shell profile
# =====================================================================
def fig3():
    rows=list(csv.DictReader(open(os.path.join(RUNS,"run2_profile_cube","profile_sheath.csv"))))
    slo=np.array([float(r["s_lo"]) for r in rows]); shi=np.array([float(r["s_hi"]) for r in rows])
    smid=0.5*(slo+shi)
    Dn=np.array([float(r["Dn"]) for r in rows]); EB=np.array([float(r["EB"]) for r in rows])
    beta=np.array([float(r["beta_med"]) for r in rows]); magf=np.array([float(r["mag_frac"]) for r in rows])

    fig, ax = plt.subplots(figsize=(8.4,5.4), constrained_layout=True)
    ax.axhline(1.0, ls="--", color=C_GREY, lw=1.0, zorder=0)
    # analysis shell shading
    ax.axvspan(0.05,0.20, color="#f6c9a8", alpha=0.35, zorder=0, label="analysis shell $s\\in[0.05,0.20)$")
    ax.axvspan(0.00,0.10, color="#d9534f", alpha=0.10, zorder=0)
    l1,=ax.plot(smid,Dn,"o-",color=C_DEN,lw=2.0,ms=6,label=r"$D_n$ (density / background)")
    l2,=ax.plot(smid,EB,"s-",color=C_FLD,lw=2.0,ms=6,label=r"$E_B$ (field / background)")
    ax.set_xlabel(r"normalised magnetosheath coordinate $s$  (0 = magnetopause, 1 = bow shock)")
    ax.set_ylabel("contrast ratio (near / background)")
    ax.set_xlim(0,1.0); ax.set_ylim(0,2.8)
    ax2=ax.twinx(); ax2.set_yscale("log")
    l3,=ax2.plot(smid,beta,"^:",color=C_BETA,lw=1.6,ms=6,label=r"$\beta$ (plasma beta)")
    ax2.set_ylabel(r"plasma $\beta$ (log)", color=C_BETA); ax2.tick_params(axis="y",colors=C_BETA,which="both")
    ax2.set_ylim(0.1, 7.0)
    ax2.set_yticks([0.1,0.2,0.5,1,2,5])
    ax2.set_yticklabels(["0.1","0.2","0.5","1","2","5"])
    ax2.minorticks_off(); ax2.grid(False)
    # contamination annotation — placed in clear space (upper area), white bbox, away from the E_B line
    # drawn on the twin (beta) axis with high zorder so it sits ABOVE the green beta line
    ax2.annotate(f"inner sub-shells contamination-prone\n($\\beta\\!<\\!0.2$: {magf[0]*100:.0f}% of $s\\in[0,0.05)$, {magf[1]*100:.0f}% of $[0.05,0.10)$)",
                 xy=(0.065,0.60), xytext=(0.285,0.965), xycoords="axes fraction", textcoords="axes fraction",
                 fontsize=8.2, color="#9c2b27", ha="left", va="top", zorder=30,
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#d8b0ad", lw=0.6, alpha=1.0),
                 arrowprops=dict(arrowstyle="->",color="#9c2b27",lw=1.1))
    ax.annotate("clean contrast\n$s\\in[0.10,0.20)$", (0.155,0.85), xytext=(0.45,0.40),
                fontsize=8.6, color="#8a4a1a", ha="left", va="center", zorder=8,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#d8c0a0", lw=0.6, alpha=1.0),
                arrowprops=dict(arrowstyle="->",color="#8a4a1a",lw=1.0))
    lns=[l1,l2,l3]; ax.legend(lns,[l.get_label() for l in lns], loc="lower right", framealpha=1.0)
    ax.set_title("Radial fine-shell profile: density falls and field rises toward the magnetopause")
    save(fig,"fig3_fineshell_profile.png")

# =====================================================================
# FIGURE 4 — selection-function funnel
# =====================================================================
def fig4():
    # committed from run10/SELECTION_FUNNEL.txt (+ run12: 60 spectrally supported)
    stages=[
        ("full substrate archive", 6248, ""),
        ("geometric magnetosheath", 3869, "window in modelled sheath"),
        ("sampled outer-sheath background", 1187, "background $s\\in[0.6,1.0]$"),
        ("contributing (membership near-shell)", 661, "screened $s\\in[0.05,0.20)$"),
        ("northward IMF ($B_z>0$)", 332, "low-shear subset"),
        ("moment-classified candidates", 107, "candidate near-magnetopause depletion"),
        ("spectrally supported (sheath-like)", 60, "ion-spectral validation"),
    ]
    labels=[s[0] for s in stages]; N=[s[1] for s in stages]; why=[s[2] for s in stages]
    fig, ax = plt.subplots(figsize=(10.6,6.6), constrained_layout=True)
    y=np.arange(len(stages))[::-1]*1.0
    wmax=N[0]   # longest bar (full archive) normalised to 1.0 so the axis fraction is honest
    DROPX, WHYX = 1.32, 1.48   # fixed right-hand columns, clear of the longest (=1.0) bar + its N label
    cmap=plt.cm.viridis(np.linspace(0.15,0.85,len(stages)))
    for i,(yi,n,c) in enumerate(zip(y,N,cmap)):
        w=n/wmax
        ax.barh(yi, w, height=0.46, color=c, edgecolor="white", lw=0.8, align="center", zorder=2)
        # stage label sits ABOVE the bar (own line) -> never overlaps the bar or the numbers
        ax.text(0.0, yi+0.34, labels[i], va="bottom", ha="left", fontsize=9.6, fontweight="bold")
        # N at the end of the bar (the honest 0-1 axis already conveys the fraction)
        ax.text(w+0.008, yi, f"$N={n}$", va="center", ha="left", fontsize=9.8,
                fontweight="bold", color=c)
        # drop + reason in fixed right-hand columns (clear of the longest =1.0 bar's N label)
        if i>0:
            ax.text(DROPX, yi, f"$-{N[i-1]-n}$", va="center", ha="left", fontsize=9, color=C_GREY)
            ax.text(WHYX, yi, why[i], va="center", ha="left", fontsize=8.6, color=C_GREY, style="italic")
    ax.text(DROPX, y[0], "drop", va="center", ha="left", fontsize=8.5, color=C_GREY, fontweight="bold")
    ax.text(WHYX, y[0], "filter applied", va="center", ha="left", fontsize=8.5, color=C_GREY, fontweight="bold")
    ax.set_xlim(0,2.18); ax.set_ylim(-0.7, y[0]+0.9); ax.set_yticks([])
    ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_xlabel("fraction of full archive (bar length)")
    ax.set_title("Observability funnel — each step is a filter, not a non-detection", pad=14)
    ax.text(0.5,-0.11,"The candidate counts are selection-limited; 107 and 60 are NOT occurrence rates.",
            transform=ax.transAxes, ha="center", fontsize=8.8, style="italic", color=C_GREY)
    ax.grid(False)
    for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
    save(fig,"fig4_selection_funnel.png")

# =====================================================================
# FIGURE 5 — spectral-tier validation
# =====================================================================
def fig5():
    rows=list(csv.DictReader(open(os.path.join(RUNS,"run12_spectral","spectral_metrics.csv"))))
    st=np.array([r["status"] for r in rows])
    pr=np.array([float(r["peak_ratio"]) for r in rows])
    sc=np.array([float(r["shape_corr"]) for r in rows])
    fr=np.array([float(r["flux_ratio"]) for r in rows])
    pr=np.clip(pr,0.08,120)
    style={"SHEATH_CONSISTENT":("sheath-like (60)",C_BETA,"o"),
           "AMBIGUOUS_SPEC":("borderline (15)","#999999","D"),
           "SHAPE_FLAG":("shape-flagged (4)","#e08a1e","^"),
           "HOT_BOUNDARY_FLAG":("hot-boundary false positive (28)",C_FLD,"X")}

    fig, axs = plt.subplots(1,2, figsize=(10.4,5.2), constrained_layout=True,
                            gridspec_kw={"width_ratios":[1.5,1]})
    ax=axs[0]
    # sheath-like acceptance box: peak in [0.4,2.5], shape>=0.8
    ax.axhspan(0.4,2.5, xmin=(0.8+1)/2, color="#2a7f3f", alpha=0.10, zorder=0)
    ax.axhline(2.5, ls="--", color=C_GREY, lw=1.0); ax.axvline(0.8, ls=":", color=C_GREY, lw=1.1)
    ax.axvline(0.6, ls=":", color="#b9a93a", lw=1.1)
    for k,(lab,col,mk) in style.items():
        m=st==k
        ec = "#5a1417" if mk=="X" else "white"   # dark edge so faint hot-boundary X's stay visible
        ax.scatter(sc[m], pr[m], s=46, c=col, marker=mk, alpha=0.9, edgecolors=ec,
                   linewidths=0.5, label=lab, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(0.08,120); ax.set_xlim(-1.05,1.05)
    ax.set_xlabel("spectral-shape correlation (near vs background)")
    ax.set_ylabel("peak-energy ratio (near / background)")
    # threshold-line labels (clear of the data)
    ax.text(-1.02, 2.85, "peak-energy ratio = 2.5", fontsize=7.2, color=C_GREY, ha="left", va="bottom")
    ax.text(0.805, 0.092, "0.8", fontsize=7.2, color=C_GREY, ha="left", va="bottom")
    ax.text(0.615, 0.092, "0.6", fontsize=7.2, color="#8a7d1f", ha="left", va="bottom")
    # region labels in clear space, white bbox, no overlap with the green cluster
    ax.text(-0.95, 38, "keV-shifted\n(hot boundary)", fontsize=8, color=C_FLD, ha="left", va="center",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#d9b3b3", lw=0.5, alpha=0.9))
    ax.annotate("sheath-like\nacceptance region", (0.9,0.9), xytext=(0.02,0.20),
                fontsize=8, color=C_BETA, ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#aacaaa", lw=0.5, alpha=0.9),
                arrowprops=dict(arrowstyle="->", color=C_BETA, lw=1.0))
    ax.legend(loc="lower left", fontsize=8.0, ncol=1, framealpha=0.95)
    ax.set_title("Spectral metrics separate sheath-like ions from hot-boundary contamination")
    # right: flux ratio by tier
    ax2=axs[1]
    order=["SHEATH_CONSISTENT","AMBIGUOUS_SPEC","SHAPE_FLAG","HOT_BOUNDARY_FLAG"]
    data=[fr[st==k] for k in order]; cols=[style[k][1] for k in order]
    bp=ax2.boxplot(data, vert=True, patch_artist=True, showfliers=False, widths=0.6,
                   medianprops=dict(color="black",lw=1.6))
    for patch,c in zip(bp["boxes"],cols): patch.set_facecolor(c); patch.set_alpha(0.6)
    ax2.set_xticks(range(1,5)); ax2.set_xticklabels(["sheath\n(60)","border\n(15)","shape\n(4)","hot-bdy\n(28)"],fontsize=8.5)
    ax2.set_ylabel("integrated flux ratio (near / background)")
    ax2.set_title("Flux suppression by tier")
    fig.suptitle("Per-event ion-spectral validation of the 107 moment-classified candidates",
                 fontsize=12.5, fontweight="bold")
    save(fig,"fig5_spectral_tiers.png")

# =====================================================================
# FIGURE 6 — dynamic-pressure ordering (re-render, house style)
# =====================================================================
def fig6():
    # committed from run16/DP_ORDERING.txt
    binlab=["low\n($D_p\\leq1.91$)","moderate\n(1.91–3.14)","high\n($D_p>3.14$)"]
    N=[223,218,220]
    Dn=[0.711,0.928,1.158]; Dn_lo=[0.621,0.847,1.046]; Dn_hi=[0.796,1.010,1.234]
    EB=[2.257,2.077,1.633]; EB_lo=[2.151,1.980,1.516]; EB_hi=[2.399,2.228,1.741]
    x=np.arange(3)
    fig, axs=plt.subplots(1,2, figsize=(9.6,5.0), constrained_layout=True)
    def panel(ax,y,lo,hi,color,mk,lab,ref):
        y=np.array(y); lo=np.array(lo); hi=np.array(hi)
        err=[y-lo, hi-y]
        ax.errorbar(x,y,yerr=err,fmt=mk+"-",color=color,lw=2.0,ms=9,capsize=4,capthick=1.2,
                    mec="white",mew=0.6)
        ax.axhline(ref, ls="--", color=C_GREY, lw=1.0)
        off=0.03*(hi.max()-lo.min())
        for xi,ni,hh in zip(x,N,hi):
            ax.text(xi, hh+off, f"$N={ni}$", fontsize=8.5, color=C_GREY, ha="center", va="bottom")
        ax.set_xticks(x); ax.set_xticklabels(binlab); ax.set_xlim(-0.4,2.4)
        ax.set_ylabel(lab)
    panel(axs[0],Dn,Dn_lo,Dn_hi,C_DEN,"o",r"median $D_n$ (near-MP density / background)",1.0)
    axs[0].set_ylim(0.55,1.30)
    panel(axs[1],EB,EB_lo,EB_hi,C_FLD,"s",r"median $E_B$ (near-MP $|B|$ / background)",1.0)
    axs[1].set_ylim(0.95,2.55)
    axs[0].set_title("Density contrast vs dynamic pressure")
    axs[1].set_title("Field contrast vs dynamic pressure")
    for ax in axs: ax.set_xlabel("solar-wind dynamic-pressure tertile")
    axs[0].text(0.03,0.04,"error bars: bootstrap 95% CI;  dashed line: ratio = 1 (no contrast)",
                transform=axs[0].transAxes, fontsize=7.6, color=C_GREY, va="bottom", ha="left")
    fig.suptitle("Contrast is ordered by dynamic pressure in the model coordinate (not a demonstrated physical driver)",
                 fontsize=11.5, fontweight="bold")
    save(fig,"fig6_dp_ordering.png")

if __name__=="__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print("DONE figs 1-6")
