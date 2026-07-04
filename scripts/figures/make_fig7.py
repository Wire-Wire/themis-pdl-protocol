#!/usr/bin/env python
"""
Figure 7 — three event archetypes, built ONLY from committed data:
  per-encounter substrate time series (substrate/<eid>.npz: t,n,bmag,beta,p_th,vmag,position,mp0,bs0,alpha)
  + committed ion-spectral classification metrics (run12_spectral/spectral_metrics.csv).
No ion spectrogram is fabricated (raw ESA spectra are not committed); each column's spectral
verdict is shown as the committed metrics (peak-energy ratio / shape correlation / flux ratio).

Archetypes:
  (a) spectrally supported sheath-like depletion  : 2019-11-02_thd  (SHEATH_CONSISTENT)
  (b) hot-boundary false positive (keV near-MP)   : 2022-02-25_the  (HOT_BOUNDARY_FLAG)
  (c) non-depleted high-shear (southward) control : 2023-03-06_tha  (run8 CONTROL, Bz<0)
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs")  # -> repo outputs/ via config.P
SUB  = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")  # -> repo data/substrate via config.P
OUT  = os.path.dirname(__file__)

plt.rcParams.update({
    "figure.dpi":200,"savefig.dpi":200,"font.family":"serif","font.serif":["DejaVu Serif"],
    "mathtext.fontset":"dejavuserif","font.size":10,"axes.titlesize":10.5,"axes.labelsize":10,
    "axes.linewidth":0.8,"axes.grid":True,"grid.alpha":0.22,"grid.linewidth":0.5,
    "xtick.direction":"in","ytick.direction":"in","savefig.bbox":"tight","savefig.pad_inches":0.06,
})
C_DEN="#1f4e9b"; C_FLD="#b4202a"; C_BETA="#2a7f3f"; C_T="#7a3fa0"; C_S="#444444"

spec={r["eid"]:r for r in csv.DictReader(open(os.path.join(RUNS,"run12_spectral","spectral_metrics.csv")))}
cand={r["eid"]:r for r in csv.DictReader(open(os.path.join(RUNS,"run8_atlas","candidates.csv")))}

EVENTS=[
    ("2019-11-02_thd","(a) sheath-like depletion","#2a7f3f"),
    ("2022-02-25_the","(b) hot-boundary false positive","#b4202a"),
    ("2023-03-06_tha","(c) non-depleted control (southward)","#555599"),
]

def load(eid):
    d=np.load(os.path.join(SUB,eid+".npz"), allow_pickle=True)
    t=np.asarray(d["t"],float); n=np.asarray(d["n"],float); B=np.asarray(d["bmag"],float)
    beta=np.asarray(d["beta"],float); pth=np.asarray(d["p_th"],float)
    x=np.asarray(d["x_re"],float); r=np.asarray(d["r_re"],float)
    mp0=float(d["mp0"]); bs0=float(d["bs0"]); alpha=float(d["alpha"]); lam=1.17
    cth=np.clip(x/r,-1,1)
    r_mp=mp0*(2.0/(1.0+cth))**alpha
    r_bs=bs0*2.0/(cth+np.sqrt(cth**2+lam**2*(1-cth**2)))
    s=(r-r_mp)/(r_bs-r_mp)
    # T in eV from ion thermal pressure: T_eV = 6240 * p_th[nPa]/n[cm^-3]  (verified below)
    with np.errstate(divide="ignore",invalid="ignore"):
        T=6240.0*pth/n
    th=(t-t[0])/3600.0
    return dict(th=th,n=n,B=B,beta=beta,s=s,T=T,pth=pth,bmag=B,
                pb=np.asarray(d["p_b"],float),beta_store=beta)

# ---- unit verification: stored beta must equal p_th / p_b, and p_b ~ 3.98e-4*B^2 (B in nT) ----
def verify(eid):
    d=np.load(os.path.join(SUB,eid+".npz"),allow_pickle=True)
    B=np.asarray(d["bmag"],float); pb=np.asarray(d["p_b"],float); pth=np.asarray(d["p_th"],float)
    beta=np.asarray(d["beta"],float)
    m=np.isfinite(B)&np.isfinite(pb)&(B>0)&(pb>0)
    pb_pred=3.98e-4*B[m]**2          # nPa if B in nT
    rb=np.nanmedian(pb[m]/pb_pred)   # ~1 if p_b in nPa & B in nT
    mb=np.isfinite(beta)&np.isfinite(pth)&np.isfinite(pb)&(pb>0)
    rbeta=np.nanmedian((pth[mb]/pb[mb])/beta[mb])  # ~1 if beta=p_th/p_b
    return rb,rbeta

print("UNIT CHECK (expect ~1.0 each):")
for eid,_,_ in EVENTS:
    rb,rbeta=verify(eid); print(f"  {eid}: p_b/(3.98e-4 B^2)={rb:.3f}   (p_th/p_b)/beta_stored={rbeta:.3f}")

fig,axs=plt.subplots(5,3,figsize=(12.2,11.4),sharex="col",constrained_layout=True)
rows=[("$n$ (cm$^{-3}$)",C_DEN,"n",True),
      ("$|B|$ (nT)",C_FLD,"B",False),
      (r"$\beta$",C_BETA,"beta",True),
      ("$T_{\\rm ion}$ (eV)",C_T,"T",True),
      ("radial $s$",C_S,"s",False)]
for j,(eid,title,col) in enumerate(EVENTS):
    D=load(eid)
    sp=spec.get(eid); cd=cand.get(eid)
    # time masks for shading: near-MP shell (orange) and magnetosphere s<0 (purple)
    in_shell=(D["s"]>=0.05)&(D["s"]<0.20)
    in_mag=D["s"]<0.0
    def shade(ax,mask,color,alpha):
        if not mask.any(): return
        th=D["th"]; edges=np.diff(mask.astype(int))
        starts=np.where(edges==1)[0]+1; ends=np.where(edges==-1)[0]+1
        if mask[0]: starts=np.r_[0,starts]
        if mask[-1]: ends=np.r_[ends,len(th)-1]
        for a,b in zip(starts,ends): ax.axvspan(th[a],th[b],color=color,alpha=alpha,lw=0,zorder=0)
    for i,(ylab,c,key,logy) in enumerate(rows):
        ax=axs[i][j]
        ax.plot(D["th"],D[key],color=c,lw=0.9)
        if logy:
            ax.set_yscale("log")
        # shade the magnetosphere (s<0) and near-MP shell intervals (all rows)
        shade(ax,in_mag,"#d9d0ec",0.55)
        shade(ax,in_shell,"#f6c9a8",0.45)
        if key=="s":
            ax.axhspan(0.05,0.20,color="#f6c9a8",alpha=0.45,lw=0)
            ax.axhspan(-0.30,0.0,color="#d9d0ec",alpha=0.55,lw=0)
            ax.axhline(0.0,ls="--",color="#3b2f6b",lw=0.9)  # MP (s=0)
            ax.set_ylim(-0.18,0.65)
        if key=="beta":
            ax.axhline(1.0,ls=":",color="#888",lw=0.8)
        if j==0: ax.set_ylabel(ylab)
        if i==0:
            # column header with classification + committed spectral metrics
            if sp is not None:
                hdr=(f"{title}\n{eid}\n"
                     f"status: {sp['status'].replace('_',' ').lower()}\n"
                     f"peak-E ratio {float(sp['peak_ratio']):.2f} · shape corr {float(sp['shape_corr']):.2f} · flux ratio {float(sp['flux_ratio']):.2f}")
            else:
                hdr=(f"{title}\n{eid}\n"
                     f"run8 control · $B_z$={cd['bz']} nT\n"
                     f"$D_n$(mem) {float(cd['Dn_mem']):.2f} · $E_B$(mem) {float(cd['EB_mem']):.2f} (non-depleted)")
            ax.set_title(hdr,fontsize=9.2,color=col,loc="left")
    axs[-1][j].set_xlabel(f"time from window start (h)")

fig.suptitle("Three near-magnetopause event archetypes (per-encounter THEMIS time series). The ion temperature $T_{\\rm ion}$ separates\n"
             "sheath-energy ions (a) from keV boundary-layer contamination (b); (c) is a non-depleted control.   "
             "Near-MP shell $s\\in[0.05,0.20)$ shaded orange; magnetosphere $s<0$ shaded purple.",
             fontsize=10.6,fontweight="bold")
p=os.path.join(OUT,"fig7_event_archetypes.png"); fig.savefig(p); plt.close(fig); print("wrote",p)
