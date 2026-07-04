import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.dirname(_cfg_os.path.abspath(__file__))))
from config import P
import numpy as np, os, sys

BASE = P(r"H:\0mssl\0mssl519\data_cache\normalized")
RE = 6371.2

PASSES = [
    ("P1", "usable_aug18_6h",     21.9),
    ("P2", "usable_sep03_6h",     13.8),
    ("P3", "usable_sep13_09_6h",  16.9),
    ("P4", "usable_sep20_09_6h",  10.4),
    ("P5", "usable_sep26_09_10h",  4.4),
    ("P6", "usable_sep27_09_10h",  3.8),
    ("P7", "usable_oct24_09_6h",  20.4),
]

# ---------------------------------------------------------------------------
# Loading / timeline construction
# ---------------------------------------------------------------------------
def _epoch_s(dt64):
    return (dt64.astype("datetime64[ns]").astype("int64")) / 1e9

def load_pass(d):
    p = os.path.join(BASE, d)
    fgm = np.load(os.path.join(p, "thd_l2_fgm.npz"), allow_pickle=True)
    mom = np.load(os.path.join(p, "thd_l2_mom.npz"), allow_pickle=True)
    st  = np.load(os.path.join(p, "thd_l1_state.npz"), allow_pickle=True)
    omni= np.load(os.path.join(p, "omni_hro_1min.npz"), allow_pickle=True)

    t_fgm = _epoch_s(fgm["_time"]); Bvec = np.asarray(fgm["thd_fgs_gsm"], float)
    Bmag  = np.linalg.norm(Bvec, axis=1)
    t_mom = _epoch_s(mom["_time"]); dens = np.asarray(mom["thd_peim_density"], float)
    t_st  = _epoch_s(st["_time"]);  pos  = np.asarray(st["thd_pos_gsm"], float) / RE  # Re
    return dict(t_fgm=t_fgm, Bmag=Bmag, t_mom=t_mom, dens=dens,
                t_st=t_st, pos=pos, omni=omni)

def interp_gap(t_new, t_src, y_src, gap=60.0):
    """Linear interp; set to NaN where nearest source sample is > gap away."""
    m = np.isfinite(y_src)
    t_src = t_src[m]; y_src = y_src[m]
    order = np.argsort(t_src); t_src=t_src[order]; y_src=y_src[order]
    y = np.interp(t_new, t_src, y_src, left=np.nan, right=np.nan)
    idx = np.searchsorted(t_src, t_new)
    idx = np.clip(idx, 1, len(t_src)-1)
    dl = np.abs(t_new - t_src[idx-1]); dr = np.abs(t_src[idx] - t_new)
    nearest = np.minimum(dl, dr)
    y[nearest > gap] = np.nan
    return y

def build_timeline(rec):
    t0 = max(rec["t_fgm"][0], rec["t_mom"][0], rec["t_st"][0])
    t1 = min(rec["t_fgm"][-1], rec["t_mom"][-1], rec["t_st"][-1])
    t = np.arange(np.ceil(t0), np.floor(t1)+1e-9, 10.0)
    Bmag = interp_gap(t, rec["t_fgm"], rec["Bmag"])
    dens = interp_gap(t, rec["t_mom"], rec["dens"])
    px = interp_gap(t, rec["t_st"], rec["pos"][:,0])
    py = interp_gap(t, rec["t_st"], rec["pos"][:,1])
    pz = interp_gap(t, rec["t_st"], rec["pos"][:,2])
    dens = np.where(dens >= 0.01, dens, np.nan)  # clip density >= 0.01
    pos = np.column_stack([px,py,pz])
    return t, Bmag, dens, pos

# ---------------------------------------------------------------------------
# Boundary parameters from OMNI medians
# ---------------------------------------------------------------------------
def boundary_params(omni):
    dp = np.nanmedian(np.asarray(omni["Pressure"], float))
    bz = np.nanmedian(np.asarray(omni["BZ_GSM"], float))
    ma = np.nanmedian(np.asarray(omni["Mach_num"], float))
    mp0 = (11.4 + (0.013 if bz >= 0 else 0.140) * bz) * dp ** (-1.0/6.6)
    alpha = (0.58 - 0.010 * bz) * (1.0 + 0.010 * dp)
    if ma > 1:
        bs0 = mp0 * (1.0 + 1.1 * (((2.0/3.0)*ma**2 + 2.0) / ((8.0/3.0)*ma**2)))
    else:
        bs0 = 1.3 * mp0
    return dp, bz, ma, mp0, alpha, bs0

# ---------------------------------------------------------------------------
# Shue surface helpers
# ---------------------------------------------------------------------------
def shue_r(r0, alpha, costh):
    return r0 * (2.0/(1.0+costh))**alpha

# ---------------------------------------------------------------------------
# Coordinate A: radial, BS flared with same alpha
# ---------------------------------------------------------------------------
def coord_A(pos, mp0, alpha, bs0):
    r = np.linalg.norm(pos, axis=1)
    costh = np.clip(pos[:,0]/r, -0.999999, 1.0)
    r_mp = shue_r(mp0, alpha, costh)
    r_bs = shue_r(bs0, alpha, costh)
    s = (r - r_mp)/(r_bs - r_mp)
    return s

# ---------------------------------------------------------------------------
# Coordinate B: radial, proper blunt bow shock (Jelinek-2012 style)
#   Jelinek et al. 2012 give MP and BS both as power-law (Shue-type) surfaces
#   r = R0 * (2/(1+cos th))^lambda, with lambda_BS ~ 1.17, lambda_MP ~ 0.78.
#   i.e. BS flares MORE than MP. We anchor BS subsolar at bs0 (as required),
#   keep MP = Shue with the pass's alpha, and give the BS its OWN larger
#   exponent. To make the BS strictly blunter than the MP for every pass we
#   scale the BS exponent by the empirical Jelinek ratio (1.17/0.78 ~ 1.5).
# ---------------------------------------------------------------------------
JEL_BS = 1.17
JEL_MP = 0.78
def coord_B(pos, mp0, alpha, bs0, bs_exp=None):
    r = np.linalg.norm(pos, axis=1)
    costh = np.clip(pos[:,0]/r, -0.999999, 1.0)
    r_mp = shue_r(mp0, alpha, costh)
    if bs_exp is None:
        bs_exp = alpha * (JEL_BS/JEL_MP)   # blunter BS, flares more
    r_bs = shue_r(bs0, bs_exp, costh)
    s = (r - r_mp)/(r_bs - r_mp)
    return s, bs_exp

# ---------------------------------------------------------------------------
# Coordinate C: local-normal distance to the Shue MP surface.
#   For each spacecraft point find the nearest point on the Shue MP surface
#   r_mp(theta)=mp0*(2/(1+cos th))^alpha (axisymmetric about the X axis).
#   Work in the meridian plane (X, rho) where rho=sqrt(y^2+z^2): the surface
#   is a 1-D curve, nearest point reduces to a 1-D minimisation over theta.
#   d_mp = signed distance (positive outside MP). Then continue along the SAME
#   outward unit normal until it crosses the BS surface (radial-Shue BS with
#   alpha, i.e. the boundary set of coord A) to get the local MP->BS normal
#   thickness L. s_normal = d_mp / L.
# ---------------------------------------------------------------------------
def _surface_xy(theta, r0, alpha):
    # point on axisymmetric surface in (X, rho) meridian plane
    costh = np.cos(theta)
    rr = r0 * (2.0/(1.0+costh))**alpha
    return rr*costh, rr*np.sin(theta)   # X, rho

def coord_C(pos, mp0, alpha, bs0):
    X = pos[:,0]
    rho = np.sqrt(pos[:,1]**2 + pos[:,2]**2)
    n = len(X)
    s = np.full(n, np.nan)

    # dense theta grid for the MP curve (avoid theta=pi singularity)
    th = np.linspace(1e-3, np.pi-1e-3, 4000)
    Sx, Srho = _surface_xy(th, mp0, alpha)

    for i in range(n):
        if not (np.isfinite(X[i]) and np.isfinite(rho[i])):
            continue
        # nearest point on MP curve to (X[i], rho[i])
        d2 = (Sx - X[i])**2 + (Srho - rho[i])**2
        j = int(np.argmin(d2))
        # refine with local quadratic / fine search around j
        lo = max(j-2,0); hi = min(j+3,len(th))
        thf = np.linspace(th[lo], th[hi-1], 200)
        sxf, srf = _surface_xy(thf, mp0, alpha)
        d2f = (sxf - X[i])**2 + (srf - rho[i])**2
        jf = int(np.argmin(d2f))
        px, prho = sxf[jf], srf[jf]
        dist = np.sqrt(d2f[jf])
        vx = X[i]-px; vrho = rho[i]-prho
        # sign: positive if spacecraft is OUTSIDE the MP. Shue is star-shaped
        # from the origin, so the robust inside/outside test is: compare the
        # spacecraft radius to the Shue radius at the spacecraft's OWN polar
        # angle (identical to coord A's inside/outside notion).
        r_sc = np.hypot(X[i], rho[i])
        cth_sc = X[i]/r_sc if r_sc > 0 else 1.0
        cth_sc = min(max(cth_sc, -0.999999), 1.0)
        sign = 1.0 if r_sc >= shue_r(mp0, alpha, cth_sc) else -1.0
        d_mp = sign*dist

        # local outward unit normal = direction of (vx,vrho) (points away from
        # surface toward sc). If sc essentially on surface, fall back to
        # gradient normal.
        vn = np.hypot(vx, vrho)
        if vn < 1e-9:
            # numerical normal from curve tangent
            t0 = max(jf-1,0); t1=min(jf+1,len(thf)-1)
            tx = sxf[t1]-sxf[t0]; trho = srf[t1]-srf[t0]
            tn = np.hypot(tx,trho); nx,nrho = trho/tn, -tx/tn
            # orient outward (away from origin)
            if nx*px+nrho*prho < 0: nx,nrho=-nx,-nrho
        else:
            nx, nrho = vx/vn, vrho/vn
        # ensure outward-pointing (component along footpoint radial >0)
        if nx*px + nrho*prho < 0:
            nx, nrho = -nx, -nrho

        # march outward along normal from footpoint to find BS crossing.
        # BS (coord-A definition): radial Shue with alpha, r_bs(thetaB)=bs0*(2/(1+cos))^a
        # Find lambda>0 s.t. footpoint + lambda*nhat lies on BS surface.
        def bs_resid(lmb):
            qx = px + lmb*nx; qrho = prho + lmb*nrho
            rq = np.hypot(qx, qrho)
            cth = qx/rq if rq>0 else 1.0
            cth = min(max(cth,-0.999999),1.0)
            return rq - shue_r(bs0, alpha, cth)
        # bracket
        lo_l, hi_l = 0.0, 1.0
        f_lo = bs_resid(lo_l)
        # expand hi until sign change or large
        hi_l = max(2.0*bs0, 50.0)
        f_hi = bs_resid(hi_l)
        L = np.nan
        if f_lo == 0:
            L = 0.0
        elif f_lo*f_hi < 0:
            for _ in range(60):
                mid=0.5*(lo_l+hi_l); fm=bs_resid(mid)
                if f_lo*fm<=0: hi_l=mid; f_hi=fm
                else: lo_l=mid; f_lo=fm
            L=0.5*(lo_l+hi_l)
        else:
            # fallback: radial thickness at footpoint angle
            cthp = px/r_fp if r_fp>0 else 1.0
            cthp=min(max(cthp,-0.999999),1.0)
            L = shue_r(bs0,alpha,cthp)-shue_r(mp0,alpha,cthp)
        if L and np.isfinite(L) and L>1e-6:
            s[i] = d_mp / L
    return s

# ---------------------------------------------------------------------------
# Dn and near-occupancy
# ---------------------------------------------------------------------------
def metrics(s, dens):
    ok = np.isfinite(s) & np.isfinite(dens)
    s=s[ok]; dens=dens[ok]
    near = (s>=0.2)&(s<0.4)
    bg   = (s>=0.6)&(s<1.0)
    n_near = int(near.sum()); n_bg=int(bg.sum())
    n_tot = int(((s>=0.0)&(s<1.0)).sum())  # samples inside sheath band
    if n_near>0 and n_bg>0:
        Dn = np.median(dens[near])/np.median(dens[bg])
    else:
        Dn = np.nan
    occ = n_near/ n_tot if n_tot>0 else np.nan
    return Dn, n_near, n_bg, n_tot, occ

# ---------------------------------------------------------------------------
def main():
    print(f"{'Pass':4} {'SZA':>5} {'dp':>5} {'bz':>6} {'ma':>5} {'mp0':>6} {'alf':>5} {'bs0':>6} {'bsexpB':>6}")
    rows={}
    for name, d, sza in PASSES:
        rec=load_pass(d)
        t,Bmag,dens,pos=build_timeline(rec)
        dp,bz,ma,mp0,alpha,bs0=boundary_params(rec["omni"])
        sA=coord_A(pos,mp0,alpha,bs0)
        sB,bsexp=coord_B(pos,mp0,alpha,bs0)
        sC=coord_C(pos,mp0,alpha,bs0)
        rows[name]=dict(sza=sza,dp=dp,bz=bz,ma=ma,mp0=mp0,alpha=alpha,bs0=bs0,bsexp=bsexp,
                        sA=sA,sB=sB,sC=sC,dens=dens)
        print(f"{name:4} {sza:5.1f} {dp:5.2f} {bz:6.2f} {ma:5.1f} {mp0:6.2f} {alpha:5.3f} {bs0:6.2f} {bsexp:6.3f}")

    print("\n=== Dn and near-occupancy by coordinate ===")
    hdr=f"{'Pass':4} {'SZA':>5} | {'Dn_A':>6} {'occA':>5} {'nN':>4}{'nB':>4} | {'Dn_B':>6} {'occB':>5} {'nN':>4}{'nB':>4} | {'Dn_C':>6} {'occC':>5} {'nN':>4}{'nB':>4}"
    print(hdr)
    summary={}
    for name,_,sza in PASSES:
        R=rows[name]
        out=[]
        m={}
        for key in ["sA","sB","sC"]:
            Dn,nN,nB,nT,occ=metrics(R[key],R["dens"])
            m[key]=(Dn,occ,nN,nB)
            out.append((Dn,occ,nN,nB))
        summary[name]=m
        def fmt(x):
            Dn,occ,nN,nB=x
            ds=f"{Dn:6.3f}" if np.isfinite(Dn) else "   nan"
            os_=f"{occ:5.2f}" if np.isfinite(occ) else "  nan"
            return f"{ds} {os_} {nN:4d}{nB:4d}"
        print(f"{name:4} {sza:5.1f} | {fmt(out[0])} | {fmt(out[1])} | {fmt(out[2])}")

    # quantify B,C vs A
    print("\n=== Differences vs A ===")
    print(f"{'Pass':4} | {'dDn(B-A)':>9} {'dDn(C-A)':>9} | {'docc(B-A)':>9} {'docc(C-A)':>9}")
    for name,_,_ in PASSES:
        m=summary[name]
        DnA,occA,_,_=m["sA"]; DnB,occB,_,_=m["sB"]; DnC,occC,_,_=m["sC"]
        def df(a,b):
            return (b-a) if (np.isfinite(a) and np.isfinite(b)) else np.nan
        print(f"{name:4} | {df(DnA,DnB):9.3f} {df(DnA,DnC):9.3f} | {df(occA,occB):9.3f} {df(occA,occC):9.3f}")

if __name__=="__main__":
    main()
