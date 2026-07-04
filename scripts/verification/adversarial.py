import numpy as np
import robustness as R

# ---------------------------------------------------------------------------
# Adversarial extensions:
#  B2 : extreme blunt BS (exponent x2.0 of MP -> far more flaring than Jelinek)
#  Bpar: parabolic dayside BS anchored at bs0 (Merka-style conic), r from
#        a paraboloid X = bs0 - rho^2/(2*L) solved for r at each direction.
#  Chyb: local-NORMAL distance to MP, normalised by the normal distance to the
#        BLUNT (Jelinek) BS  -> most adversarial "normal" variant.
# Goal: see whether any choice REFILLS the near bin [0.2,0.4] for P3/P4/P7 or
#       restores deep (<<1) depletion for the surviving passes.
# ---------------------------------------------------------------------------

def coord_B_exp(pos, mp0, alpha, bs0, factor):
    return R.coord_B(pos, mp0, alpha, bs0, bs_exp=alpha*factor)[0]

def coord_B_parabolic(pos, mp0, alpha, bs0):
    # Paraboloid of revolution anchored at subsolar bs0, focus on +X axis.
    # Standard form: rho^2 = 2*p*(bs0 - X)  with p the semi-latus rectum.
    # Choose p so the flank flaring is BLUNTER than Shue: use p = bs0
    # (a moderately blunt paraboloid; tested also p=1.5*bs0 below).
    p = bs0
    r = np.linalg.norm(pos, axis=1)
    costh = np.clip(pos[:,0]/r, -0.999999, 1.0)
    sinth = np.sqrt(np.clip(1-costh**2,0,1))
    r_mp = R.shue_r(mp0, alpha, costh)
    # solve for r_bs along direction (costh,sinth): rho=r*sinth, X=r*costh
    # (r*sinth)^2 = 2p(bs0 - r*costh) -> r^2 sin^2 + 2p costh r - 2p bs0 =0
    a = sinth**2
    b = 2*p*costh
    c = -2*p*bs0
    r_bs = np.where(a>1e-9,
                    (-b + np.sqrt(b**2 - 4*a*c))/(2*a),
                    bs0/np.maximum(costh,1e-9))  # subsolar limit -> bs0
    s = (r - r_mp)/(r_bs - r_mp)
    return s

def coord_C_hybrid(pos, mp0, alpha, bs0, bs_factor):
    """Normal distance to Shue MP, normalised by normal distance to a BLUNT
    Shue BS with exponent alpha*bs_factor (the Jelinek-style blunt BS)."""
    bs_exp = alpha*bs_factor
    X = pos[:,0]; rho = np.sqrt(pos[:,1]**2+pos[:,2]**2)
    n=len(X); s=np.full(n,np.nan)
    th=np.linspace(1e-3,np.pi-1e-3,4000)
    Sx,Srho=R._surface_xy(th,mp0,alpha)
    for i in range(n):
        if not(np.isfinite(X[i]) and np.isfinite(rho[i])): continue
        d2=(Sx-X[i])**2+(Srho-rho[i])**2; j=int(np.argmin(d2))
        lo=max(j-2,0); hi=min(j+3,len(th))
        thf=np.linspace(th[lo],th[hi-1],200); sxf,srf=R._surface_xy(thf,mp0,alpha)
        d2f=(sxf-X[i])**2+(srf-rho[i])**2; jf=int(np.argmin(d2f))
        px,prho=sxf[jf],srf[jf]; dist=np.sqrt(d2f[jf])
        r_sc=np.hypot(X[i],rho[i]); cth_sc=min(max(X[i]/r_sc,-0.999999),1.0)
        sign=1.0 if r_sc>=R.shue_r(mp0,alpha,cth_sc) else -1.0
        d_mp=sign*dist
        vx=X[i]-px; vrho=rho[i]-prho; vn=np.hypot(vx,vrho)
        if vn<1e-9:
            t0=max(jf-1,0); t1=min(jf+1,len(thf)-1)
            tx=sxf[t1]-sxf[t0]; trho=srf[t1]-srf[t0]; tn=np.hypot(tx,trho)
            nx,nrho=trho/tn,-tx/tn
        else:
            nx,nrho=vx/vn,vrho/vn
        if nx*px+nrho*prho<0: nx,nrho=-nx,-nrho
        def bs_resid(lmb):
            qx=px+lmb*nx; qrho=prho+lmb*nrho; rq=np.hypot(qx,qrho)
            cth=min(max(qx/rq,-0.999999),1.0) if rq>0 else 1.0
            return rq - R.shue_r(bs0,bs_exp,cth)
        lo_l=0.0; f_lo=bs_resid(lo_l); hi_l=max(2.0*bs0,50.0); f_hi=bs_resid(hi_l)
        L=np.nan
        if f_lo*f_hi<0:
            for _ in range(60):
                mid=0.5*(lo_l+hi_l); fm=bs_resid(mid)
                if f_lo*fm<=0: hi_l=mid; f_hi=fm
                else: lo_l=mid; f_lo=fm
            L=0.5*(lo_l+hi_l)
        if L and np.isfinite(L) and L>1e-6: s[i]=d_mp/L
    return s

print(f"{'Pass':4} {'SZA':>5} | "
      f"{'A':>16} | {'B_jel':>16} | {'B_x2':>16} | {'B_parab':>16} | {'Chyb_jel':>16}")
print(f"{'':4} {'':>5} |  Dn   occ  nN nB |  Dn   occ  nN nB |  Dn   occ  nN nB |  Dn   occ  nN nB |  Dn   occ  nN nB")
for name,d,sza in R.PASSES:
    rec=R.load_pass(d); t,B,dens,pos=R.build_timeline(rec)
    dp,bz,ma,mp0,alpha,bs0=R.boundary_params(rec["omni"])
    variants={
        "A":   R.coord_A(pos,mp0,alpha,bs0),
        "Bjel":R.coord_B(pos,mp0,alpha,bs0)[0],
        "Bx2": coord_B_exp(pos,mp0,alpha,bs0,2.0),
        "Bpar":coord_B_parabolic(pos,mp0,alpha,bs0),
        "Chyb":coord_C_hybrid(pos,mp0,alpha,bs0,R.JEL_BS/R.JEL_MP),
    }
    cells=[]
    for k in ["A","Bjel","Bx2","Bpar","Chyb"]:
        Dn,nN,nB,nT,occ=R.metrics(variants[k],dens)
        ds=f"{Dn:5.3f}" if np.isfinite(Dn) else "  nan"
        cells.append(f"{ds} {occ:4.2f} {nN:3d}{nB:3d}")
    print(f"{name:4} {sza:5.1f} | " + " | ".join(cells))
