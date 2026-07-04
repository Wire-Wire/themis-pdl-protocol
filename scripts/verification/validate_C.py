import numpy as np
import robustness as R

# Validate coord C internals on a synthetic point well off-axis, where radial
# and normal MUST differ, plus instrument the real passes for fallback usage.

mp0, alpha, bs0 = 10.0, 0.58, 13.0

# ---- synthetic check: a point at large theta (flank) ----
# place a spacecraft midway (radially) between MP and BS at theta=70deg
th0 = np.radians(70.0)
r_mp = R.shue_r(mp0, alpha, np.cos(th0))
r_bs = R.shue_r(bs0, alpha, np.cos(th0))
r_sc = 0.5*(r_mp+r_bs)
X = r_sc*np.cos(th0); rho = r_sc*np.sin(th0)
pos = np.array([[X, rho, 0.0]])
sC = R.coord_C(pos, mp0, alpha, bs0)
# radial s for comparison
sA = R.coord_A(pos, mp0, alpha, bs0)
print("FLANK synthetic theta=70deg:")
print(f"  r_mp(rad)={r_mp:.3f}  r_bs(rad)={r_bs:.3f}  r_sc={r_sc:.3f}")
print(f"  s_radial(A) = {sA[0]:.4f}   s_normal(C) = {sC[0]:.4f}   diff={sA[0]-sC[0]:+.4f}")
print("  (at a true flank these SHOULD differ noticeably -> proves C != A structurally)\n")

# instrument: re-implement the normal-direction calc to report angle between
# the local MP normal and the radial direction at each real spacecraft sample.
def normal_vs_radial_angles(pos, mp0, alpha):
    X = pos[:,0]; rho = np.sqrt(pos[:,1]**2+pos[:,2]**2)
    th = np.linspace(1e-3, np.pi-1e-3, 4000)
    Sx, Srho = R._surface_xy(th, mp0, alpha)
    ang=[]
    for i in range(len(X)):
        if not (np.isfinite(X[i]) and np.isfinite(rho[i])): continue
        d2=(Sx-X[i])**2+(Srho-rho[i])**2; j=int(np.argmin(d2))
        px,prho=Sx[j],Srho[j]
        # outward normal from gradient of surface curve
        t0=max(j-1,0); t1=min(j+1,len(th)-1)
        tx=Sx[t1]-Sx[t0]; trho=Srho[t1]-Srho[t0]
        tn=np.hypot(tx,trho); nx,nrho=trho/tn,-tx/tn
        if nx*px+nrho*prho<0: nx,nrho=-nx,-nrho
        # radial dir at footpoint
        rfp=np.hypot(px,prho); rx,rrho=px/rfp,prho/rfp
        cosang=np.clip(nx*rx+nrho*rrho,-1,1)
        ang.append(np.degrees(np.arccos(cosang)))
    return np.array(ang)

print("Per-pass: angle between local MP normal and radial at footpoints (deg)")
for name,d,sza in R.PASSES:
    rec=R.load_pass(d); t,B,dens,pos=R.build_timeline(rec)
    dp,bz,ma,mp0p,alpha,bs0=R.boundary_params(rec["omni"])
    a=normal_vs_radial_angles(pos,mp0p,alpha)
    print(f"  {name} SZA{sza:4.1f}: normal-radial angle min{np.nanmin(a):4.1f} med{np.nanmedian(a):4.1f} max{np.nanmax(a):4.1f}")
