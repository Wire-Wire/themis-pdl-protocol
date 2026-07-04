import numpy as np
import robustness as R

# Final probes:
# (1) Can ANY BS exponent factor (0.3..3.0) refill the near bin for P3/P4/P7
#     or drive a surviving pass to deep depletion (Dn<0.5)? Sweep it.
# (2) Occupancy robustness: report near-occupancy with TWO denominators:
#     (a) in-band 0<=s<1 (used in report)  (b) all finite-s samples.

print("=== (1) BS-exponent sweep: factor = bs_exp/alpha ===")
print("    For each pass show whether near bin [0.2,0.4] is occupied (occ) and Dn.")
factors=[0.3,0.5,0.8,1.0,1.5,2.0,3.0]
for name,d,sza in R.PASSES:
    rec=R.load_pass(d); t,B,dens,pos=R.build_timeline(rec)
    dp,bz,ma,mp0,alpha,bs0=R.boundary_params(rec["omni"])
    line=f"{name} SZA{sza:4.1f}: "
    for f in factors:
        s=R.coord_B(pos,mp0,alpha,bs0,bs_exp=alpha*f)[0]
        Dn,nN,nB,nT,occ=R.metrics(s,dens)
        ds=f"{Dn:4.2f}" if np.isfinite(Dn) else "nan"
        line+=f" [x{f:.1f}:{ds},occ{occ:.2f}]"
    print(line)

print("\n=== (2) near-occupancy under two denominators (coord A) ===")
print(f"{'Pass':4} {'occ_inband':>11} {'occ_allfinite':>14}  nN  nFiniteS  nInBand")
for name,d,sza in R.PASSES:
    rec=R.load_pass(d); t,B,dens,pos=R.build_timeline(rec)
    dp,bz,ma,mp0,alpha,bs0=R.boundary_params(rec["omni"])
    s=R.coord_A(pos,mp0,alpha,bs0)
    ok=np.isfinite(s)&np.isfinite(dens); ss=s[ok]
    near=((ss>=0.2)&(ss<0.4)).sum()
    inband=((ss>=0.0)&(ss<1.0)).sum()
    occ_in=near/inband if inband else np.nan
    occ_all=near/len(ss) if len(ss) else np.nan
    print(f"{name:4} {occ_in:11.3f} {occ_all:14.3f} {near:4d} {len(ss):8d} {inband:8d}")
