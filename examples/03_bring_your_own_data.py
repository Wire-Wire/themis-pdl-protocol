'''Example 3: run the protocol on YOUR data (synthetic encounter matching docs/DATA.md).'''
import os, sys
import numpy as np
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol import load_encounter, member_mask, shell_contrast

rng = np.random.default_rng(1)
N = 2000
x = np.linspace(10.2, 13.0, N); y = np.zeros(N); z = np.zeros(N)
n = 20 + rng.normal(0, 1, N); b = 20 + rng.normal(0, 1, N); v = np.full(N, 100.0)
inner = x < 10.8
n[inner] *= 0.5
b[inner] *= 2.0
KB = 1.602e-4
pth = n * KB * 200.0
d = dict(t=np.arange(N, dtype=float), x_re=x, y_re=y, z_re=z, n=n, bmag=b, vmag=v,
         beta=pth / (3.98e-4 * b ** 2), p_th=pth, p_b=3.98e-4 * b ** 2, p_tot=pth,
         dp=np.float64(2.0), bz=np.float64(1.0), mp0=np.float64(10.0), alpha=np.float64(0.58),
         bs0=np.float64(13.5), ma=np.float64(8.0), cone_deg=np.float64(70.0),
         clock_deg=np.float64(20.0), sza_deg=np.float64(5.0), r_re=np.float64(11.5),
         x_mean=np.float64(11.5), probe="syn", date="2026-01-01")
e = load_encounter(d)
dn, eb = shell_contrast(e, member_mask(e, 3.0, 0.2, 0.3), 0.05, 0.20)
print(f"synthetic encounter: D_n = {dn:.2f}, E_B = {eb:.2f}  (injected: 0.5 / 2.0)")
print("-> drop real NPZ files with this schema into data/substrate/ and every script runs on them")
