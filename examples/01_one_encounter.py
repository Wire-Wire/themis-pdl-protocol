'''Example 1: analyse ONE encounter end-to-end with the library (shipped demo data).'''
import os, sys
import numpy as np
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol import load_encounter, member_mask, shell_contrast
from pdl_protocol.config import SUBSTRATE

d = np.load(os.path.join(SUBSTRATE, "2020-12-12_tha.npz"), allow_pickle=True)
e = load_encounter(d)
m = member_mask(e, TBAND=3.0, VFRAC=0.2, NFLOOR=0.3)
dn, eb = shell_contrast(e, m, 0.05, 0.20)
print(f"THEMIS-A 2020-12-12: D_n = {dn:.3f}, E_B = {eb:.3f}")
print("committed values: 0.554 / 2.279 -> a depleted-density, enhanced-field candidate interval")
