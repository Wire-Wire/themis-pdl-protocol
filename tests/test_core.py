'''Unit tests anchored to committed values; runs on the shipped 3-encounter demo subset.'''
import os, sys
import numpy as np
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol import load_encounter, member_mask, shell_contrast, compute_s, list_files

EXPECT = {"2020-12-12_tha": (0.554, 2.279), "2011-12-18_the": (0.526, 1.998)}
fails = 0
files = list_files()
assert len(files) >= 3, "demo substrate missing"
for f in files:
    eid = os.path.basename(f)[:-4]
    d = np.load(f, allow_pickle=True)
    e = load_encounter(d)
    assert e is not None, eid
    s = compute_s(d["x_re"], d["y_re"], d["z_re"], float(d["mp0"]), float(d["alpha"]), float(d["dp"]))
    assert np.nanmin(s) < 0.2 and np.nanmax(s) > 0.6, "s range unexpected for " + eid
    if eid in EXPECT:
        dn, eb = shell_contrast(e, member_mask(e, 3.0, 0.2, 0.3), 0.05, 0.20)
        edn, eeb = EXPECT[eid]
        ok = abs(dn - edn) < 0.01 and abs(eb - eeb) < 0.01
        print(("PASS " if ok else "FAIL ") + eid + f": D_n={dn:.3f} (exp {edn}), E_B={eb:.3f} (exp {eeb})")
        fails += (not ok)
print("test_core:", "ALL PASS" if not fails else str(fails) + " FAIL")
sys.exit(1 if fails else 0)
