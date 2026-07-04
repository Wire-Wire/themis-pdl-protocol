'''Example 4: the one-command finder on the shipped demo encounters.

Equivalent CLI:  pdl-find            (after pip install -e .)
'''
import os, sys
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol import find_candidates

rows = find_candidates()
for r in rows:
    print("%-18s D_n=%.3f E_B=%.3f N=%3d  %s" % (r["eid"], r["D_n"], r["E_B"], r["n_samples"], r["status"]))
n = sum(r["candidate"] for r in rows)
print("-> %d of %d encounters are candidates (selection-limited indications, not confirmed layers)"
      % (n, len(rows)))
