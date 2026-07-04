'''Example 5: the one-command contrast analysis on the shipped demo encounters.

Equivalent CLI:  pdl-contrast --paired-check

On three demo encounters the statistics are illustrative only; point the same calls
at a full substrate (docs/TUTORIAL.md level 6) for population-level results.
'''
import os, sys
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol import shell_profile, population_contrast, paired_coordinate_check

for r in shell_profile():
    print("s [%.2f, %.2f)  N=%d  D_n=%.3f  E_B=%.3f" % (r["s_lo"], r["s_hi"], r["N"], r["D_n"], r["E_B"]))
print("population:", population_contrast())
print("paired coordinate check:", paired_coordinate_check())
