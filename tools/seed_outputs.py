"""Seed outputs/ from committed_outputs/ so any pipeline stage can run without recomputing its inputs."""
import os, shutil, sys
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol.config import OUTPUTS, COMMITTED
n = 0
for root, _, files in os.walk(COMMITTED):
    rel = os.path.relpath(root, COMMITTED)
    for f in files:
        d = os.path.join(OUTPUTS, rel)
        os.makedirs(d, exist_ok=True)
        shutil.copy(os.path.join(root, f), d)
        n += 1
print(f"seeded {n} committed files into outputs/")
