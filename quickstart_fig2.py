"""Quickstart: reproduce the headline figure (paired coordinate test) from committed_outputs/ in seconds."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "committed_outputs", "run23_paired", "paired_1d_vs_radial.csv"), newline="") as f:
    rows = list(csv.DictReader(f))
d1 = np.array([float(r["Dn_1d"]) for r in rows]); dr = np.array([float(r["Dn_radial"]) for r in rows])
fig, ax = plt.subplots(figsize=(6, 4.2))
bp = ax.boxplot([d1, dr], labels=["fixed Sun-Earth axis", "radial"], showfliers=False, widths=0.5)
ax.axhline(1.0, color="0.4", ls="--", lw=1)
ax.set_ylabel(r"near-boundary density contrast $D_n$")
ax.set_title(f"Same {len(rows)} encounters, two coordinate constructions:\n"
             f"median {np.median(d1):.3f} -> {np.median(dr):.3f} (the apparent depletion is a coordinate artefact)")
out = os.path.join(HERE, "outputs"); os.makedirs(out, exist_ok=True)
fig.tight_layout(); fig.savefig(os.path.join(out, "quickstart_fig2.png"), dpi=150)
print("wrote outputs/quickstart_fig2.png   (compare with figures/fig2_paired_artefact.png)")
