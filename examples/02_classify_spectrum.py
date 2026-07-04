'''Example 2: classify the shipped ion spectrogram of the deep-dive event.'''
import os, sys
import numpy as np
try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from pdl_protocol.spectral import spectral_metrics, classify_spectrum
from pdl_protocol.config import EVENTS, SUBSTRATE
from pdl_protocol import compute_s

esa = np.load(os.path.join(EVENTS, "esa_2020-12-12.npz"))
d = np.load(os.path.join(SUBSTRATE, "2020-12-12_tha.npz"), allow_pickle=True)
s_sub = compute_s(d["x_re"], d["y_re"], d["z_re"], float(d["mp0"]), float(d["alpha"]), float(d["dp"]))
s_at = np.interp(esa["tunix"], d["t"].astype(float), s_sub)
near = np.nanmedian(esa["eflux"][(s_at >= 0.05) & (s_at < 0.20)], axis=0)
bg = np.nanmedian(esa["eflux"][(s_at >= 0.6) & (s_at <= 1.0)], axis=0)
peak, shape, flux = spectral_metrics(near, bg, esa["energy"])
print(f"peak ratio {peak:.2f}, shape corr {shape:.2f}, flux ratio {flux:.2f} -> {classify_spectrum(peak, shape)}")
print("(pipeline-committed metrics: 1.00 / 0.96 / 0.39; the example's simplified shell windows give a slightly higher flux ratio - the classification is identical)")
