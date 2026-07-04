"""Data acquisition. Modes:
  --check           : report what is present (demo subset ships with the repo)
  --event EID       : fetch 1-min OMNI + the ESA ion spectrogram for one encounter (e.g. 2020-12-12_tha)
  --omni            : fetch 1-min OMNI for every encounter present in data/substrate (feeds run27/run29)
  --full-substrate  : rebuild the full 6,248-encounter substrate from CDAWeb (hours; see docs/DATA.md)
Bring-your-own-data: drop NPZ files matching the schema in docs/DATA.md into data/substrate/.
"""
import argparse, os, sys, glob
import numpy as np

try:
    import pdl_protocol  # installed: pip install -e .
except ImportError:      # source checkout without install
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "reproduce"))
from pdl_protocol.config import SUBSTRATE, OMNI_CACHE, PIPE, EVENTS

ap = argparse.ArgumentParser()
ap.add_argument("--check", action="store_true")
ap.add_argument("--event")
ap.add_argument("--omni", action="store_true")
ap.add_argument("--full-substrate", action="store_true")
a = ap.parse_args()

if a.check or not any([a.event, a.omni, a.full_substrate]):
    subs = glob.glob(os.path.join(SUBSTRATE, "*.npz"))
    print(f"substrate encounters present: {len(subs)} in {SUBSTRATE}")
    for s in subs[:10]:
        print("  ", os.path.basename(s))
    print(f"omni cache: {len(glob.glob(os.path.join(OMNI_CACHE, '*.npz')))} windows")
    sys.exit(0)

if a.event:
    from datetime import datetime, timezone
    from cdasws import CdasWs
    cdas = CdasWs()
    eid = a.event; date, probe = eid.rsplit("_", 1)
    d = np.load(os.path.join(SUBSTRATE, eid + ".npz"), allow_pickle=True)
    t0, t1 = float(d["t"][0]), float(d["t"][-1])
    iso = lambda ts: datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    st, om = cdas.get_data("OMNI_HRO_1MIN", ["Pressure", "BZ_GSM", "flow_speed", "proton_density"], iso(t0), iso(t1))
    np.savez_compressed(os.path.join(EVENTS, f"omni_{date}.npz"),
                        Pressure=np.array(om["Pressure"], float), BZ=np.array(om["BZ_GSM"], float))
    var = f"{probe}_peir_en_eflux"
    st, esa = cdas.get_data(f"{probe.upper()}_L2_ESA", [var], iso(t0), iso(t1))
    da = esa[var]
    tunix = esa[da.dims[0]].values.astype("datetime64[s]").astype(float)
    keep = (tunix >= t0) & (tunix <= t1)
    np.savez_compressed(os.path.join(EVENTS, f"esa_{date}.npz"), eflux=np.array(da.values, float)[keep],
                        tunix=tunix[keep], energy=np.logspace(np.log10(5.6), np.log10(25300.0), 32))
    print(f"fetched omni_{date}.npz + esa_{date}.npz into data/events/ (used by make_fig8_deepdive.py)")
elif a.omni:
    sys.argv = [sys.argv[0]]
    import fetch_omni_contributing as F  # P-wrapped: reads data/substrate, writes data/omni_cache
    F.main()
elif a.full_substrate:
    print("Full rebuild: python pipeline/build_substrate.py  (fetches THEMIS FGM/ESA + OMNI for 2007-2025;")
    print("takes hours and ~2 GB; see docs/DATA.md for the encounter definition and schema)")
