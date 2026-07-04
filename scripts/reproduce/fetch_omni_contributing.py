"""Batch-fetch OMNI 1-min time series for the 661 contributing encounters (review).

Feeds TWO analyses with one fetch:
  run27 — within-window boundary-motion QC (reviewer request: "check the MP/BS positions
          do not change significantly during the interval; exclude such events")
  run29 — THEMIS-to-OMNI density-ratio sanity check (the reviewer's in-situ-to-OMNI
          region-classification check)

For each contributing eid: window = the encounter's own substrate time range; fetch
OMNI_HRO_1MIN [Pressure, BZ_GSM, flow_speed, proton_density]; save omni_cache/{eid}.npz.
Resumable (skips existing). Fill values -> NaN.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys
from datetime import datetime, timezone
import numpy as np
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
CACHE = P(r"H:\0mssl\review\repair\option3\omni_cache")
CONTRIB = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run10_selection\funnel_contributing.csv")
os.makedirs(CACHE, exist_ok=True)

import csv
with open(CONTRIB, newline="") as f:
    EIDS = [r["eid"] for r in csv.DictReader(f)]

from cdasws import CdasWs
cdas = CdasWs()
VARS = ["Pressure", "BZ_GSM", "flow_speed", "proton_density"]
FILL = {"Pressure": 90.0, "BZ_GSM": 9000.0, "flow_speed": 90000.0, "proton_density": 900.0}


def iso(ts):
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fetch_one(eid):
    out = os.path.join(CACHE, eid + ".npz")
    if os.path.exists(out):
        return "skip"
    try:
        d = np.load(os.path.join(SUB, eid + ".npz"), allow_pickle=True)
        t = d["t"].astype(float)
        status, data = cdas.get_data("OMNI_HRO_1MIN", VARS, iso(t[0]), iso(t[-1]))
        if data is None:
            return "nodata"
        tc = data[VARS[0]].dims[0]
        tunix = data[tc].values.astype("datetime64[s]").astype(float)
        arrs = {"tunix": tunix}
        for v in VARS:
            a = np.array(data[v].values, dtype=float)
            a[np.abs(a) >= FILL[v]] = np.nan
            arrs[v] = a
        np.savez_compressed(out, **arrs)
        return "ok"
    except Exception as ex:
        return f"err:{type(ex).__name__}"


def main():
    print(f"fetching OMNI for {len(EIDS)} contributing encounters -> {CACHE}", flush=True)
    counts = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        for i, r in enumerate(ex.map(fetch_one, EIDS)):
            counts[r] = counts.get(r, 0) + 1
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(EIDS)}  {counts}", flush=True)
    print(f"DONE  {counts}", flush=True)


if __name__ == "__main__":
    main()
