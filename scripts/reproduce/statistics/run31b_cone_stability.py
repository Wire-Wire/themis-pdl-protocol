"""run31b - within-window OMNI cone-angle stability, for the cone-preference robustness re-check.

The surviving environmental result (A cone 69.3 vs NN 63.3, run28; survives the run30 null) rests
on OMNI-propagated IMF DIRECTION, which carries propagation error (Vokhmyanin 2019: prediction
quality degrades for transverse-distant monitors). The frozen omni_cache stored only Pressure/BZ/
flow/density - NOT Bx/By - so within-window cone(t) was never available. Here we re-fetch the IMF
vector over each encounter's EXACT cached time window (tunix from omni_cache, so the windows match
run27/run28), compute cone(t) = arccos(|Bx|/|B|), and record the per-encounter median and IQR.

Output feeds run31's "does the cone preference survive on the upstream-steady subset?" test.
Resumable: skips eids already in the output CSV. Continue-on-failure with coverage reporting.
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, time, socket
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
socket.setdefaulttimeout(30)  # so a stalled OMNI call raises instead of hanging the whole run
CACHE = P(r"H:\0mssl\review\repair\option3\omni_cache")
R28 = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run28_candidate_conditions\run28_per_encounter.csv")
OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run31_candidate_context")
os.makedirs(OUT, exist_ok=True)
CSV = os.path.join(OUT, "run31_cone_stability.csv")

with open(R28, newline="") as f:
    recs = list(csv.DictReader(f))
# universe the cone test lives in: northward contributing encounters
north = [r for r in recs if r["north"] == "True"]
print(f"northward contributing in run28: {len(north)}", flush=True)

done = {}
if os.path.exists(CSV):
    with open(CSV, newline="") as f:
        for r in csv.DictReader(f):
            done[r["eid"]] = r
    print(f"resume: {len(done)} already fetched", flush=True)

from cdasws import CdasWs
cdas = CdasWs()


def fetch_cone(t0, t1):
    """median, iqr, n of cone(t)=acos(|Bx|/|B|) in deg over [t0,t1] unix seconds."""
    import datetime as dt
    s = dt.datetime.utcfromtimestamp(t0).strftime("%Y-%m-%dT%H:%M:%SZ")
    e = dt.datetime.utcfromtimestamp(t1).strftime("%Y-%m-%dT%H:%M:%SZ")
    st, data = cdas.get_data("OMNI_HRO_1MIN", ["BX_GSE", "BY_GSM", "BZ_GSM"], s, e)
    if data is None:
        return None
    bx = np.asarray(data["BX_GSE"], float)
    by = np.asarray(data["BY_GSM"], float)
    bz = np.asarray(data["BZ_GSM"], float)
    bad = (np.abs(bx) > 900) | (np.abs(by) > 900) | (np.abs(bz) > 900)
    bx[bad] = np.nan; by[bad] = np.nan; bz[bad] = np.nan
    bmag = np.sqrt(bx**2 + by**2 + bz**2)
    fin = np.isfinite(bmag) & (bmag > 0)
    if fin.sum() < 60:
        return (np.nan, np.nan, int(fin.sum()))
    cone = np.degrees(np.arccos(np.clip(np.abs(bx[fin]) / bmag[fin], 0, 1)))
    return (float(np.median(cone)),
            float(np.percentile(cone, 75) - np.percentile(cone, 25)),
            int(fin.sum()))


rows = list(done.values())
have = set(done)
n_new = n_fail = 0
for i, r in enumerate(north):
    eid = r["eid"]
    if eid in have:
        continue
    p = os.path.join(CACHE, eid + ".npz")
    if not os.path.exists(p):
        continue
    try:
        t = np.asarray(np.load(p)["tunix"], float)
        t0, t1 = float(np.nanmin(t)), float(np.nanmax(t))
        res = None
        for attempt in range(3):
            try:
                res = fetch_cone(t0, t1)
                break
            except Exception:
                time.sleep(2)
        time.sleep(0.2)  # be gentle on CDAWeb to avoid throttling
        if res is None:
            n_fail += 1
            continue
        med, iqr, n = res
        rows.append(dict(eid=eid, grp=r["grp"], cone28=r["cone"],
                         cone_omni_med=med, cone_iqr=iqr, n_min=n))
        n_new += 1
        if n_new % 25 == 0:
            print(f"  fetched {n_new} new (fail {n_fail})", flush=True)
            with open(CSV, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["eid", "grp", "cone28", "cone_omni_med", "cone_iqr", "n_min"])
                w.writeheader(); w.writerows(rows)
    except Exception as ex:
        n_fail += 1
        continue

with open(CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["eid", "grp", "cone28", "cone_omni_med", "cone_iqr", "n_min"])
    w.writeheader(); w.writerows(rows)
print(f"DONE: {len(rows)} encounters with cone-stability ({n_new} new, {n_fail} failed) -> {CSV}", flush=True)
