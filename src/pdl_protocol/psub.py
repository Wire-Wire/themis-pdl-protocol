"""psub — parallel substrate mapper (ThreadPoolExecutor).

Substrate analysis is local (no network), but the two bottleneck ops both RELEASE THE GIL:
  - np.load on a compressed .npz (zlib decompression), and
  - numpy median / percentile / boolean-masking (C loops).
So threads overlap disk I/O + numpy compute for a solid speedup, without the Windows
multiprocessing spawn/pickle fragility. Same thread model as build_substrate.py.

Usage:
    from psub import pmap, list_files
    results = pmap(profile_one)            # profile_one(npz_dict) -> result or None
    results = pmap(fn, with_name=True)     # fn(filename, npz_dict) -> result or None
"""
import os as _cfg_os, sys as _cfg_sys
from .config import P
import os, glob
from concurrent.futures import ThreadPoolExecutor
import numpy as np

SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
# Powerful machine: use up to 32 (build_substrate's proven-safe ceiling). Load-bound work
# (np.load zlib) scales with threads since decompression releases the GIL.
DEFAULT_WORKERS = min(32, (os.cpu_count() or 8))


def list_files(sub=SUB):
    return sorted(glob.glob(os.path.join(sub, "*.npz")))


def pmap(fn, files=None, workers=DEFAULT_WORKERS, with_name=False):
    """Apply fn to every substrate npz in threads; return list of non-None results.
    fn receives the loaded npz (or (basename, npz) if with_name=True). Result order
    is not guaranteed (aggregation must be order-independent)."""
    files = files if files is not None else list_files()

    def work(f):
        try:
            d = np.load(f, allow_pickle=True)
        except Exception:
            return None
        try:
            return fn(os.path.basename(f), d) if with_name else fn(d)
        except Exception:
            return None

    out = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(work, files):
            if r is not None:
                out.append(r)
    return out
