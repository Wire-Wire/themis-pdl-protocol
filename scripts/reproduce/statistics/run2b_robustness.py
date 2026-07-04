"""RUN 2b — robustness: boundary-standoff SYSTEMATIC sensitivity + bootstrap.
Loads the substrate ONCE (threaded via psub), then computes all settings in memory
(was 5 full reloads; now 1 threaded load + in-memory passes)."""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys
import numpy as np
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run2_profile_cube")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0)
EDGES = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00]
MIN_SAMP_BIN = 3; MIN_BG = 5; MIN_NENC = 30; NB = len(EDGES) - 1
RNG = np.random.default_rng(42)


def load_enc(d):
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    if not (np.isfinite(mp0) and np.isfinite(alpha) and np.isfinite(dp) and dp > 0):
        return None
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    n = d['n'].astype(float); b = d['bmag'].astype(float)
    base = np.isfinite(n) & (n > 0) & np.isfinite(b) & (b > 0) & np.isfinite(x)
    if base.sum() < MIN_BG + MIN_SAMP_BIN:
        return None
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    flare = (2.0 / (1.0 + ct)) ** alpha
    r_bs0 = jelinek_bs_r(ct, dp)
    return (r[base], flare[base], r_bs0[base], mp0, n[base], b[base])


def prof_mem(enc, off, sc):
    r, flare, r_bs0, mp0, n, b = enc
    r_mp = (mp0 + off) * flare; r_bs = r_bs0 * sc
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    s = np.where(denom > 0, d_mp / denom, np.nan)
    valid = np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
    sv = s[valid]; nv = n[valid]; bv = b[valid]
    if len(sv) < MIN_BG + MIN_SAMP_BIN:
        return None
    bgm = (sv >= BG[0]) & (sv <= BG[1])
    if bgm.sum() < MIN_BG:
        return None
    n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm])
    if not (n_bg > 0 and b_bg > 0):
        return None
    out = {}
    for i in range(NB):
        lo, hi = EDGES[i], EDGES[i + 1]
        m = (sv >= lo) & (sv <= hi) if i == NB - 1 else (sv >= lo) & (sv < hi)
        if m.sum() < MIN_SAMP_BIN:
            continue
        out[i] = (float(np.median(nv[m]) / n_bg), float(np.median(bv[m]) / b_bg))
    return out


def aggregate(data, off, sc):
    perDn = {i: [] for i in range(NB)}; perEB = {i: [] for i in range(NB)}
    for enc in data:
        p = prof_mem(enc, off, sc)
        if not p:
            continue
        for i, (dn, eb) in p.items():
            perDn[i].append(dn); perEB[i].append(eb)
    return perDn, perEB


def main():
    data = [e for e in pmap(load_enc) if e is not None]
    print(f"loaded {len(data)} encounters (threaded, once)", flush=True)
    perDn, perEB = aggregate(data, 0.0, 1.0)
    settings = [('mp0-0.8', -0.8, 1.0), ('mp0+0.8', 0.8, 1.0), ('bs*0.9', 0.0, 0.9), ('bs*1.1', 0.0, 1.1)]
    sens = {name: aggregate(data, off, sc) for name, off, sc in settings}

    def block(per, idx, label):
        out = [f"=== {label} robustness (sheath_geom) ===",
               f"{'shell':13s} {'N':>5s} {'nom':>7s} {'boot95':>16s} " + " ".join(f'{n:>8s}' for n, _, _ in settings)]
        for i in range(NB):
            v = per[i]
            if len(v) < MIN_NENC:
                out.append(f"[{EDGES[i]:.2f},{EDGES[i+1]:.2f}]  {len(v):5d}   (N<{MIN_NENC})"); continue
            arr = np.array(v); med = np.median(arr)
            boot = np.array([np.median(RNG.choice(arr, len(arr), replace=True)) for _ in range(1000)])
            lo, hi = np.percentile(boot, [2.5, 97.5])
            sline = " ".join(
                f'{(np.median(sens[n][idx][i]) if len(sens[n][idx][i]) >= MIN_NENC else float("nan")):8.3f}'
                for n, _, _ in settings)
            out.append(f"[{EDGES[i]:.2f},{EDGES[i+1]:.2f}]  {len(v):5d} {med:7.3f} [{lo:6.3f},{hi:6.3f}] {sline}")
        return "\n".join(out)

    txt = block(perDn, 0, "Dn") + "\n\n" + block(perEB, 1, "EB")
    print("\n" + txt, flush=True)
    with open(os.path.join(OUT, "RUN2b_robustness.txt"), "w") as f:
        f.write("RUN 2b — standoff-error + bootstrap robustness (load-once, threaded)\n")
        f.write("mp0 offset in R_E; bs scale multiplies the Jelinek standoff; boot95 = encounter-bootstrap 95% CI.\n\n")
        f.write(txt + "\n")
    print("\nsaved -> RUN2b_robustness.txt", flush=True)


if __name__ == "__main__":
    main()
