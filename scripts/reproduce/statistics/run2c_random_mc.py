"""RUN 2c — realistic model-error band (per-encounter random standoff MC).
Loads substrate ONCE (threaded via psub); the NMC Monte-Carlo iterations run in
parallel threads (numpy median/masking release the GIL)."""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys
import numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
from radial_models import jelinek_bs_r
from psub import pmap

OUT = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run2_profile_cube")
os.makedirs(OUT, exist_ok=True)
BG = (0.6, 1.0)
EDGES = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00]
MIN_SAMP_BIN = 3; MIN_BG = 5; MIN_NENC = 30; NB = len(EDGES) - 1
NMC = 200; SIG_MP = 0.8; SIG_BS = 0.10
WORKERS = min(16, (os.cpu_count() or 8))


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


def one_iter(data, rng):
    perDn = [[] for _ in range(NB)]; perEB = [[] for _ in range(NB)]
    for (r, flare, r_bs0, mp0, n, b) in data:
        r_mp = (mp0 + rng.normal(0, SIG_MP)) * flare
        r_bs = r_bs0 * max(0.5, rng.normal(1.0, SIG_BS))
        d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
        s = np.where(denom > 0, d_mp / denom, np.nan)
        valid = np.isfinite(s) & (d_mp > 0) & (d_bs > 0)
        sv = s[valid]; nv = n[valid]; bv = b[valid]
        if len(sv) < MIN_BG + MIN_SAMP_BIN:
            continue
        bgm = (sv >= BG[0]) & (sv <= BG[1])
        if bgm.sum() < MIN_BG:
            continue
        n_bg = np.median(nv[bgm]); b_bg = np.median(bv[bgm])
        if not (n_bg > 0 and b_bg > 0):
            continue
        for i in range(NB):
            lo, hi = EDGES[i], EDGES[i + 1]
            m = (sv >= lo) & (sv <= hi) if i == NB - 1 else (sv >= lo) & (sv < hi)
            if m.sum() < MIN_SAMP_BIN:
                continue
            perDn[i].append(np.median(nv[m]) / n_bg); perEB[i].append(np.median(bv[m]) / b_bg)
    Dn = [np.median(perDn[i]) if len(perDn[i]) >= MIN_NENC else np.nan for i in range(NB)]
    EB = [np.median(perEB[i]) if len(perEB[i]) >= MIN_NENC else np.nan for i in range(NB)]
    return Dn, EB


def run_chunk(seed, niter, data):
    rng = np.random.default_rng(seed)
    return [one_iter(data, rng) for _ in range(niter)]


def main():
    data = [e for e in pmap(load_enc) if e is not None]   # threaded load (np.load releases the GIL -> helps)
    print(f"loaded {len(data)} encounters; NMC={NMC} (SIG_MP={SIG_MP} R_E, SIG_BS={SIG_BS})", flush=True)
    # MC runs SERIALLY: the inner per-encounter loop is Python/GIL-bound (many small numpy ops),
    # so threading it only adds GIL contention (measured ~2.5x SLOWER). The load above is threaded.
    rng = np.random.default_rng(7)
    allit = []
    for k in range(NMC):
        allit.append(one_iter(data, rng))
        if (k + 1) % 50 == 0:
            print(f"  iter {k+1}/{NMC}", flush=True)
    Dn = np.array([it[0] for it in allit]); EB = np.array([it[1] for it in allit])

    def band(A, label):
        lines = [f"=== {label} — realistic per-encounter model-error band (random mp0/BS) ===",
                 f"{'shell':13s} {'median':>8s} {'p2.5':>7s} {'p97.5':>7s}"]
        for i in range(NB):
            col = A[:, i]; col = col[np.isfinite(col)]
            if len(col) < len(allit) * 0.5:
                lines.append(f"[{EDGES[i]:.2f},{EDGES[i+1]:.2f}]   (insufficient)"); continue
            lines.append(f"[{EDGES[i]:.2f},{EDGES[i+1]:.2f}]  {np.median(col):8.3f} {np.percentile(col,2.5):7.3f} {np.percentile(col,97.5):7.3f}")
        return "\n".join(lines)

    txt = band(Dn, "Dn") + "\n\n" + band(EB, "EB")
    print("\n" + txt, flush=True)
    with open(os.path.join(OUT, "RUN2c_model_error_band.txt"), "w") as f:
        f.write(f"RUN 2c — realistic per-encounter random standoff MC (NMC={NMC}, SIG_MP={SIG_MP} R_E, SIG_BS={SIG_BS}; threaded)\n\n" + txt + "\n")
    print("\nsaved -> RUN2c_model_error_band.txt", flush=True)


if __name__ == "__main__":
    main()
