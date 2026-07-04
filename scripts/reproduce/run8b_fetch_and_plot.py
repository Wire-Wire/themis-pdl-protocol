"""RUN 8b — per-event atlas: fetch ESA ion energy spectra + plot diagnostic panels.

The substrate already holds the moments (n, |B|, beta, T, v, position); we re-fetch ONLY the
ESA reduced-ion energy spectrogram (th{p}_peir_en_eflux, 32 energies) for each candidate window.
A real magnetosheath PDL must keep a SHOCKED-SHEATH ion spectrum (broad, ~100s eV) as s->0;
magnetosphere/boundary contamination would show a different spectrum.

Panels per event: (1) ion energy spectrogram, (2) n, (3) |B|, (4) beta, (5) T, (6) radial s(t)
with the PDL shell shaded. PNG per event -> runs/run8_atlas/panels/.

Usage: python run8b_fetch_and_plot.py [start] [count]   (smoke test: defaults 0 10)
"""
import os as _cfg_os, sys as _cfg_sys
_cfg_sys.path.insert(0, _cfg_os.path.dirname(_cfg_os.path.abspath(__file__)))
from config import P
import os, sys, csv, datetime as dt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, P(r"H:\0mssl\review\repair\option3"))
sys.path.insert(0, P(r"H:\0mssl\0mssl519\src"))
from radial_models import jelinek_bs_r
from cdasws import CdasWs

SUB = P(r"H:\0mssl\review\01_CURRENT__rebuild\substrate")
ATLAS = P(r"H:\0mssl\review\01_CURRENT__rebuild\runs\run8_atlas")
PANELS = os.path.join(ATLAS, "panels")
os.makedirs(PANELS, exist_ok=True)
KB = 1.602e-4
_tls = CdasWs()


def e2u(a):
    if hasattr(a, 'dtype') and np.issubdtype(a.dtype, np.datetime64):
        return (a.astype('datetime64[ns]') - np.datetime64('1970-01-01T00:00:00', 'ns')).astype(np.float64) / 1e9
    return np.asarray(a, dtype=np.float64)


def s_jbs(x, y, z, mp0, alpha, dp):
    r = np.sqrt(x * x + y * y + z * z)
    ct = np.clip(np.divide(x, r, out=np.ones_like(r), where=r > 0), -1.0, 1.0)
    r_mp = mp0 * (2.0 / (1.0 + ct)) ** alpha; r_bs = jelinek_bs_r(ct, dp)
    d_mp = r - r_mp; d_bs = r_bs - r; denom = d_mp + d_bs
    return np.where(denom > 0, d_mp / denom, np.nan)


def fetch_spec(probe, t0iso, t1iso):
    pl = probe[-1]; var = f"th{pl}_peir_en_eflux"
    for attempt in range(3):
        try:
            st, d = _tls.get_data(f"TH{pl.upper()}_L2_ESA", [var], t0iso, t1iso)
            if d is None or var not in getattr(d, 'data_vars', {}):
                return None
            spec = np.asarray(d[var].values, float)            # (Nt, 32)
            ec = [c for c in d.coords if 'epoch' in c.lower()]
            et = e2u(d.coords[ec[0]].values) if ec else None
            yc = [c for c in d.coords if 'yaxis' in c.lower()]
            energy = np.asarray(d.coords[yc[0]].values, float) if yc else None
            if energy is not None and energy.ndim == 2:
                energy = np.nanmedian(energy, axis=0)
            return et, energy, spec
        except Exception:
            pass
    return None


def panel(meta):
    eid = meta['eid']; probe = meta['probe']
    d = np.load(os.path.join(SUB, eid + ".npz"), allow_pickle=True)
    t = d['t'].astype(float); n = d['n'].astype(float); b = d['bmag'].astype(float)
    beta = d['beta'].astype(float); pth = d['p_th'].astype(float); v = d['vmag'].astype(float)
    x = d['x_re'].astype(float); y = d['y_re'].astype(float); z = d['z_re'].astype(float)
    mp0 = float(d['mp0']); alpha = float(d['alpha']); dp = float(d['dp'])
    s = s_jbs(x, y, z, mp0, alpha, dp)
    T = np.divide(pth, n * KB, out=np.full_like(pth, np.nan), where=(n > 0))
    t0 = float(t.min())
    hh = (t - t0) / 3600.0
    sp = fetch_spec(probe, meta['t_start'], meta['t_end'])

    fig, ax = plt.subplots(6, 1, figsize=(11, 13), sharex=True)
    cat = meta.get('category', 'CONF'); rk = meta.get('rank', '0')
    extra = f" n_near={meta['n_near']} beta_near={meta['beta_near']} T_near={meta['T_near']}" if 'T_near' in meta else ""
    fig.suptitle(f"{cat} #{rk}  {eid}   Dn_mem={meta['Dn_mem']} EB_mem={meta['EB_mem']} "
                 f"n={meta['n_mem']}  Bz={meta['bz']} cone={meta['cone']} SZA={meta['sza']} Dp={meta['dp']}{extra}",
                 fontsize=10)
    # (1) spectrogram — coords (time x, energy y) MUST be finite for pcolormesh; mask only the colour array
    spec_ok = False
    if sp is not None and sp[0] is not None and sp[1] is not None:
        et, energy, spec = sp
        try:
            emask = np.isfinite(energy) & (energy > 0)
            tmask = np.isfinite(et)
            energy2 = energy[emask]; sh = (et[tmask] - t0) / 3600.0
            spec2 = spec[np.ix_(tmask, emask)]
            order = np.argsort(energy2); energy2 = energy2[order]; spec2 = spec2[:, order]
            with np.errstate(invalid='ignore', divide='ignore'):
                logC = np.log10(np.where(spec2 > 0, spec2, np.nan)).T   # (Ne, Nt)
            fin = logC[np.isfinite(logC)]
            if fin.size > 10 and energy2.size > 2 and sh.size > 2:
                vmin, vmax = np.percentile(fin, [5, 99])
                pm = ax[0].pcolormesh(sh, energy2, np.ma.masked_invalid(logC), shading='auto', cmap='jet', vmin=vmin, vmax=vmax)
                ax[0].set_yscale('log'); ax[0].set_ylabel('ion E [eV]')
                fig.colorbar(pm, ax=ax[0], pad=0.01, label='log eflux')
                ax[0].set_title('ESA reduced-ion energy spectrogram', fontsize=9)
                spec_ok = True
            else:
                ax[0].text(0.5, 0.5, 'spectra all-invalid', transform=ax[0].transAxes, ha='center', color='red')
        except Exception as ex:
            ax[0].text(0.5, 0.5, f'spec plot failed: {ex}', transform=ax[0].transAxes, ha='center', fontsize=7)
    else:
        ax[0].text(0.5, 0.5, 'NO ESA SPECTRA RETURNED', transform=ax[0].transAxes, ha='center', color='red')
    # (2-5) moments
    ax[1].plot(hh, n, '.', ms=2); ax[1].set_yscale('log'); ax[1].set_ylabel('n [cm^-3]')
    bgm = (s >= 0.6) & (s <= 1.0) & np.isfinite(n)
    if bgm.any():
        ax[1].axhline(np.median(n[bgm]), color='g', ls='--', lw=1, label='bg median'); ax[1].legend(fontsize=7)
    ax[2].plot(hh, b, '.', ms=2, color='purple'); ax[2].set_ylabel('|B| [nT]')
    ax[3].plot(hh, beta, '.', ms=2, color='brown'); ax[3].set_yscale('log'); ax[3].set_ylabel('beta'); ax[3].axhline(1, color='k', ls=':', lw=1)
    ax[4].plot(hh, T, '.', ms=2, color='orange'); ax[4].set_yscale('log'); ax[4].set_ylabel('T [eV]')
    # (6) radial s with PDL shell
    ax[5].plot(hh, s, '.', ms=2, color='navy'); ax[5].set_ylabel('s (radial)'); ax[5].set_ylim(-0.05, 1.05)
    ax[5].axhspan(0.05, 0.20, color='gold', alpha=0.3, label='PDL shell'); ax[5].axhline(0, color='r', ls='--', lw=1, label='model MP')
    ax[5].legend(fontsize=7); ax[5].set_xlabel(f'hours from {meta["t_start"]}')
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    outp = os.path.join(PANELS, f"{cat}_{int(rk):02d}_{eid}.png")
    fig.savefig(outp, dpi=90); plt.close(fig)
    return spec_ok, outp


def main():
    start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    csvpath = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ATLAS, "candidates.csv")
    rows = list(csv.DictReader(open(csvpath)))
    todo = rows[start:start + count]
    print(f"plotting {len(todo)} events (rows {start}..{start+len(todo)-1})", flush=True)
    ok = 0; nospec = []
    for r in todo:
        try:
            got, outp = panel(r)
            ok += 1
            tag = "spec OK" if got else "NO SPEC"
            if not got:
                nospec.append(r['eid'])
            print(f"  {r.get('category','CONF'):8s} {r['eid']:22s} -> {tag}  {os.path.basename(outp)}", flush=True)
        except Exception as ex:
            print(f"  {r['eid']:22s} FAILED: {ex!r}", flush=True)
    print(f"\ndone: {ok}/{len(todo)} plotted; no-spectra: {nospec}", flush=True)
    print("panels ->", PANELS, flush=True)


if __name__ == "__main__":
    main()
