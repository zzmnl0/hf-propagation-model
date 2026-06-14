"""
plot_ne_background.py – 2-D electron-density background visualisation.

Generates four independent figures saved to output/:
    ne_baseline.png   IRI background only
    ne_tid.png        IRI + TID (MSTID)
    ne_es.png         IRI + Es layer
    ne_bubble.png     IRI + plasma bubble

Colour axis : log10(Ne)  [m^-3], range 9 ~ 12
Contours    : white lines every 0.5 decades, labelled

Usage:
    conda run -n pytorch_cpu python viz/plot_ne_background.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import config as cfg
from models.ionosphere_model import IonosphereModel

# Output directory
_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(_ROOT, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Colormap / range ──────────────────────────────────────────────────────────

CMAP       = 'plasma'
LOG10_VMIN = 9      # lower bound: 10^9 m^-3 (E-layer onset)
LOG10_VMAX = 12     # upper bound: 10^12 m^-3

# Contour levels every 0.5 decades
CONTOUR_LEVELS = np.arange(9.0, 12.01, 0.5)   # 9.0 9.5 10.0 ... 12.0


def _log10Ne(Ne_2d):
    return np.log10(np.maximum(Ne_2d, 1.0))


# ── Core plot ─────────────────────────────────────────────────────────────────

def _plot_ne(Ne_2d, x_km, z_km, title, save_path,
             annot_hlines=(), annot_points=()):
    """
    2-D log10(Ne) colour map with equal-density contours.

    Parameters
    ----------
    annot_hlines : [(z_km, label, color), ...]
    annot_points : [(x_km, z_km, label, color), ...]
    """
    fig, ax = plt.subplots(figsize=(12, 5))

    log_Ne = _log10Ne(Ne_2d)

    # ── Colour mesh ───────────────────────────────────────────────────────────
    pcm = ax.pcolormesh(
        x_km, z_km, log_Ne.T,
        cmap=CMAP, vmin=LOG10_VMIN, vmax=LOG10_VMAX,
        shading='auto'
    )
    cb = plt.colorbar(pcm, ax=ax, extend='both', pad=0.02)
    cb.set_label(r'$\log_{10}\,N_e\ \ [\mathrm{m}^{-3}]$', fontsize=11)
    cb.locator   = mticker.MultipleLocator(1)
    cb.formatter = mticker.FormatStrFormatter('%d')
    cb.update_ticks()

    # ── Contour lines ─────────────────────────────────────────────────────────
    cs = ax.contour(
        x_km, z_km, log_Ne.T,
        levels=CONTOUR_LEVELS,
        colors='white', linewidths=0.6, alpha=0.7
    )
    ax.clabel(cs, inline=True, fontsize=7, fmt='%.1f', colors='white')

    # ── TX / RX ───────────────────────────────────────────────────────────────
    ax.axvline(cfg.TX_POS[0], color='white', ls='--', lw=1.0, label='TX')
    ax.axvline(cfg.RX_POS[0], color='cyan',  ls='--', lw=1.0, label='RX')

    # ── Optional annotations ──────────────────────────────────────────────────
    for z_ann, lbl, col in annot_hlines:
        ax.axhline(z_ann, color=col, ls=':', lw=1.2, alpha=0.85, label=lbl)
    for xp, zp, lbl, col in annot_points:
        ax.plot(xp, zp, marker='+', ms=14, mew=2.2,
                color=col, label=lbl, zorder=5)

    ax.set_xlabel('Horizontal distance  (km)', fontsize=11)
    ax.set_ylabel('Height  (km)', fontsize=11)
    ax.set_xlim(x_km[0], x_km[-1])
    ax.set_ylim(z_km[0], z_km[-1])
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=9, loc='upper right', framealpha=0.6)

    plt.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    print("  Saved:", save_path)
    return fig, ax


# ── Four scenarios ────────────────────────────────────────────────────────────

def plot_baseline():
    print("[1] IRI baseline ...")
    Ne, _ = IonosphereModel().build_Ne_field(cfg.BG_X, cfg.BG_Z)

    mid  = len(cfg.BG_X) // 2
    hmF2 = cfg.BG_Z[Ne[mid].argmax()]
    foF2 = cfg.K_FP * np.sqrt(Ne[mid].max()) / 1e6

    title = (
        'IRI Background'
        '  |  {dt}  {lat}N {lon}E'
        '  |  hmF2={h:.0f} km,  foF2={f:.1f} MHz'
    ).format(
        dt=cfg.IRI_DT.strftime('%Y-%m-%d %H:%M UT'),
        lat=cfg.IRI_LAT, lon=cfg.IRI_LON,
        h=hmF2, f=foF2
    )
    _plot_ne(Ne, cfg.BG_X, cfg.BG_Z, title,
             os.path.join(OUT_DIR, 'ne_baseline.png'),
             annot_hlines=[(hmF2, 'hmF2={:.0f} km'.format(hmF2), 'white')])


def plot_tid():
    print("[2] IRI + TID ...")
    Ne, _ = IonosphereModel(
        tid_params={**cfg.TID, 'enable': True}
    ).build_Ne_field(cfg.BG_X, cfg.BG_Z)

    title = (
        'IRI + TID'
        '  |  lambda={lam:.0f} km,  T={T:.0f} min'
        '  |  dNe/Ne={amp:.0f}%,  I_dip={I:.0f} deg'
    ).format(
        lam=cfg.TID['lambda_h_km'],
        T=cfg.TID['T_s'] / 60.0,
        amp=cfg.TID['amplitude'] * 100,
        I=cfg.TID['I_dip_deg']
    )
    _plot_ne(Ne, cfg.BG_X, cfg.BG_Z, title,
             os.path.join(OUT_DIR, 'ne_tid.png'))


def plot_es():
    print("[3] IRI + Es ...")
    Ne, _ = IonosphereModel(
        es_params={**cfg.ES, 'enable': True}
    ).build_Ne_field(cfg.BG_X, cfg.BG_Z)

    title = (
        'IRI + Es Layer'
        '  |  foEs={foes:.1f} MHz,  h_Es={h:.0f} km,  half-thickness={dh:.0f} m'
    ).format(
        foes=cfg.ES['foEs_MHz'],
        h=cfg.ES['h_Es_km'],
        dh=cfg.ES['delta_h_m']
    )
    _plot_ne(Ne, cfg.BG_X, cfg.BG_Z, title,
             os.path.join(OUT_DIR, 'ne_es.png'),
             annot_hlines=[(cfg.ES['h_Es_km'],
                            'h_Es={:.0f} km'.format(cfg.ES['h_Es_km']),
                            'yellow')])


def plot_bubble():
    print("[4] IRI + plasma bubble ...")
    Ne, _ = IonosphereModel(
        bubble_params={**cfg.BUBBLE, 'enable': True}
    ).build_Ne_field(cfg.BG_X, cfg.BG_Z)

    x0, z0 = cfg.BUBBLE['x0_km'], cfg.BUBBLE['z0_km']
    title = (
        'IRI + Plasma Bubble'
        '  |  delta_max={d:.0%},  centre=({x:.0f}, {z:.0f}) km'
        '  |  Lx={lx:.0f} km,  Lz={lz:.0f} km'
    ).format(
        d=cfg.BUBBLE['delta_max'], x=x0, z=z0,
        lx=cfg.BUBBLE['Lx_km'], lz=cfg.BUBBLE['Lz_km']
    )
    _plot_ne(Ne, cfg.BG_X, cfg.BG_Z, title,
             os.path.join(OUT_DIR, 'ne_bubble.png'),
             annot_points=[(x0, z0,
                            'centre ({:.0f},{:.0f}) km'.format(x0, z0),
                            'white')])


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    plot_baseline()
    plot_tid()
    plot_es()
    plot_bubble()
    print("\nAll figures saved to output/.")
