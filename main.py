"""
main.py - Entry point for the HF hybrid propagation model.

Usage:
    conda run -n pytorch_cpu python main.py

Scenarios (each saves ray_fan_*.png and pd_*.png to output/):
    run_baseline()    IRI only          -> ray_fan_baseline, pd_baseline
    run_with_tid()    IRI + TID         -> ray_fan_tid,      pd_tid
    run_with_es()     IRI + Es          -> ray_fan_es,       pd_es
    run_with_bubble() IRI + bubble      -> ray_fan_bubble,   pd_bubble
    run_full()        TID + Es + bubble -> ray_fan_full,     pd_full
                      (full uses P2P mode-path overlay instead of fan)

Test link: TX @ 30N 120E,  RX @ 1169 km,  f = 10 MHz.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use('Agg')          # non-interactive backend; avoids GUI windows
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import config as cfg
from viz.plot_utils import print_mode_table, plot_ray_fan, plot_pd_spectrum

_OUT = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(_OUT, exist_ok=True)


# ── Shared helpers ────────────────────────────────────────────────────────────

def _save(fig, name):
    """Save figure to output/ and close it."""
    fig.savefig(os.path.join(_OUT, name), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: output/{}".format(name))


def _make_iono_params(enable_tid=False, enable_es=False, enable_bubble=False):
    """Assemble iono_params dict from config defaults with selected flags."""
    return {
        'iri_params':    {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':    {**cfg.TID,    'enable': enable_tid},
        'es_params':     {**cfg.ES,     'enable': enable_es},
        'bubble_params': {**cfg.BUBBLE, 'enable': enable_bubble},
    }


def _make_radar_params():
    return {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W,
            'Gt': cfg.GT, 'Gr': cfg.GR}


def _print_main(main):
    if main:
        pr_str = ("{:.1f} dBW".format(main['Pr_dBW'])
                  if 'Pr_dBW' in main else "n/a")
        print("  Main mode: {}  tau = {:.3f} ms  Pr = {}".format(
            main['label'], main['tau_ms'], pr_str))
    else:
        print("  Main mode: None")


def _add_free_space_power(modes):
    """
    Fill Pr_W, Pr_dBW, delta_tau_ms into mode dicts returned by
    find_all_rays_p2p (which does not compute power).
    Used in baseline and TID scenarios where HybridModel is not called.
    """
    from utils import free_space_loss_dB, to_dBW
    Pt = cfg.PT_W; Gt = cfg.GT; Gr = cfg.GR
    for m in modes:
        L_dB       = free_space_loss_dB(m['group_path_km'], cfg.FREQ_MHZ)
        m['Pr_W']          = Pt * Gt * Gr * 10.0 ** (-L_dB / 10.0)
        m['Pr_dBW']        = to_dBW(m['Pr_W'])
        m['delta_tau_ms']  = 0.0


def _plot_mode_paths(modes, Ne_2d, x_km, z_km, title='P2P Mode Paths'):
    """
    Plot converged P2P variational paths (one per mode) over the Ne field.
    Used in run_full() where the full fan is replaced by the P2P solutions.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.pcolormesh(x_km, z_km, Ne_2d.T / 1e10,
                  cmap='Blues', shading='auto', alpha=0.35)
    colors = cm.tab10(np.linspace(0, 1, max(len(modes), 1)))
    for m, c in zip(modes, colors):
        pts = m['points']          # (n_ctrl+2, 2) control-point array
        ax.plot(pts[:, 0], pts[:, 1], '-', color=c, lw=1.8,
                label='{} ({:.2f} ms)'.format(m['label'], m['tau_ms']))
    ax.axvline(cfg.TX_POS[0], color='g', lw=1.2, ls='--', label='TX')
    ax.axvline(cfg.RX_POS[0], color='b', lw=1.2, ls='--', label='RX')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0, 500])
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    return fig, ax


# ── Environment check ─────────────────────────────────────────────────────────

def check_env() -> bool:
    """Verify that all required packages are installed and print versions."""
    ok = True
    print("=" * 55)
    print("  Environment check")
    print("=" * 55)
    for name in ('numpy', 'scipy', 'matplotlib', 'iri2016'):
        try:
            mod = __import__(name)
            ver = getattr(mod, '__version__', 'unknown')
            print("  OK      {:<14} {}".format(name, ver))
        except ImportError:
            print("  MISSING {:<14}  <- pip install {}".format(name, name))
            ok = False
    print("=" * 55)
    print("  Test link : TX @ {}N {}E".format(cfg.TX_LAT, cfg.TX_LON))
    print("  RX range  : {} km    freq = {} MHz".format(cfg.RX_RANGE, cfg.FREQ_MHZ))
    print()
    return ok


# ── Scenario 0: IRI baseline ──────────────────────────────────────────────────

def run_baseline():
    """IRI background only. Saves ray_fan_baseline.png and pd_baseline.png."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.point_to_point import find_all_rays_p2p
    from models.hybrid_model import build_pd_spectrum, identify_main_mode

    print("=" * 55)
    print("  Scenario 0: Baseline (IRI only)")
    print("=" * 55)

    iono = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    print("  Ne field: shape {}  max = {:.2e} m^-3".format(
        Ne_2d.shape, float(Ne_2d.max())))

    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

    # Ray fan
    rays = shoot_rays_fan(cfg.TX_POS, n_model)
    print("  Ray fan : {} rays".format(len(rays)))
    fig, _ = plot_ray_fan(rays, Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI baseline')
    _save(fig, 'ray_fan_baseline.png')

    # P2P modes + free-space power
    modes = find_all_rays_p2p(cfg.TX_POS, cfg.RX_POS, n_model, cfg.FREQ_MHZ)
    _add_free_space_power(modes)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    # P-D spectrum
    tau_ax, pd = build_pd_spectrum(modes)
    main, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D - IRI baseline')
    _save(fig2, 'pd_baseline.png')


# ── Scenario 1: IRI + TID ────────────────────────────────────────────────────

def run_with_tid():
    """IRI + MSTID. Saves ray_fan_tid.png and pd_tid.png."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.point_to_point import find_all_rays_p2p
    from models.hybrid_model import build_pd_spectrum, identify_main_mode

    print("=" * 55)
    print("  Scenario 1: IRI + TID")
    print("=" * 55)

    iono = IonosphereModel(tid_params={**cfg.TID, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

    # Ray fan
    rays = shoot_rays_fan(cfg.TX_POS, n_model)
    print("  Ray fan : {} rays".format(len(rays)))
    fig, _ = plot_ray_fan(rays, Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + TID')
    _save(fig, 'ray_fan_tid.png')

    # P2P modes + free-space power
    modes = find_all_rays_p2p(cfg.TX_POS, cfg.RX_POS, n_model, cfg.FREQ_MHZ)
    _add_free_space_power(modes)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    # P-D spectrum
    tau_ax, pd = build_pd_spectrum(modes)
    main, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D - IRI + TID')
    _save(fig2, 'pd_tid.png')


# ── Scenario 2: IRI + Es ─────────────────────────────────────────────────────

def run_with_es():
    """
    IRI + sporadic-E layer. Saves ray_fan_es.png and pd_es.png.
    Ray fan uses IRI+Es Ne field.  Power computation via HybridPropagationModel.
    """
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.hybrid_model import HybridPropagationModel

    print("=" * 55)
    print("  Scenario 2: IRI + Es")
    print("=" * 55)
    print("  foEs = {} MHz  f = {} MHz  foEs/f = {:.2f}".format(
        cfg.ES['foEs_MHz'], cfg.FREQ_MHZ,
        cfg.ES['foEs_MHz'] / cfg.FREQ_MHZ))

    # Build Ne field for ray fan visualisation
    iono = IonosphereModel(es_params={**cfg.ES, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

    rays = shoot_rays_fan(cfg.TX_POS, n_model)
    print("  Ray fan : {} rays".format(len(rays)))
    fig, _ = plot_ray_fan(rays, Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + Es')
    _save(fig, 'ray_fan_es.png')

    # Full power pipeline via HybridPropagationModel (Decision 1: Option A)
    model = HybridPropagationModel(
        _make_iono_params(enable_es=True), _make_radar_params())
    modes, tau_ax, pd, main = model.compute(cfg.TX_POS, cfg.RX_POS)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    _print_main(main)

    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D - IRI + Es')
    _save(fig2, 'pd_es.png')


# ── Scenario 3: IRI + plasma bubble ──────────────────────────────────────────

def run_with_bubble():
    """
    IRI + plasma bubble (PE/SSF). Saves ray_fan_bubble.png and pd_bubble.png.
    Ray fan uses bubble-modified Ne field (depletion visible in background).
    """
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.hybrid_model import HybridPropagationModel

    print("=" * 55)
    print("  Scenario 3: IRI + Plasma Bubble")
    print("=" * 55)
    print("  Bubble: x0={} km  z0={} km  Lx={} km  Lz={} km  dmax={}".format(
        cfg.BUBBLE['x0_km'], cfg.BUBBLE['z0_km'],
        cfg.BUBBLE['Lx_km'], cfg.BUBBLE['Lz_km'], cfg.BUBBLE['delta_max']))

    # Build Ne field for ray fan visualisation (bubble-modified background)
    iono = IonosphereModel(bubble_params={**cfg.BUBBLE, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

    rays = shoot_rays_fan(cfg.TX_POS, n_model)
    print("  Ray fan : {} rays".format(len(rays)))
    fig, _ = plot_ray_fan(rays, Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + Plasma Bubble')
    _save(fig, 'ray_fan_bubble.png')

    # Full power pipeline via HybridPropagationModel (PE/SSF for bubble)
    model = HybridPropagationModel(
        _make_iono_params(enable_bubble=True), _make_radar_params())
    modes, tau_ax, pd, main = model.compute(cfg.TX_POS, cfg.RX_POS)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    _print_main(main)

    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - IRI + Plasma Bubble')
    _save(fig2, 'pd_bubble.png')


# ── Scenario 4: Full model ────────────────────────────────────────────────────

def run_full():
    """
    IRI + TID + Es + bubble. Saves ray_fan_full.png and pd_full.png.
    ray_fan_full shows P2P mode paths (Decision 2: Option C) so each
    labelled mode is individually visible against the composite Ne field.
    """
    from models.hybrid_model import HybridPropagationModel

    print("=" * 55)
    print("  Scenario 4: Full model (TID + Es + Bubble)")
    print("=" * 55)

    model = HybridPropagationModel(
        _make_iono_params(enable_tid=True, enable_es=True, enable_bubble=True),
        _make_radar_params())
    modes, tau_ax, pd, main = model.compute(cfg.TX_POS, cfg.RX_POS)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    _print_main(main)

    # Mode-path overlay: rebuild Ne from the model's iono (all effects active)
    Ne_2d, _ = model.iono.build_Ne_field(model.x_array, model.z_array)
    fig, _ = _plot_mode_paths(
        modes, Ne_2d, model.x_array, model.z_array,
        title='P2P Mode Paths - Full model (TID + Es + Bubble)')
    _save(fig, 'ray_fan_full.png')

    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - Full model (TID + Es + Bubble)')
    _save(fig2, 'pd_full.png')


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    if not check_env():
        sys.exit(1)

    run_baseline()
    run_with_tid()
    run_with_es()
    run_with_bubble()
    run_full()
