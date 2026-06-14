"""
main.py – Entry point for the HF hybrid propagation model.

Usage:
    conda run -n pytorch_cpu python main.py

Scenarios (uncomment to run):
    check_env()          – verify all packages are importable
    run_baseline()       – IRI only, ray fan, print mode table
    run_with_tid()       – IRI + TID
    run_with_es()        – IRI + Es
    run_with_bubble()    – IRI + plasma bubble
    run_full()           – all three irregularities combined

Each scenario uses the test link defined in config.py:
    TX @ 30N 120E,  RX @ 1169 km,  f = 10 MHz.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))   # make sure local imports work

import numpy as np
import config as cfg
from viz.plot_utils import print_mode_table, plot_ne_field, plot_ray_fan, plot_pd_spectrum

_OUT = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(_OUT, exist_ok=True)


# ── Environment check ─────────────────────────────────────────────────────────

def check_env() -> bool:
    """Verify that all required packages are installed and print versions."""
    ok = True
    packages = {
        'numpy':      'np.__version__',
        'scipy':      'scipy.__version__',
        'matplotlib': 'matplotlib.__version__',
        'iri2016':    'iri2016.__version__',
    }
    print("=" * 55)
    print("  Environment check")
    print("=" * 55)
    for name, ver_attr in packages.items():
        try:
            mod = __import__(name)
            ver = getattr(mod, '__version__', 'unknown')
            print(f"  {'OK':>4}  {name:<14} {ver}")
        except ImportError:
            print(f"  {'MISSING':>4}  {name:<14}  <- pip install {name}")
            ok = False
    print("=" * 55)
    if ok:
        print("  All packages present.\n")
    else:
        print("  Please install missing packages before continuing.\n")

    # Show test-link summary
    print(f"  Test link:  TX @ {cfg.TX_LAT}N {cfg.TX_LON}E")
    print(f"              RX @ {cfg.RX_RANGE} km,  f = {cfg.FREQ_MHZ} MHz")
    print(f"  BG grid  :  x [{cfg.BG_X_MIN}, {cfg.BG_X_MAX}] km  dx={cfg.BG_DX} km  "
          f"({len(cfg.BG_X)} pts)")
    print(f"              z [{cfg.BG_Z_MIN}, {cfg.BG_Z_MAX}] km  dz={cfg.BG_DZ} km  "
          f"({len(cfg.BG_Z)} pts)")
    lam_m = cfg.C_MS / (cfg.FREQ_MHZ * 1e6)
    print(f"  lambda = {lam_m:.1f} m   ->   PE dz = {cfg.PE['dz_m']:.1f} m  "
          f"(lambda/4),  PE dx = {cfg.PE['dx_km']*1000:.0f} m")
    print()
    return ok


# ── Scenario helpers ──────────────────────────────────────────────────────────

def _make_iono_params(enable_tid=False, enable_es=False, enable_bubble=False):
    """Assemble iono_params dict from config defaults with selected flags."""
    tid    = {**cfg.TID,    'enable': enable_tid}
    es     = {**cfg.ES,     'enable': enable_es}
    bubble = {**cfg.BUBBLE, 'enable': enable_bubble}
    return {
        'iri_params':    {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':    tid,
        'es_params':     es,
        'bubble_params': bubble,
    }

def _make_radar_params():
    return {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W, 'Gt': cfg.GT, 'Gr': cfg.GR}


# ── Runnable scenarios ────────────────────────────────────────────────────────

def run_baseline():
    """
    Scenario 0 – IRI background only, no irregularities.
    Expected results:
        - 2 dominant modes: 1F_low and 1F_high
        - Group delays ~ 5-8 ms (1169 km at 10 MHz)
        - Group path error vs. equivalent-path theorem < 1 %
    """
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.point_to_point import find_all_rays_p2p

    print("=" * 55)
    print("  Scenario: Baseline (IRI only)")
    print("=" * 55)

    # 1. Build density field
    iono = IonosphereModel()
    Ne_2d, n_2d = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    print(f"  Ne field built:  shape {Ne_2d.shape},  "
          f"max = {Ne_2d.max():.2e} m^-3")

    # 2. Build refractive-index model
    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

    # 3. Ray fan (free-end, visualisation)
    rays = shoot_rays_fan(cfg.TX_POS, n_model)
    print(f"  Ray fan: {len(rays)} rays,  "
          f"tau range [{min(r['tau_ms'] for r in rays):.2f}, "
          f"{max(r['tau_ms'] for r in rays):.2f}] ms")
    fig, _ = plot_ray_fan(rays, Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan – IRI baseline')
    fig.savefig(os.path.join(_OUT, 'ray_fan_baseline.png'), dpi=150, bbox_inches='tight')

    # 4. Point-to-point: find all modes
    modes = find_all_rays_p2p(cfg.TX_POS, cfg.RX_POS, n_model, cfg.FREQ_MHZ)
    print(f"\n  P2P modes found: {len(modes)}")
    print_mode_table(modes)

    # 5. P-D spectrum
    from models.hybrid_model import build_pd_spectrum, identify_main_mode
    tau_ax, pd = build_pd_spectrum(modes)
    main, _ = identify_main_mode(tau_ax, pd, modes)
    print(f"\n  Main mode: {main.get('label', '?') if main else 'None'}  "
          f"tau = {main['tau_ms']:.3f} ms" if main else "  Main mode: None")
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D – IRI baseline')
    fig2.savefig(os.path.join(_OUT, 'pd_baseline.png'), dpi=150, bbox_inches='tight')
    print("  Saved: output/ray_fan_baseline.png, output/pd_baseline.png\n")


def run_with_tid():
    """Scenario 1 – IRI + TID (MSTID, lambda=300 km, T=40 min, 10% amplitude)."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.point_to_point import find_all_rays_p2p

    print("=" * 55)
    print("  Scenario: IRI + TID")
    print("=" * 55)
    iono = IonosphereModel(tid_params={**cfg.TID, 'enable': True})
    Ne_2d, n_2d = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)
    rays = shoot_rays_fan(cfg.TX_POS, n_model)
    fig, _ = plot_ray_fan(rays, Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan – IRI + TID')
    fig.savefig(os.path.join(_OUT, 'ray_fan_tid.png'), dpi=150, bbox_inches='tight')
    modes = find_all_rays_p2p(cfg.TX_POS, cfg.RX_POS, n_model, cfg.FREQ_MHZ)
    print(f"  Modes found: {len(modes)}")
    print_mode_table(modes)
    print("  Saved: output/ray_fan_tid.png\n")


def run_with_es():
    """Scenario 2 – IRI + Es  (foEs = 5 MHz, f = 10 MHz  ->  foEs/f = 0.5 -> reflect)."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer import RefractiveIndex, shoot_rays_fan
    from models.point_to_point import find_all_rays_p2p, extract_es_params
    from models.es_model import EsLayerModel

    print("=" * 55)
    print("  Scenario: IRI + Es")
    print("=" * 55)
    iono = IonosphereModel(es_params={**cfg.ES, 'enable': True})
    Ne_2d, n_2d = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n_model = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)
    es_model = EsLayerModel()

    mode, alpha = es_model.classify(cfg.FREQ_MHZ)
    print(f"  foEs/f = {cfg.ES['foEs_MHz']}/{cfg.FREQ_MHZ} = "
          f"{cfg.ES['foEs_MHz']/cfg.FREQ_MHZ:.2f}  ->  {mode}  (alpha={alpha:.2f})")

    modes = find_all_rays_p2p(cfg.TX_POS, cfg.RX_POS, n_model, cfg.FREQ_MHZ)
    print(f"  Modes found: {len(modes)}")
    print_mode_table(modes)
    print()


def run_with_bubble():
    """Scenario 3 – IRI + plasma bubble (PE/SSF path)."""
    from models.hybrid_model import HybridPropagationModel, build_pd_spectrum

    print("=" * 55)
    print("  Scenario: IRI + Plasma Bubble")
    print("=" * 55)
    iono_p = _make_iono_params(enable_bubble=True)
    model  = HybridPropagationModel(iono_p, _make_radar_params())
    modes, tau_ax, pd, main = model.compute(cfg.TX_POS, cfg.RX_POS)
    print(f"  Modes found: {len(modes)}")
    print_mode_table(modes)
    print(f"  Main mode: {main.get('label','?') if main else 'None'}\n")


def run_full():
    """Scenario 4 – IRI + TID + Es + plasma bubble (full model)."""
    from models.hybrid_model import HybridPropagationModel

    print("=" * 55)
    print("  Scenario: Full model  (TID + Es + Bubble)")
    print("=" * 55)
    iono_p = _make_iono_params(enable_tid=True, enable_es=True, enable_bubble=True)
    model  = HybridPropagationModel(iono_p, _make_radar_params())
    modes, tau_ax, pd, main = model.compute(cfg.TX_POS, cfg.RX_POS)
    print(f"  Modes found: {len(modes)}")
    print_mode_table(modes)
    fig, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D – Full model')
    fig.savefig(os.path.join(_OUT, 'pd_full.png'), dpi=150, bbox_inches='tight')
    print("  Saved: output/pd_full.png\n")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    if not check_env():
        sys.exit(1)

    # ── Select scenario here ──────────────────────────────────────────────────
    # run_baseline()       # Part 2 complete -> can run
    # run_with_tid()       # Part 1 TID complete -> can run
    # run_with_es()        # Part 4 complete -> can run
    # run_with_bubble()    # Part 5 complete -> can run
    # run_full()           # Part 6 complete -> can run

    print("Scenarios are commented out.  Uncomment the desired scenario above.")
    print("Start with run_baseline() once Part 2 is implemented.")
