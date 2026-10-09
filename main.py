"""
main.py - Entry point for the HF hybrid propagation model.
All scenarios use 3-D Haselgrove ray tracing (Jones & Stephenson 1975).
TID irregularities: 3-component multi-direction superposition (config.TID_3D).
Plasma bubble: 3-D field-aligned Gaussian depletion (config.BUBBLE_3D).

Usage:
    conda run -n pytorch_cpu python main.py

Scenarios save tx_fan_*.png, p2p_modes_*.png, and pd_*.png to output/:
    run_baseline()        IRI only, 3-D RT
    run_with_tid()        IRI + 3-dir TID (3-D), lateral deflection plot
    run_with_es()         IRI + Es, 3-D RT + Hao 2017 power
    run_with_bubble()     IRI + 3-D field-aligned bubble
    run_full()            3-dir TID + Es + 3-D bubble (mode paths overlay)
    run_with_OX()         3-dir TID; AH3D if enabled, otherwise diagnostic post-proc
    run_radar_point_target() monostatic point-target radar, two-way observables
    run_radar_clutter()   blocked until legacy flux-tube coordinates are migrated
run_with_spreadf()    IRI + spread-F, 3-D RT
    run_tube_radar()      blocked until legacy tube coordinates are migrated

Test link: TX @ 30N 120E, RX @ 1169 km, f = 10 MHz.
"""
import sys
import os
import copy
import argparse
from contextlib import contextmanager
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import config as cfg
from viz.plot_utils import (print_mode_table, plot_tx_ray_fan, plot_p2p_modes,
                             plot_pd_spectrum, save_modes_csv, build_ray_info_text)
from models.pd_spectrum import build_pd_spectrum, identify_main_mode
from models.scenario_provider_factory import ScenarioProviderFactory

_OUT = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(_OUT, exist_ok=True)

RUN_PRODUCTS = {
    'ray_fan': True,
    'ray_fan_x': False,
    'p2p': False,
    'pd': False,
    'csv': False,
}
_PROVIDER_FACTORY = ScenarioProviderFactory(cfg)


# -- Shared helpers -------------------------------------------------------------

def _run_product(name):
    return bool(RUN_PRODUCTS.get(name, True))


def _save(fig, name):
    """Save figure to output/ and close it."""
    if name.startswith(('tx_fan', 'ray_fan')) and not _run_product('ray_fan'):
        plt.close(fig)
        print("  Skipped ray fan output: {}".format(name))
        return
    if name.startswith('p2p_modes') and not _run_product('p2p'):
        plt.close(fig)
        print("  Skipped P2P output: {}".format(name))
        return
    if name.startswith('pd_') and not _run_product('pd'):
        plt.close(fig)
        print("  Skipped P-D output: {}".format(name))
        return
    fig.savefig(os.path.join(_OUT, name), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: output/{}".format(name))


def _merge_dict(dst, src):
    for key, val in src.items():
        if isinstance(val, dict) and isinstance(dst.get(key), dict):
            _merge_dict(dst[key], val)
        else:
            dst[key] = val


@contextmanager
def _temporary_update2(overrides):
    old = copy.deepcopy(cfg.UPDATE2)
    try:
        _merge_dict(cfg.UPDATE2, overrides)
        yield
    finally:
        cfg.UPDATE2.clear()
        cfg.UPDATE2.update(old)


def _enable_main_reference_physics():
    """Enable physically richer runtime switches for main.py only."""
    from models.unified_hybrid_pipeline import reject_archived_physics
    reject_archived_physics()
    _merge_dict(cfg.UPDATE2, {
        'AH3D_ENABLE': True,
        'OX_MODE': 'ah3d',
        'RAY_FAN_INCLUDE_X': bool(RUN_PRODUCTS.get('ray_fan_x', False)),
        'ABSORPTION_COMPLEX': {
            'enable': True,
        },
        'PE_COUPLING': {
            'enable': False,
            'apply_to_power': False,
            'mode': 'local_2p5d',
        },
        'SCATTER_PD': {
            'apply_to_pd': False,
            'kernel': 'diagnostic_only',
        },
        'RADAR_CLUTTER': {
            'enable': True,
            'apply_to_power': True,
            'surface_type': 'sea',
            'sigma0_dB': cfg.RADAR.get('sigma0_ground_dB', -20.0),
            'radial_velocity_ms': 0.0,
            'clutter_class': 'surface_bragg',
        },
    })
    print("  Reference physics switches: AH3D and complex absorption enabled; PE power and scatter-PD disabled")


# -- 3-D model helpers ----------------------------------------------------------

def _scenario_provider(kind, *, components=None):
    configured_seed = cfg.SPREAD_F.get('seed', 0)
    seed = 0 if configured_seed is None else int(configured_seed)
    return _PROVIDER_FACTORY.create(
        kind, components=components, seed=seed, t_s=0.0, validation=True)


def _provider_slice(provider, z_km=None):
    z_axis = cfg.BG_Z if z_km is None else np.asarray(z_km, dtype=float)
    return provider.sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, z_axis, t=0.0)


def _make_n3d(provider):
    """Return the formal geocentric Provider consumed by ECEF AH3D."""
    from models.ne_provider_geocentric import GeocentricNeProvider3D
    if not isinstance(provider, GeocentricNeProvider3D):
        raise TypeError('formal main entry requires GeocentricNeProvider3D')
    return provider


def _use_ah3d_reference():
    return (
        bool(cfg.UPDATE2.get('AH3D_ENABLE', False)) and
        cfg.UPDATE2.get('OX_MODE', 'postproc') == 'ah3d'
    )


def _shoot_rays_fan_reference_3d(n3d, rt_config=None):
    """Build an O-wave AH3D fan; include X only through its explicit gate."""
    if not _run_product('ray_fan'):
        return []
    rt_config = rt_config or cfg.RT_3D
    if not _use_ah3d_reference():
        raise RuntimeError('formal main entry requires ECEF AH3D; no legacy fallback')
    from models.unified_hybrid_pipeline import (
        UnifiedPropagationRequest, run_unified_pipeline)
    request = UnifiedPropagationRequest(
        scenario='main_validation_fan', provider=n3d,
        tx_geocentric=cfg.TX_GEOCENTRIC, rx_geocentric=cfg.RX_GEOCENTRIC,
        bearing_deg=cfg.LINK_BEARING_DEG,
        earth_radius_km=cfg.EARTH_MODEL['radius_km'], solver='ah3d_ecef',
        rt_config=dict(rt_config), make_ray_fan=True, make_p2p=False, make_pd=False)
    return run_unified_pipeline(request).rays3d


def _find_all_rays_reference_3d(n3d, rays_fan=None, environment=None):
    """Use formal geocentric/ECEF AH3D O/X P2P without legacy fallback."""
    if not _run_product('p2p'):
        return []
    if not _use_ah3d_reference():
        raise RuntimeError('formal main entry requires ECEF AH3D; no legacy fallback')
    from models.unified_hybrid_pipeline import (
        UnifiedPropagationRequest, run_unified_pipeline)
    request = UnifiedPropagationRequest(
        scenario='main_validation_p2p', provider=n3d,
        tx_geocentric=cfg.TX_GEOCENTRIC, rx_geocentric=cfg.RX_GEOCENTRIC,
        bearing_deg=cfg.LINK_BEARING_DEG,
        earth_radius_km=cfg.EARTH_MODEL['radius_km'], solver='ah3d_ecef',
        environment=environment,
        make_ray_fan=False, make_p2p=True, make_pd=False)
    return run_unified_pipeline(request).modes


def _rays3d_to_fan(rays3d):
    """Formal plots consume explicit geocentric trajectories directly."""
    return rays3d


def _tag_communication_modes(modes):
    """Mark one-way communication observables explicitly."""
    for m in modes:
        m['tau_oneway_ms'] = float(m.get('tau_ms', 0.0))
        m['tau_2way_ms'] = None
        m['range_km'] = float(m.get('group_path_km', 0.0))
        m['radar_observable'] = 'none'
        m.setdefault('physics_basis', '3d_haselgrove_friis')
    return modes


def _tag_point_target_radar_modes(modes):
    """Mark monostatic point-target radar observables explicitly."""
    for m in modes:
        tau_oneway = float(m.get('tau_ms', 0.0))
        m['tau_oneway_ms'] = tau_oneway
        m['tau_2way_ms'] = 2.0 * tau_oneway
        m['range_km'] = float(m.get('group_path_km', 0.0))
        m['radar_observable'] = 'two_way_point_target'
        m['radar_equation'] = 'monostatic_point_target_R4'
        m['physics_basis'] = 'Skolnik_2008_point_target_radar_equation'
    return modes


def _build_pd(modes, method='gaussian', delay_key='tau_ms', scatter_pd=False):
    """Build P-D spectrum; radar callers pass delay_key='tau_2way_ms'."""
    if not _run_product('pd') or not modes:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    scatter_cfg = cfg.UPDATE2.get('SCATTER_PD', {})
    if scatter_pd and bool(scatter_cfg.get('apply_to_pd', False)):
        with _temporary_update2({'SCATTER_PD': {'apply_to_pd': True, 'kernel': 'physical'}}):
            return _build_pd(modes, method=method, delay_key=delay_key, scatter_pd=False)
    if delay_key == 'tau_ms':
        return build_pd_spectrum(modes, method=method)
    shifted = []
    for m in modes:
        mm = dict(m)
        mm['tau_ms'] = float(m.get(delay_key, m.get('tau_ms', 0.0)))
        if delay_key == 'tau_2way_ms':
            if 'tau_spread_ms' in mm:
                mm['tau_spread_ms'] = 2.0 * float(mm.get('tau_spread_ms', 0.0))
            if 'delta_tau_ms' in mm:
                mm['delta_tau_ms'] = 2.0 * float(mm.get('delta_tau_ms', 0.0))
        shifted.append(mm)
    return build_pd_spectrum(shifted, method=method)


def _identify_main(tau_ax, pd, modes, delay_key='tau_ms'):
    """Identify main mode on the same delay axis used by the plotted P-D spectrum."""
    if delay_key == 'tau_ms':
        return identify_main_mode(tau_ax, pd, modes)
    shifted = []
    for m in modes:
        mm = dict(m)
        mm['tau_ms'] = float(m.get(delay_key, m.get('tau_ms', 0.0)))
        shifted.append(mm)
    main_s, ranked_s, summary = identify_main_mode(tau_ax, pd, shifted)
    if main_s is None:
        return None, [], summary
    main = min(modes, key=lambda m: abs(float(m.get(delay_key, m.get('tau_ms', 0.0))) - main_s['tau_ms']))
    ranked = [
        min(modes, key=lambda m: abs(float(m.get(delay_key, m.get('tau_ms', 0.0))) - r['tau_ms']))
        for r in ranked_s
    ]
    summary['main_tau_ms'] = float(main.get(delay_key, main.get('tau_ms', 0.0)))
    summary['main_label'] = main.get('label', '')
    summary['main_Pr_dBW'] = main.get('Pr_dBW', -999.0)
    return main, ranked, summary


def _save_modes_csv(scenario, modes, tau_ax, pd, pd_delay_key='tau_ms'):
    return save_modes_csv(
        scenario, modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT,
        pd_delay_key=pd_delay_key)


def _plot_pd(tau_ax, pd, modes, title, filename, delay_key='tau_ms'):
    if not _run_product('pd') or len(tau_ax) == 0 or len(pd) == 0:
        print("  Skipped P-D plot: {}".format(filename))
        return
    xlabel = 'Two-way Group Delay  (ms)' if delay_key == 'tau_2way_ms' else 'One-way Group Delay  (ms)'
    info_text = _make_pd_info(modes, delay_key=delay_key)
    fig, _ = plot_pd_spectrum(
        tau_ax, pd, modes, title=title, delay_key=delay_key, xlabel=xlabel,
        info_text=info_text)
    _save(fig, filename)


def _add_power_3d(modes, n3d=None, radar_mode=False, sigma_rcs=1.0):
    """Fill Pr_W, Pr_dBW, delta_tau_ms in 3-D P2P mode dicts."""
    from utils import free_space_loss_dB, to_dBW, radar_equation_W
    from models.unified_hybrid_pipeline import reject_archived_physics
    reject_archived_physics()
    for m in modes:
        if 'absorption_applied_once' in m:
            continue
        gp = float(m['group_path_km'])
        if radar_mode:
            m['Pr_W'] = radar_equation_W(cfg.PT_W, cfg.GT, cfg.GR,
                                          cfg.FREQ_MHZ, gp, sigma_rcs)
        else:
            L_dB      = free_space_loss_dB(gp, cfg.FREQ_MHZ)
            m['Pr_W'] = cfg.PT_W * cfg.GT * cfg.GR * 10.0**(-L_dB / 10.0)

        abs_oneway_dB = 0.0
        abs_total_dB = 0.0
        abs_model = 'none'
        abs_diag = {'abs_status': 'disabled'}
        complex_cfg = cfg.UPDATE2.get('ABSORPTION_COMPLEX', {})
        if complex_cfg.get('enable', False):
            from models.absorption_collision import complex_absorption_ecef_ray_dB
            if n3d is None:
                raise ValueError('ECEF absorption requires the formal geocentric Provider')
            abs_diag = complex_absorption_ecef_ray_dB(
                m, n3d, complex_cfg, two_way=bool(radar_mode))
            abs_oneway_dB = float(abs_diag.get('abs_complex_dB', 0.0))
            abs_total_dB = float(abs_diag.get('abs_total_dB', abs_oneway_dB))
            abs_model = 'appleton_hartree_collision_path'
            m['Pr_W'] *= 10.0**(-abs_total_dB / 10.0)
        m['abs_dB']        = float(abs_oneway_dB)
        m['abs_total_dB']  = float(abs_total_dB)
        m['abs_model']     = abs_model
        m['abs_complex_dB'] = float(abs_diag.get('abs_complex_dB', 0.0))
        m['abs_status']    = abs_diag.get('abs_status', '')
        m['collision_model'] = abs_diag.get('collision_model', '')
        m['abs_index_model'] = abs_diag.get('abs_index_model', '')
        m['Pr_dBW']       = to_dBW(m['Pr_W'])
        m['delta_tau_ms'] = 0.0
    if radar_mode:
        _tag_point_target_radar_modes(modes)
    else:
        _tag_communication_modes(modes)

def _add_ray_pe_reference_3d(modes):
    """Record the disabled legacy PE boundary without importing legacy code."""
    from models.unified_hybrid_pipeline import reject_archived_physics
    reject_archived_physics()
    for mode in modes:
        mode['pe_solver_dimensionality'] = 'off'
        mode['pe_validity'] = 'archived_local_2p5d_not_in_production_chain'
        mode['pe_apply_to_power'] = False
    return modes


def _add_scatter_3d(modes, n3d):
    """Attach Born/Rytov/MPS scatter diagnostics to 3-D modes."""
    if not cfg.SCATTER.get('enable', True):
        return modes
    from models.scatter_coupler import apply_scatter_regime_3d
    return apply_scatter_regime_3d(
        modes, n3d, cfg.SCATTER, cfg.FREQ_MHZ, cfg.PT_W, cfg.GT, cfg.GR,
        enable_mps=cfg.SCATTER.get('enable_mps', False))


def _extract_es_grazing_3d(ray3d, h_Es_km):
    """
    Find the group-ray grazing angle at the first upward crossing of h_Es_km
    using spherical physical height and the local radial Up direction.
    Returns None if no crossing is found.
    """
    traj = np.asarray(ray3d['trajectory_lat_lon_height'], dtype=float)
    tangents = np.asarray(ray3d['trajectory_group_tangent_ecef'], dtype=float)
    from utils import enu_basis, enu_basis_wgs84
    for i in range(1, len(traj)):
        s0, s1 = traj[i-1], traj[i]
        z0p, z1p = float(s0[2]), float(s1[2])
        if z0p < h_Es_km <= z1p:
            t   = (h_Es_km - z0p) / max(z1p - z0p, 1e-6)
            direction = tangents[i - 1] + t * (tangents[i] - tangents[i - 1])
            lat = float(s0[0] + t * (s1[0] - s0[0]))
            lon = float(s0[1] + t * (s1[1] - s0[1]))
            up = (enu_basis_wgs84(lat, lon)[:, 2]
                  if ray3d.get('earth_model') == 'WGS84'
                  else enu_basis(lat, lon, cfg.RE_KM)[:, 2])
            return float(np.degrees(np.arcsin(np.clip(
                abs(float(np.dot(direction, up)))
                / max(float(np.linalg.norm(direction)), 1e-30),
                0.0, 1.0))))
    return None


# -- Iono / radar / link helpers ------------------------------------------------

def _make_iono_params(enable_tid=False, enable_es=False,
                       enable_bubble=False, enable_spread_f=False):
    """
    Assemble iono_params dict.
    TID always uses multi-direction 3-component TID_3D.
    All environment fields are sampled from the shared 3-D Provider Factory.
    """
    return {
        'iri_params':      {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':      {**cfg.TID_3D,     'enable': enable_tid},
        'es_params':       {**cfg.ES,         'enable': enable_es},
        'bubble_params':   {**cfg.BUBBLE,     'enable': enable_bubble},
        'spread_f_params': {**cfg.SPREAD_F,   'enable': enable_spread_f},
    }


def _make_link_info():
    return {
        'TX_LAT':   cfg.TX_LAT,
        'TX_LON':   cfg.TX_LON,
        'freq_MHz': cfg.FREQ_MHZ,
        'Pt_W':     cfg.PT_W,
        'RX_km':    cfg.RX_RANGE,
        'IRI_DT':   cfg.IRI_DT,
        'IRI_LAT':  cfg.IRI_LAT,
        'IRI_LON':  cfg.IRI_LON,
        'earth_model': cfg.EARTH_MODEL['name'],
        'earth_radius_km': cfg.EARTH_MODEL['radius_km'],
        'coordinate_model': 'geocentric lat/lon/height -> ECEF state',
        'oneway_or_twoway': 'one_way',
        'reflection_coupling': 'on',
        'ray_pe_coupling': 'archived_local_2p5d_not_in_production_chain',
        'scatter_coupling': ('on' if cfg.UPDATE2.get(
            'SCATTER_PD', {}).get('apply_to_pd', False) else 'diagnostic_only'),
    }


def _make_pert_info(enable_tid=False, enable_es=False, enable_bubble=False):
    """Annotation dict for ray-fan plots (active perturbations only)."""
    p = {}
    if enable_tid:
        # Show dominant (first) component; n_comp count goes in title
        p['TID'] = {
            'amplitude':   cfg.TID_3D['amplitude_list'][0],
            'lambda_h_km': cfg.TID_3D['lambda_h_km_list'][0],
            'T_s':         cfg.TID_3D['period_s_list'][0],
        }
    if enable_es:
        p['Es'] = {'foEs_MHz': cfg.ES['foEs_MHz'],
                   'h_Es_km':  cfg.ES['h_Es_km']}
    if enable_bubble:
        # sigma_perp_km serves as cross-field scale, analogous to Lx_km
        p['Bubble'] = {
            'z0_km':     cfg.BUBBLE_3D['z0_km'],
            'Lx_km':     cfg.BUBBLE_3D['sigma_perp_km'],
            'delta_max': cfg.BUBBLE_3D['delta_max'],
        }
    return p


def _make_spreadf_pert_info():
    p = _make_pert_info()
    p['SpreadF'] = {
        'Cs': cfg.SPREAD_F.get('Cs', 0.0),
        'p': cfg.SPREAD_F.get('p', 0.0),
        'h_screen_km': cfg.SPREAD_F.get('h_screen_km', 0.0),
    }
    return p


def _make_pd_info(modes, delay_key='tau_ms', extra=None):
    """Compact semantic annotation for P-D plots."""
    if delay_key == 'tau_2way_ms':
        obs = 'two-way radar observable'
    else:
        obs = 'one-way communication observable'
    sample = next((m for m in modes if isinstance(m, dict)), {})
    pe_state = 'on' if any(bool(m.get('pe_apply_to_power', False)) for m in modes) else 'off'
    scatter_state = 'on' if any(float(m.get('sigma_tau_ms', 0.0) or 0.0) > 0.0 for m in modes) else 'diagnostic/off'
    abs_state = sample.get('abs_model', 'none')
    radar_obs = sample.get('radar_observable', 'none')
    ox_solver = sample.get('ox_solver', 'n/a')
    basis = sample.get('physics_basis', '3d_haselgrove_friis')
    lines = [
        obs,
        'delay_key={}'.format(delay_key),
        'earth=sphere R={:.1f}km'.format(cfg.EARTH_MODEL['radius_km']),
        'coordinates=geocentric input / ECEF state',
        'basis={}'.format(basis),
        'absorption={}'.format(abs_state),
        'Ray-PE power={}'.format(pe_state),
        'scatter P-D={}'.format(scatter_state),
    ]
    if radar_obs != 'none':
        lines.append('radar={}'.format(radar_obs))
    if sample.get('clutter_pd_model'):
        lines.append('clutter_pd={}'.format(sample.get('clutter_pd_model')))
        lines.append('Bragg fd={:.2f}/{:.2f} Hz'.format(
            float(sample.get('bragg_doppler_minus_hz', 0.0)),
            float(sample.get('bragg_doppler_plus_hz', 0.0))))
        lines.append('gate={:.0f} km bins={}/{}'.format(
            float(sample.get('range_gate_span_km', 0.0)),
            int(sample.get('n_range_bins', 0)),
            int(sample.get('n_doppler_bins', 0))))
    if ox_solver != 'n/a':
        lines.append('O/X={}'.format(ox_solver))
    if extra:
        lines.extend(str(x) for x in extra)
    return '\n'.join(lines)


def _print_main(main):
    if main:
        pr_str = ("{:.1f} dBW".format(main['Pr_dBW'])
                  if 'Pr_dBW' in main else "n/a")
        if main.get('radar_observable', 'none') != 'none' and main.get('tau_2way_ms') is not None:
            print("  Main mode: {}  tau2 = {:.3f} ms  Pr = {}".format(
                main['label'], main['tau_2way_ms'], pr_str))
        else:
            print("  Main mode: {}  tau = {:.3f} ms  Pr = {}".format(
                main['label'], main['tau_ms'], pr_str))
    else:
        print("  Main mode: None")


def _print_summary(summary):
    print("  n_modes={} main={} tau={:.3f}ms Pr={:.1f}dBW F_focus={:.2f}".format(
        summary['n_modes'], summary['main_label'],
        summary['main_tau_ms'], summary['main_Pr_dBW'],
        summary['main_F_focus']))
    if summary.get('mode_pairs'):
        for p in summary['mode_pairs']:
            print("    O/X: {}/{} dTau={:.3f}ms dPr={:.1f}dB".format(
                p['label_O'], p['label_X'],
                p['delta_tau_OX_ms'], p['delta_Pr_OX_dB']))


def _plot_mode_paths(modes, Ne_2d, x_km, z_km, title='3-D Mode Paths',
                     link_info=None, pert_info=None):
    """
    Mode path overlay on Ne background.
    Accepts formal geocentric/ECEF modes and historical 2-D reference modes.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    pcm = ax.pcolormesh(x_km, z_km, Ne_2d.T / 1e10,
                        cmap='jet', shading='auto', alpha=0.5)
    plt.colorbar(pcm, ax=ax, label='Ne  [x10^10 m^-3]',
                 fraction=0.046, pad=0.04)
    colors = cm.tab10(np.linspace(0, 1, max(len(modes), 1)))
    for m, c in zip(modes, colors):
        if m.get('trajectory_lat_lon_height') is not None:
            from viz.plot_utils import _mode_xz
            xs, zs = _mode_xz(m)
        elif 'points' in m:
            pts = np.asarray(m['points'])
            xs, zs = pts[:, 0].tolist(), pts[:, 1].tolist()
        else:
            continue
        ax.plot(xs, zs, '-', color=c, lw=1.8,
                label='{} ({:.2f} ms)'.format(m.get('label', '?'), m['tau_ms']))
    ax.axvline(0.0, color='g', lw=1.2, ls='--', label='TX')
    ax.axvline(cfg.RX_RANGE, color='b', lw=1.2, ls='--', label='RX')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0, 500])
    ax.legend(fontsize=7, ncol=2, loc='upper left')
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    if link_info is not None:
        info_text = build_ray_info_text(link_info, Ne_2d, z_km, pert_info)
        ax.text(0.99, 0.97, info_text,
                transform=ax.transAxes, fontsize=6.5, va='top', ha='right',
                family='monospace',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          alpha=0.88, edgecolor='gray', lw=0.8))
    plt.tight_layout()
    return fig, ax


# -- Environment check ---------------------------------------------------------

def check_env() -> bool:
    """Verify required packages and print versions."""
    ok = True
    print("=" * 55)
    print("  Environment check (3-D mode)")
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
    print("  TID_3D    : {} components  az={}".format(
        cfg.TID_3D['n_components'], cfg.TID_3D['az_deg_list']))
    print()
    return ok


def _plot_x_limits():
    return (float(cfg.BG_X[0]), float(cfg.BG_X[-1]))


def _save_p2p_modes_plot(scenario, modes, Ne_2d, pert_info=None, title=None):
    link_info = _make_link_info()
    if any(m.get('tau_2way_ms') is not None for m in modes):
        link_info['oneway_or_twoway'] = 'two_way'
    if modes:
        sample = modes[0]
        link_info['ray_pe_coupling'] = sample.get(
            'ray_pe_dimensionality', link_info['ray_pe_coupling'])
        link_info['scatter_coupling'] = sample.get(
            'scatter_status', link_info['scatter_coupling'])
    fig, _ = plot_p2p_modes(
        modes, Ne_2d, cfg.BG_X, cfg.BG_Z,
        title=title or 'P2P Modes - {}'.format(scenario),
        link_info=link_info, pert_info=pert_info,
        xlim=_plot_x_limits())
    _save(fig, 'p2p_modes_{}.png'.format(scenario))


# -- Scenario 0: IRI baseline ---------------------------------------------------

def run_baseline():
    """IRI background only, 3-D ray tracer.
    Saves tx_fan_baseline.png, pd_baseline.png."""
    print("=" * 55)
    print("  Scenario 0: Baseline (IRI only, 3-D RT)")
    print("=" * 55)

    provider = _scenario_provider('baseline')
    Ne_2d = _provider_slice(provider)
    print("  Ne field: shape {}  max = {:.2e} m^-3".format(
        Ne_2d.shape, float(Ne_2d.max())))

    n3d    = _make_n3d(provider)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI baseline (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(),
                          xlim=_plot_x_limits())
    _save(fig, 'tx_fan_baseline.png')

    modes = _find_all_rays_reference_3d(n3d, environment='baseline')
    _add_power_3d(modes, n3d=n3d)
    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('baseline', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('baseline', modes, Ne_2d)
    _plot_pd(tau_ax, pd, modes, title='P-D - IRI baseline (3-D)',
             filename='pd_baseline.png')


# -- Scenario 1: IRI + multi-direction TID -------------------------------------

def run_with_tid():
    """IRI + 3-component multi-direction TID, 3-D ray tracer.
    Saves tx_fan_tid.png, pd_tid.png, ray_deflect_tid.png."""
    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    print("  Scenario 1: IRI + {}-dir TID (3-D RT)".format(n_comp))
    print("=" * 55)
    print("  az_list={} deg  amp_list={}".format(
        cfg.TID_3D['az_deg_list'], cfg.TID_3D['amplitude_list']))

    provider = _scenario_provider('tid')
    Ne_2d = _provider_slice(provider)

    # TID baked into Ne_2d; grad_n2() fast path active (no runtime TID interp)
    n3d    = _make_n3d(provider)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + TID ({}-dir, 3-D)'.format(n_comp),
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_tid=True),
                          xlim=_plot_x_limits())
    _save(fig, 'tx_fan_tid.png')

    # Mini-fan with 3-D TID active to show lateral y-deflection (cross-path TID)
    _rt_mini  = {**cfg.RT_3D, 'n_fan': 7, 'max_steps': 3000}
    n3d_tid   = n3d
    rays_mini = _shoot_rays_fan_reference_3d(n3d_tid, _rt_mini)
    betas = [float(r['beta_deg_tx']) for r in rays_mini]
    ys = []
    for ray in rays_mini:
        landing = ray.get('landing_lat_lon_height')
        ys.append(float('nan') if landing is None else float(
            n3d_tid.derive_link_coordinates(*landing)[1]))
    fig3, ax3 = plt.subplots(figsize=(8, 4))
    ax3.plot(betas, ys, 'b.-', lw=1.2)
    ax3.axhline(0.0, color='k', lw=0.6, ls=':')
    ax3.set_xlabel('TX Elevation (deg)')
    ax3.set_ylabel('Lateral landing y (km)')
    ax3.set_title('3-D TID lateral deflection ({}-component)'.format(n_comp))
    ax3.grid(True, alpha=0.3)
    plt.tight_layout()
    _save(fig3, 'ray_deflect_tid.png')

    modes = _find_all_rays_reference_3d(n3d, environment='tid')
    _add_power_3d(modes, n3d=n3d)
    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('tid', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('tid', modes, Ne_2d, _make_pert_info(enable_tid=True))
    _plot_pd(tau_ax, pd, modes,
             title='P-D - IRI + TID ({}-dir, 3-D)'.format(n_comp),
             filename='pd_tid.png')


# -- Scenario 2: IRI + Es -------------------------------------------------------

def run_with_es():
    """IRI + Es layer. 3-D RT + Es three-segment power correction (Hao 2017).
    Saves tx_fan_es.png, pd_es.png."""
    from models.es_model import EsLayerModel
    from utils import to_dBW

    print("=" * 55)
    print("  Scenario 2: IRI + Es (3-D RT + Hao 2017)")
    print("=" * 55)
    print("  foEs = {} MHz  f = {} MHz  foEs/f = {:.2f}".format(
        cfg.ES['foEs_MHz'], cfg.FREQ_MHZ,
        cfg.ES['foEs_MHz'] / cfg.FREQ_MHZ))

    provider = _scenario_provider('es')
    Ne_2d = _provider_slice(provider)
    n3d    = _make_n3d(provider)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + Es (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_es=True),
                          xlim=_plot_x_limits())
    _save(fig, 'tx_fan_es.png')

    modes = _find_all_rays_reference_3d(n3d, environment='es')
    _add_power_3d(modes, n3d=n3d)

    es    = EsLayerModel(**{k: v for k, v in cfg.ES.items() if k != 'enable'})
    h_Es  = float(cfg.ES['h_Es_km'])
    for m in modes:
        # Only apply Es correction to rays that reflect at or below Es layer
        if m.get('h_reflect_km', 999.0) > h_Es + 20.0:
            continue
        theta_deg = _extract_es_grazing_3d(m, h_Es)
        if theta_deg is None:
            continue
        theta_rad = np.radians(max(theta_deg, 0.5))
        es_res    = es.compute_power(m['Pr_W'], cfg.GT, cfg.GR,
                                      cfg.FREQ_MHZ, m['group_path_km'] / 2.0,
                                      theta_rad)
        m['Pr_W']         = float(es_res['Pr_W'])
        m['Pr_dBW']       = to_dBW(m['Pr_W'])
        m['delta_tau_ms'] = float(es_res['delta_tau_ms'])
        m['label']        = 'Es_' + es_res['mode']

    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('es', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('es', modes, Ne_2d, _make_pert_info(enable_es=True))
    _plot_pd(tau_ax, pd, modes, title='P-D - IRI + Es (3-D)',
             filename='pd_es.png')


# -- Scenario 3: IRI + 3-D field-aligned bubble --------------------------------

def run_with_bubble():
    """IRI + 3-D field-aligned plasma bubble (Ma et al. 2023, Martinis 2023).
    Saves tx_fan_bubble.png, pd_bubble.png."""
    print("=" * 55)
    print("  Scenario 3: IRI + 3-D Field-Aligned Bubble")
    print("=" * 55)
    b = cfg.BUBBLE_3D
    print("  x0={} km  z0={} km  s_perp={} km  s_para={} km  dmax={}".format(
        b['x0_km'], b['z0_km'], b['sigma_perp_km'], b['sigma_para_km'], b['delta_max']))

    provider = _scenario_provider('bubble')
    Ne_2d = _provider_slice(provider)
    n3d     = _make_n3d(provider)
    rays3d  = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + 3-D Field-Aligned Bubble',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_bubble=True),
                          xlim=_plot_x_limits())
    _save(fig, 'tx_fan_bubble.png')

    modes = _find_all_rays_reference_3d(n3d, environment='bubble')
    _add_power_3d(modes, n3d=n3d)
    _add_ray_pe_reference_3d(modes)
    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('bubble', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('bubble', modes, Ne_2d, _make_pert_info(enable_bubble=True))
    _plot_pd(tau_ax, pd, modes, title='P-D - IRI + 3-D Bubble',
             filename='pd_bubble.png')


# -- Scenario 4: Full model -----------------------------------------------------

def run_full():
    """Full 3-D model: 3-dir TID + Es + 3-D field-aligned bubble.
    Saves ray_fan_full.png (mode paths), pd_full.png."""
    from models.es_model import EsLayerModel
    from utils import to_dBW

    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    print("  Scenario 4: Full 3-D ({}-dir TID + Es + Bubble)".format(n_comp))
    print("=" * 55)


    provider = _scenario_provider(
        'baseline', components=('tid', 'es', 'bubble'))
    Ne_2d = _provider_slice(provider)

    n3d          = _make_n3d(provider)
    _rays_seed   = _shoot_rays_fan_reference_3d(n3d)
    fig_fan, _ = plot_tx_ray_fan(_rays3d_to_fan(_rays_seed), Ne_2d, cfg.BG_X, cfg.BG_Z,
                                 title='TX Ray Fan - Full 3-D',
                                 link_info=_make_link_info(),
                                 pert_info=_make_pert_info(enable_tid=True, enable_es=True, enable_bubble=True),
                                 xlim=_plot_x_limits())
    _save(fig_fan, 'tx_fan_full.png')
    modes  = _find_all_rays_reference_3d(n3d, environment='full')
    _add_power_3d(modes, n3d=n3d)

    # Es three-segment correction
    es   = EsLayerModel(**{k: v for k, v in cfg.ES.items() if k != 'enable'})
    h_Es = float(cfg.ES['h_Es_km'])
    for m in modes:
        if m.get('h_reflect_km', 999.0) > h_Es + 20.0:
            continue
        theta_deg = _extract_es_grazing_3d(m, h_Es)
        if theta_deg is None:
            continue
        theta_rad = np.radians(max(theta_deg, 0.5))
        es_res = es.compute_power(m['Pr_W'], cfg.GT, cfg.GR,
                                   cfg.FREQ_MHZ, m['group_path_km'] / 2.0,
                                   theta_rad)
        m['Pr_W']         = float(es_res['Pr_W'])
        m['Pr_dBW']       = to_dBW(m['Pr_W'])
        m['delta_tau_ms'] = float(es_res['delta_tau_ms'])
        m['label']        = 'Es_' + es_res['mode']

    _add_ray_pe_reference_3d(modes)
    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    tau_ax, pd = _build_pd(modes, scatter_pd=True)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('full', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('full', modes, Ne_2d, _make_pert_info(enable_tid=True, enable_es=True, enable_bubble=True),
                         title='P2P Modes - Full 3-D')

    _plot_pd(tau_ax, pd, modes,
             title='P-D - Full 3-D ({}-dir TID + Es + Bubble)'.format(n_comp),
             filename='pd_full.png')


def run_with_OX():
    """
    IRI + TID with separate 3-D Appleton-Hartree O/X trajectories.
    Saves tx_fan_OX.png, pd_OX.png, modes_OX.csv.
    """
    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    if not _use_ah3d_reference():
        raise RuntimeError('OX validation requires AH3D; obsolete 2-D AH post-processing was removed')
    ox_basis = 'AH3D'
    print("  Scenario OX: {}-dir TID + O/X ({})".format(n_comp, ox_basis))
    print("=" * 55)
    print("  fH = {} MHz   dip = {} deg   decl = {} deg".format(
        cfg.GEOMAG['fH_MHz'], cfg.GEOMAG['dip_deg'], cfg.GEOMAG['decl_deg']))

    provider = _scenario_provider('tid')
    Ne_2d = _provider_slice(provider)
    n3d    = _make_n3d(provider)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))
    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - O/X {} (3-D + {}-dir TID)'.format(ox_basis, n_comp),
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_tid=True),
                          xlim=_plot_x_limits())
    _save(fig, 'tx_fan_OX.png')

    modes = _find_all_rays_reference_3d(n3d, environment='OX')
    _add_power_3d(modes, n3d=n3d)
    _add_scatter_3d(modes, n3d)

    print("  P2P modes: {} ({})".format(len(modes), ox_basis))
    print_mode_table(modes)
    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = _save_modes_csv('OX_ah3d', modes, tau_ax, pd)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('OX_ah3d', modes, Ne_2d, _make_pert_info(enable_tid=True))
    _plot_pd(tau_ax, pd, modes,
             title='P-D - O/X {} ({}-dir TID)'.format(ox_basis, n_comp),
             filename='pd_OX_ah3d.png')


# -- Scenario R1/R2: OTH radar observables ------------------------------------

def run_radar_point_target():
    """Monostatic point-target radar, IRI background, 3-D RT.

    This is not surface clutter. It uses the standard monostatic point-target
    R^-4 radar equation (Skolnik, Radar Handbook, 3rd ed., 2008).
    Saves tx_fan_radar_point.png, pd_radar_point.png.
    """
    print("=" * 55)
    print("  Scenario R1: Monostatic point-target radar (3-D RT)")
    print("=" * 55)
    print("  sigma_rcs = {} m^2   bearing = {} deg".format(
        cfg.RADAR['sigma_rcs_m2'], cfg.LINK_BEARING_DEG))

    provider = _scenario_provider('baseline')
    Ne_2d = _provider_slice(provider)
    n3d    = _make_n3d(provider)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))
    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - Monostatic Point-Target Radar (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(),
                          xlim=_plot_x_limits())
    _save(fig, 'tx_fan_radar_point.png')

    modes = _find_all_rays_reference_3d(
        n3d, rays_fan=rays3d, environment='radar_point_target')
    _add_power_3d(modes, n3d=n3d, radar_mode=True,
                  sigma_rcs=float(cfg.RADAR['sigma_rcs_m2']))
    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = _build_pd(modes, delay_key='tau_2way_ms')
    main, _, _ = _identify_main(tau_ax, pd, modes, delay_key='tau_2way_ms')
    _print_main(main)
    csv_path = _save_modes_csv('radar_point_target', modes, tau_ax, pd,
                               pd_delay_key='tau_2way_ms')
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('radar_point_target', modes, Ne_2d,
                         title='P2P Modes - Monostatic Point-Target Radar')
    _plot_pd(tau_ax, pd, modes,
             title='P-D - Monostatic Point-Target Radar (two-way delay)',
             filename='pd_radar_point.png', delay_key='tau_2way_ms')


def run_radar_clutter():
    raise RuntimeError(
        'C0R formal entry rejects the unported link-local flux-tube solver; '
        'migrate it only when A4 validation requires this product')


def run_radar_baseline():
    """Backward-compatible alias for the point-target radar scenario."""
    return run_radar_point_target()


# -- Scenario SF: Spread-F ------------------------------------------------------

def run_with_spreadf():
    """IRI + spread-F only (Rino 1979 Ne perturbation), 3-D RT.
    Saves ray_fan_spreadf.png (mode paths), pd_spreadf.png."""
    print("=" * 55)
    print("  Scenario SF: Spread-F only (3-D RT)")
    print("=" * 55)
    sf = cfg.SPREAD_F
    print("  Cs={} p={} h_screen={} km L0={} km".format(
        sf['Cs'], sf['p'], sf['h_screen_km'], sf['L0_km']))

    provider = _scenario_provider('spreadf')
    Ne_2d = _provider_slice(provider)
    n3d          = _make_n3d(provider)
    _rays_seed   = _shoot_rays_fan_reference_3d(n3d)
    fig_fan, _ = plot_tx_ray_fan(_rays3d_to_fan(_rays_seed), Ne_2d, cfg.BG_X, cfg.BG_Z,
                                 title='TX Ray Fan - Spread-F',
                                 link_info=_make_link_info(),
                                 pert_info=_make_spreadf_pert_info(),
                                 xlim=_plot_x_limits())
    _save(fig_fan, 'tx_fan_spreadf.png')
    modes  = _find_all_rays_reference_3d(n3d, environment='spreadf')
    _add_power_3d(modes, n3d=n3d)
    _add_ray_pe_reference_3d(modes)
    _add_scatter_3d(modes, n3d)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = _build_pd(modes, scatter_pd=True)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('spreadf', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('spreadf', modes, Ne_2d, _make_spreadf_pert_info(),
                         title='P2P Modes - Spread-F')

    _plot_pd(tau_ax, pd, modes, title='P-D - Spread-F only (3-D scattering)',
             filename='pd_spreadf.png')


# -- Scenario TR: 3-D flux-tube OTH radar --------------------------------------

def run_tube_radar():
    raise RuntimeError(
        'C0R formal entry rejects the unported link-local tube solver; '
        'migrate it only when A4 validation requires this product')


# -- Entry point ---------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AH3D scenario validation runner')
    parser.add_argument(
        '--ray-fan-x', action='store_true',
        help='Include and distinguish X-mode rays in ray fans; default is O only')
    args = parser.parse_args()
    RUN_PRODUCTS['ray_fan_x'] = bool(args.ray_fan_x)

    if not check_env():
        sys.exit(1)

    _enable_main_reference_physics()

    run_baseline()
    run_with_tid()
    run_with_es()
    run_with_bubble()
    run_full()
    run_with_OX()
    run_radar_point_target()
    run_with_spreadf()




