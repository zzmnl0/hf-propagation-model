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
    run_radar_clutter()   single-station OTH surface clutter, two-way observables
run_with_spreadf()    IRI + spread-F, 3-D RT
    run_tube_radar()      OTH radar, 3-D flux-tube tracer (Coleman 1997/1998)

Test link: TX @ 30N 120E, RX @ 1169 km, f = 10 MHz.
"""
import sys
import os
import copy
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
from models.hybrid_model import build_pd_spectrum, identify_main_mode

_OUT = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(_OUT, exist_ok=True)


# -- Shared helpers -------------------------------------------------------------

def _save(fig, name):
    """Save figure to output/ and close it."""
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
    _merge_dict(cfg.UPDATE2, {
        'AH3D_ENABLE': True,
        'OX_MODE': 'ah3d',
        'ABSORPTION_COMPLEX': {
            'enable': True,
            'mode': 'reference',
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

def _make_n3d(Ne_2d, tid_params=None, bubble_3d_params=None):
    """Build RefractiveIndex3D with optional TID and 3-D field-aligned bubble."""
    from models.ray_tracer_3d import RefractiveIndex3D
    return RefractiveIndex3D(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                              tid_params=tid_params,
                              bubble_3d_params=bubble_3d_params,
                              earth_flat=True)


def _use_ah3d_reference():
    return (
        bool(cfg.UPDATE2.get('AH3D_ENABLE', False)) and
        cfg.UPDATE2.get('OX_MODE', 'postproc') == 'ah3d'
    )


def _shoot_rays_fan_reference_3d(n3d, rt_config=None):
    """Use AH3D O/X ray fan when enabled; otherwise keep isotropic legacy."""
    rt_config = rt_config or cfg.RT_3D
    if _use_ah3d_reference():
        from models.ray_tracer_ah_3d import RefractiveIndexAH3D, shoot_rays_fan_ah_3d
        rays = []
        for wave_mode in ('O', 'X'):
            n_ah = RefractiveIndexAH3D(
                n3d, cfg.FREQ_MHZ, wave_mode=wave_mode, geomag=cfg.GEOMAG)
            rays.extend(shoot_rays_fan_ah_3d(
                cfg.TX_POS, n_ah, cfg.FREQ_MHZ, rt_config))
        return rays
    from models.ray_tracer_3d import shoot_rays_fan_3d
    return shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, rt_config)


def _find_all_rays_reference_3d(n3d, rays_fan=None):
    """Use AH3D O/X P2P as the research mainline; iso3d is legacy fallback."""
    if _use_ah3d_reference():
        from models.point_to_point_ah_3d import find_all_rays_ah_3d
        modes = find_all_rays_ah_3d(
            cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
            wave_mode='both', geomag=cfg.GEOMAG, rays_fan=rays_fan)
        for m in modes:
            m['ox_solver'] = 'ah3d_hamiltonian'
            m['physics_basis'] = '3d_appleton_hartree_hamiltonian'
            m['freq_MHz'] = float(cfg.FREQ_MHZ)
        return modes
    from models.point_to_point_3d import find_all_rays_3d
    modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D)
    for m in modes:
        m['ox_solver'] = 'iso3d_legacy'
        m['physics_basis'] = m.get('physics_basis', '3d_haselgrove_isotropic_legacy')
        m['freq_MHz'] = float(cfg.FREQ_MHZ)
    return modes


def _rays3d_to_fan(rays3d):
    """
    Project 3-D trajectory_3d (list of [x,y,z,...] states) onto the x-z
    plane so that plot_ray_fan (which reads r['trajectory']) can display it.
    """
    out = []
    for r in rays3d:
        r2 = dict(r)
        r2['trajectory'] = [(float(s[0]), float(s[2])) for s in r['trajectory_3d']]
        out.append(r2)
    return out


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


def _tag_clutter_radar_modes(modes):
    """Mark monostatic surface-clutter radar observables explicitly."""
    for m in modes:
        tau_oneway = float(m.get('tau_oneway_ms', m.get('tau_ms', 0.0)))
        m['tau_oneway_ms'] = tau_oneway
        m['tau_2way_ms'] = float(m.get('tau_2way_ms', 2.0 * tau_oneway))
        m['range_km'] = float(m.get('range_km', m.get('group_path_km', 0.0)))
        m['radar_observable'] = 'two_way_surface_clutter'
        m['radar_equation'] = 'monostatic_sigma0_area_R4'
        m['physics_basis'] = 'Coleman_1997_1998_flux_tube_Barrick_1972_Bragg'
    return modes


def _build_pd(modes, method='gaussian', delay_key='tau_ms', scatter_pd=False):
    """Build P-D spectrum; radar callers pass delay_key='tau_2way_ms'."""
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
    xlabel = 'Two-way Group Delay  (ms)' if delay_key == 'tau_2way_ms' else 'One-way Group Delay  (ms)'
    info_text = _make_pd_info(modes, delay_key=delay_key)
    fig, _ = plot_pd_spectrum(
        tau_ax, pd, modes, title=title, delay_key=delay_key, xlabel=xlabel,
        info_text=info_text)
    _save(fig, filename)


def _add_power_3d(modes, n3d=None, radar_mode=False, sigma_rcs=1.0):
    """Fill Pr_W, Pr_dBW, delta_tau_ms in 3-D P2P mode dicts."""
    from utils import (free_space_loss_dB, to_dBW, radar_equation_W,
                       d_layer_absorption_dB)
    for m in modes:
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
            from models.absorption_collision import complex_absorption_path_dB
            path = m.get('trajectory_3d', m.get('points', []))
            sampler = None if n3d is None else (lambda x, y, z: n3d.Ne_at(x, y, z))
            abs_diag = complex_absorption_path_dB(
                path, cfg.FREQ_MHZ, sampler, complex_cfg,
                two_way=bool(radar_mode),
                wave_mode=m.get('wave_mode', 'O'), geomag=cfg.GEOMAG)
            abs_oneway_dB = float(abs_diag.get('abs_complex_dB', 0.0))
            abs_total_dB = float(abs_diag.get('abs_total_dB', abs_oneway_dB))
            abs_model = 'appleton_hartree_collision_path'
            m['Pr_W'] *= 10.0**(-abs_total_dB / 10.0)
        elif cfg.ABSORPTION.get('enable', False):
            beta_deg = float(m.get('beta_deg_tx',
                             m.get('elevation_tx_deg',
                             m.get('beta_deg', 30.0))))
            fH_L = (float(cfg.GEOMAG.get('fH_MHz', 1.197)) *
                    float(np.cos(np.radians(cfg.GEOMAG.get('dip_deg', 48.7)))))
            ab_kw = {k: v for k, v in cfg.ABSORPTION.items() if k != 'enable'}
            ab_kw['fH_L_MHz'] = fH_L
            abs_oneway_dB = d_layer_absorption_dB(
                cfg.FREQ_MHZ, beta_deg, **ab_kw)
            n_cross = 2 if radar_mode else 1
            abs_total_dB = abs_oneway_dB * n_cross
            abs_model = 'legacy_d_layer'
            abs_diag = {'abs_status': 'ok'}
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

def _add_ray_pe_3d(modes, Ne_2d, irregularity_params):
    pe_cfg = cfg.UPDATE2.get('PE_COUPLING', {})
    if not pe_cfg.get('enable', False):
        return modes
    from models.pe_propagator import PEPropagator
    from models.hybrid_ray_pe import apply_ray_pe_coupling_3d
    pe = PEPropagator(freq_MHz=cfg.FREQ_MHZ)
    for m in modes:
        Pr_in = float(m.get('Pr_W', 0.0))
        Pr_out, diag = apply_ray_pe_coupling_3d(
            pe, m, Ne_2d, cfg.BG_X, cfg.BG_Z, Pr_in,
            pe.k0, irregularity_params, pe_cfg)
        m.update(diag)
        if bool(diag.get('pe_apply_to_power', False)):
            m['Pr_W'] = float(Pr_out)
            from utils import to_dBW
            m['Pr_dBW'] = to_dBW(m['Pr_W'])
        m['delta_tau_ms'] = float(m.get('delta_tau_ms', 0.0)) + float(diag.get('pe_delta_tau_ms', 0.0))
    return modes


def _add_ray_pe_reference_3d(modes, Ne_2d, irregularity_params):
    """Apply local 2.5-D Ray-PE/SSF coupling as a reference-physics output."""
    if not cfg.UPDATE2.get('PE_COUPLING', {}).get('enable', False):
        for m in modes:
            m['pe_solver_dimensionality'] = 'off'
            m['pe_validity'] = 'disabled_not_full_3d_pe'
            m['pe_apply_to_power'] = False
        return modes
    with _temporary_update2({
        'PE_COUPLING': {
            'enable': True,
            'apply_to_power': False,
            'mode': 'local_2p5d',
            'q_boundary': 1.0,
        }
    }):
        return _add_ray_pe_3d(modes, Ne_2d, irregularity_params)


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
    Find grazing angle [deg] at the first upward crossing of h_Es_km in a
    3-D trajectory (physical height = z_cart - (x^2+y^2)/(2*R_E)).
    Returns None if no crossing is found.
    """
    traj = ray3d['trajectory_3d']
    for i in range(1, len(traj)):
        s0, s1 = traj[i-1], traj[i]
        z0p = float(s0[2]) - (float(s0[0])**2 + float(s0[1])**2) / (2.0*cfg.RE_KM)
        z1p = float(s1[2]) - (float(s1[0])**2 + float(s1[1])**2) / (2.0*cfg.RE_KM)
        if z0p < h_Es_km <= z1p:
            t   = (h_Es_km - z0p) / max(z1p - z0p, 1e-6)
            kx  = float(s0[3]) + t*(float(s1[3])-float(s0[3]))
            ky  = float(s0[4]) + t*(float(s1[4])-float(s0[4]))
            kz  = float(s0[5]) + t*(float(s1[5])-float(s0[5]))
            km  = max((kx**2 + ky**2 + kz**2)**0.5, 1e-30)
            return float(np.degrees(np.arcsin(np.clip(abs(kz)/km, 0.0, 1.0))))
    return None


# -- Iono / radar / link helpers ------------------------------------------------

def _make_iono_params(enable_tid=False, enable_es=False,
                       enable_bubble=False, enable_spread_f=False):
    """
    Assemble iono_params dict.
    TID always uses multi-direction 3-component TID_3D.
    Bubble uses 2-D Gaussian (BUBBLE) for IonosphereModel background field;
    3-D field-aligned depletion is applied via RefractiveIndex3D directly.
    """
    return {
        'iri_params':      {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':      {**cfg.TID_3D,     'enable': enable_tid},
        'es_params':       {**cfg.ES,         'enable': enable_es},
        'bubble_params':   {**cfg.BUBBLE,     'enable': enable_bubble},
        'spread_f_params': {**cfg.SPREAD_F,   'enable': enable_spread_f},
    }


def _make_radar_params():
    return {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W,
            'Gt': cfg.GT, 'Gr': cfg.GR}


def _make_tube_radar_params():
    return {
        'freq_MHz'        : cfg.FREQ_MHZ,
        'Pt_W'            : cfg.PT_W,
        'Gt'              : cfg.GT,
        'Gr'              : cfg.GR,
        'sigma_rcs_m2'    : cfg.RADAR['sigma_rcs_m2'],
        'bearing_deg'     : cfg.LINK_BEARING_DEG,
        'sigma0_ground_dB': cfg.RADAR['sigma0_ground_dB'],
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
    Accepts 3-D modes (key 'trajectory_3d') and 2-D P2P modes (key 'points').
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    pcm = ax.pcolormesh(x_km, z_km, Ne_2d.T / 1e10,
                        cmap='jet', shading='auto', alpha=0.5)
    plt.colorbar(pcm, ax=ax, label='Ne  [x10^10 m^-3]',
                 fraction=0.046, pad=0.04)
    colors = cm.tab10(np.linspace(0, 1, max(len(modes), 1)))
    for m, c in zip(modes, colors):
        if 'trajectory_3d' in m:
            xs = [float(s[0]) for s in m['trajectory_3d']]
            zs = [float(s[2]) for s in m['trajectory_3d']]
        elif 'points' in m:
            pts = np.asarray(m['points'])
            xs, zs = pts[:, 0].tolist(), pts[:, 1].tolist()
        else:
            continue
        ax.plot(xs, zs, '-', color=c, lw=1.8,
                label='{} ({:.2f} ms)'.format(m.get('label', '?'), m['tau_ms']))
    ax.axvline(cfg.TX_POS[0], color='g', lw=1.2, ls='--', label='TX')
    ax.axvline(cfg.RX_POS[0], color='b', lw=1.2, ls='--', label='RX')
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


def _save_p2p_modes_plot(scenario, modes, Ne_2d, pert_info=None, title=None):
    fig, _ = plot_p2p_modes(
        modes, Ne_2d, cfg.BG_X, cfg.BG_Z,
        title=title or 'P2P Modes - {}'.format(scenario),
        link_info=_make_link_info(), pert_info=pert_info)
    _save(fig, 'p2p_modes_{}.png'.format(scenario))


def _annotate_p2p_causes(modes, environment=None):
    from models.point_to_point_3d import annotate_mode_causes_3d
    annotate_mode_causes_3d(modes, cfg.P2P_3D, environment=environment)


# -- Scenario 0: IRI baseline ---------------------------------------------------

def run_baseline():
    """IRI background only, 3-D ray tracer.
    Saves tx_fan_baseline.png, pd_baseline.png."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    print("=" * 55)
    print("  Scenario 0: Baseline (IRI only, 3-D RT)")
    print("=" * 55)

    iono  = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    print("  Ne field: shape {}  max = {:.2e} m^-3".format(
        Ne_2d.shape, float(Ne_2d.max())))

    n3d    = _make_n3d(Ne_2d)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI baseline (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info())
    _save(fig, 'tx_fan_baseline.png')

    modes = _find_all_rays_reference_3d(n3d)
    _add_power_3d(modes, n3d=n3d)
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'baseline')
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
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    print("  Scenario 1: IRI + {}-dir TID (3-D RT)".format(n_comp))
    print("=" * 55)
    print("  az_list={} deg  amp_list={}".format(
        cfg.TID_3D['az_deg_list'], cfg.TID_3D['amplitude_list']))

    tid_p = {**cfg.TID_3D, 'enable': True}
    iono  = IonosphereModel(tid_params=tid_p)
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)

    # TID baked into Ne_2d; grad_n2() fast path active (no runtime TID interp)
    n3d    = _make_n3d(Ne_2d)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + TID ({}-dir, 3-D)'.format(n_comp),
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_tid=True))
    _save(fig, 'tx_fan_tid.png')

    # Mini-fan with 3-D TID active to show lateral y-deflection (cross-path TID)
    _rt_mini  = {**cfg.RT_3D, 'n_fan': 7, 'max_steps': 3000}
    n3d_tid   = _make_n3d(Ne_2d, tid_params=tid_p)
    rays_mini = _shoot_rays_fan_reference_3d(n3d_tid, _rt_mini)
    betas = [float(r['beta_deg_tx']) for r in rays_mini]
    ys    = [float(r['y_land_km'])   for r in rays_mini]
    fig3, ax3 = plt.subplots(figsize=(8, 4))
    ax3.plot(betas, ys, 'b.-', lw=1.2)
    ax3.axhline(0.0, color='k', lw=0.6, ls=':')
    ax3.set_xlabel('TX Elevation (deg)')
    ax3.set_ylabel('Lateral landing y (km)')
    ax3.set_title('3-D TID lateral deflection ({}-component)'.format(n_comp))
    ax3.grid(True, alpha=0.3)
    plt.tight_layout()
    _save(fig3, 'ray_deflect_tid.png')

    modes = _find_all_rays_reference_3d(n3d)
    _add_power_3d(modes, n3d=n3d)
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'tid')
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
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d
    from models.es_model import EsLayerModel
    from utils import to_dBW

    print("=" * 55)
    print("  Scenario 2: IRI + Es (3-D RT + Hao 2017)")
    print("=" * 55)
    print("  foEs = {} MHz  f = {} MHz  foEs/f = {:.2f}".format(
        cfg.ES['foEs_MHz'], cfg.FREQ_MHZ,
        cfg.ES['foEs_MHz'] / cfg.FREQ_MHZ))

    iono  = IonosphereModel(es_params={**cfg.ES, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + Es (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_es=True))
    _save(fig, 'tx_fan_es.png')

    modes = _find_all_rays_reference_3d(n3d)
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
    _annotate_p2p_causes(modes, 'es')
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
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    print("=" * 55)
    print("  Scenario 3: IRI + 3-D Field-Aligned Bubble")
    print("=" * 55)
    b = cfg.BUBBLE_3D
    print("  x0={} km  z0={} km  s_perp={} km  s_para={} km  dmax={}".format(
        b['x0_km'], b['z0_km'], b['sigma_perp_km'], b['sigma_para_km'], b['delta_max']))

    # 2-D Gaussian bubble in IonosphereModel for Ne background visualisation;
    # 3-D field-aligned depletion applied by RefractiveIndex3D during ray tracing.
    iono   = IonosphereModel(bubble_params={**cfg.BUBBLE, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)

    bub3d_p = {**cfg.BUBBLE_3D, 'enable': True}
    n3d     = _make_n3d(Ne_2d, bubble_3d_params=bub3d_p)
    rays3d  = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + 3-D Field-Aligned Bubble',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_bubble=True))
    _save(fig, 'tx_fan_bubble.png')

    modes = _find_all_rays_reference_3d(n3d)
    _add_power_3d(modes, n3d=n3d)
    _add_ray_pe_reference_3d(modes, Ne_2d, bub3d_p)
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'bubble')
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
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d
    from models.es_model import EsLayerModel
    from utils import to_dBW

    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    print("  Scenario 4: Full 3-D ({}-dir TID + Es + Bubble)".format(n_comp))
    print("=" * 55)

    tid_p   = {**cfg.TID_3D,   'enable': True}
    bub3d_p = {**cfg.BUBBLE_3D,'enable': True}

    iono = IonosphereModel(
        tid_params=tid_p,
        es_params={**cfg.ES,    'enable': True},
        bubble_params={**cfg.BUBBLE, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)

    n3d          = _make_n3d(Ne_2d, bubble_3d_params=bub3d_p)
    _rays_seed   = _shoot_rays_fan_reference_3d(n3d)
    fig_fan, _ = plot_tx_ray_fan(_rays3d_to_fan(_rays_seed), Ne_2d, cfg.BG_X, cfg.BG_Z,
                                 title='TX Ray Fan - Full 3-D',
                                 link_info=_make_link_info(),
                                 pert_info=_make_pert_info(enable_tid=True, enable_es=True, enable_bubble=True))
    _save(fig_fan, 'tx_fan_full.png')
    modes  = _find_all_rays_reference_3d(n3d)
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

    _add_ray_pe_reference_3d(modes, Ne_2d, bub3d_p)
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'full')
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


# -- Scenario OX: O/X splitting ------------------------------------------------

def _apply_OX_split_3d(modes_iso, Ne_2d, n3d):
    """
    Post-process 3-D isotropic modes into O/X pairs via Appleton-Hartree.

    For each mode the 3-D trajectory is used as-is (isotropic path approximation).
    O and X group paths are computed by integrating ds/n_{O,X} along that path,
    where n_{O,X}(x,z) is the AH phase index from RefractiveIndexAH.n_batch().

    Group path formula (Budden 1961, plasma approximation):
        P'_{O,X} = integral(ds / n_{O,X})   ->  tau_{O,X} = P'_{O,X} / c

    Limitation: O and X rays physically follow different paths; proper 3-D AH
    ray tracing (separate Haselgrove trajectories per mode) is a future upgrade.
    """
    from models.ray_tracer import RefractiveIndexAH
    from models.point_to_point_3d import classify_mode_3d
    from utils import group_delay_ms

    nm_O = RefractiveIndexAH(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                              wave_mode='O', geomag=cfg.GEOMAG)
    nm_X = RefractiveIndexAH(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                              wave_mode='X', geomag=cfg.GEOMAG)

    result = []
    for m in modes_iso:
        traj = m['trajectory_3d']
        n_seg = len(traj) - 1
        if n_seg < 1:
            continue

        # Collect midpoint (x, z_phys) for every segment — one vectorized call each
        xs   = np.empty(n_seg)
        zps  = np.empty(n_seg)
        segs = np.empty(n_seg)

        for i in range(n_seg):
            s0, s1 = traj[i], traj[i + 1]
            x_m  = (float(s0[0]) + float(s1[0])) * 0.5
            y_m  = (float(s0[1]) + float(s1[1])) * 0.5
            z_cm = (float(s0[2]) + float(s1[2])) * 0.5
            xs[i]   = x_m
            zps[i]  = max(z_cm - (x_m**2 + y_m**2) / (2.0 * cfg.RE_KM),
                           float(cfg.BG_Z[0]))
            dx = float(s1[0]) - float(s0[0])
            dy = float(s1[1]) - float(s0[1])
            dz = float(s1[2]) - float(s0[2])
            segs[i] = (dx*dx + dy*dy + dz*dz) ** 0.5

        pts = np.column_stack([xs, zps])          # (n_seg, 2)
        n_O_arr = nm_O.n_batch(pts)               # AH phase index (not n^2)
        n_X_arr = nm_X.n_batch(pts)

        n_O_arr = np.maximum(n_O_arr, 1e-4)
        n_X_arr = np.maximum(n_X_arr, 1e-4)

        # Group paths
        gp_O = float(np.sum(segs / n_O_arr))
        gp_X = float(np.sum(segs / n_X_arr))
        gp_iso = float(m.get('group_path_km',
                             np.sum(segs / np.maximum(n_O_arr, 1e-4))))

        # O-mode copy
        m_O = dict(m)
        m_O['wave_mode']     = 'O'
        m_O['group_path_km'] = gp_O
        m_O['tau_ms']        = group_delay_ms(gp_O)
        m_O['OX_gp_ratio']   = gp_O / max(gp_iso, 1.0)
        m_O['label']         = classify_mode_3d(m_O)
        result.append(m_O)

        # X-mode copy — skip only if X is evanescent over most of path
        x_valid = float(np.mean(n_X_arr > 1e-3))
        if x_valid > 0.5:
            m_X = dict(m)
            m_X['wave_mode']     = 'X'
            m_X['group_path_km'] = gp_X
            m_X['tau_ms']        = group_delay_ms(gp_X)
            m_X['OX_gp_ratio']   = gp_X / max(gp_iso, 1.0)
            m_X['label']         = classify_mode_3d(m_X)
            result.append(m_X)

    result.sort(key=lambda r: r['tau_ms'])
    return result


def run_with_OX():
    """
    IRI + TID + Appleton-Hartree O/X splitting, full 3-D RT.
    Fan and P2P: 3-D Haselgrove (isotropic n^2).
    O/X correction: post-processing via RefractiveIndexAH group-path integration
    along 3-D isotropic trajectories.
    Saves tx_fan_OX.png, pd_OX.png, modes_OX.csv.
    """
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    use_ah3d = bool(cfg.UPDATE2.get('AH3D_ENABLE', False)) and cfg.UPDATE2.get('OX_MODE', 'postproc') == 'ah3d'
    ox_basis = 'AH3D' if use_ah3d else 'AH post-proc diagnostic'
    print("  Scenario OX: {}-dir TID + O/X ({})".format(n_comp, ox_basis))
    print("=" * 55)
    print("  fH = {} MHz   dip = {} deg   decl = {} deg".format(
        cfg.GEOMAG['fH_MHz'], cfg.GEOMAG['dip_deg'], cfg.GEOMAG['decl_deg']))

    tid_p = {**cfg.TID_3D, 'enable': True}
    iono  = IonosphereModel(tid_params=tid_p)
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))
    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - O/X {} (3-D + {}-dir TID)'.format(ox_basis, n_comp),
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_tid=True))
    _save(fig, 'tx_fan_OX.png')

    if use_ah3d:
        from models.point_to_point_ah_3d import find_all_rays_ah_3d
        modes = find_all_rays_ah_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                                    wave_mode='both', geomag=cfg.GEOMAG)
    else:
        # 3-D P2P for isotropic modes (same pipeline as all other 3-D scenarios)
        modes_iso = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D)

        # Appleton-Hartree O/X post-processing on 3-D trajectories
        modes = _apply_OX_split_3d(modes_iso, Ne_2d, n3d)
    for m in modes:
        m['ox_solver'] = 'ah3d_hamiltonian' if use_ah3d else 'postproc_diagnostic'
        m['physics_basis'] = (
            '3d_appleton_hartree_hamiltonian'
            if use_ah3d else
            'diagnostic_only_isotropic_path_appleton_hartree_group_integral')
    _add_power_3d(modes, n3d=n3d)
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'OX')

    print("  P2P modes: {} ({})".format(len(modes), ox_basis))
    print_mode_table(modes)
    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = _save_modes_csv('OX_ah3d' if use_ah3d else 'OX_diagnostic', modes, tau_ax, pd)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _save_p2p_modes_plot('OX_ah3d' if use_ah3d else 'OX_diagnostic', modes, Ne_2d, _make_pert_info(enable_tid=True))
    _plot_pd(tau_ax, pd, modes,
             title='P-D - O/X {} ({}-dir TID)'.format(ox_basis, n_comp),
             filename='pd_OX_ah3d.png' if use_ah3d else 'pd_OX_diagnostic.png')


# -- Scenario R1/R2: OTH radar observables ------------------------------------

def run_radar_point_target():
    """Monostatic point-target radar, IRI background, 3-D RT.

    This is not surface clutter. It uses the standard monostatic point-target
    R^-4 radar equation (Skolnik, Radar Handbook, 3rd ed., 2008).
    Saves tx_fan_radar_point.png, pd_radar_point.png.
    """
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    print("=" * 55)
    print("  Scenario R1: Monostatic point-target radar (3-D RT)")
    print("=" * 55)
    print("  sigma_rcs = {} m^2   bearing = {} deg".format(
        cfg.RADAR['sigma_rcs_m2'], cfg.LINK_BEARING_DEG))

    iono   = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)
    rays3d = _shoot_rays_fan_reference_3d(n3d)
    print("  Ray fan : {} rays".format(len(rays3d)))
    fig, _ = plot_tx_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - Monostatic Point-Target Radar (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info())
    _save(fig, 'tx_fan_radar_point.png')

    modes = _find_all_rays_reference_3d(n3d, rays_fan=rays3d)
    _add_power_3d(modes, n3d=n3d, radar_mode=True,
                  sigma_rcs=float(cfg.RADAR['sigma_rcs_m2']))
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'radar_point_target')
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
    """Single-station OTH surface clutter using 3-D flux tubes.

    The propagation/tube geometry follows Coleman-style OTH flux-tube
    backscatter. Surface return uses sigma0 times illuminated tube area, with
    Bragg/Doppler diagnostics following the HF sea-clutter convention of
    Barrick (1972). No multi-station assumption is made here.
    Saves pd_radar_clutter.png, modes_radar_clutter.csv.
    """
    from models.ionosphere_model import IonosphereModel
    from models.tube_tracer_3d import TubeRayTracer3D
    from models.point_to_point_3d import classify_mode_3d
    from models.radar_clutter import apply_clutter_to_mode, build_surface_clutter_pd_spectrum

    print("=" * 55)
    print("  Scenario R2: Single-station OTH surface clutter")
    print("=" * 55)
    print("  target_range={} km  surface={}  sigma0={} dB".format(
        cfg.RADAR['target_range_km'],
        cfg.UPDATE2.get('RADAR_CLUTTER', {}).get('surface_type', 'sea'),
        cfg.RADAR['sigma0_ground_dB']))

    iono = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d = _make_n3d(Ne_2d)

    tracer = TubeRayTracer3D(n3d, cfg.FREQ_MHZ,
                              tube_params=cfg.TUBE_3D, rt_params=cfg.RT_3D)
    # Initial sigma0 only seeds the tube power; the explicit clutter model below
    # overwrites Pr_W from sigma0, grazing angle, and illuminated area.
    modes = tracer.compute(
        cfg.TX_POS, float(cfg.RADAR['target_range_km']),
        Pt_W=cfg.PT_W, Gt=cfg.GT, Gr=cfg.GR, sigma0=1.0)

    clutter_gate = cfg.UPDATE2.get('RADAR_CLUTTER', {})
    clutter_params = {
        'surface_type': clutter_gate.get('surface_type', 'sea'),
        'sigma0_dB': (clutter_gate.get('sigma0_dB')
                      if clutter_gate.get('sigma0_dB') is not None
                      else cfg.RADAR.get('sigma0_ground_dB', -20.0)),
        'radial_velocity_ms': clutter_gate.get('radial_velocity_ms', 0.0),
        'grazing_angle_deg': clutter_gate.get('grazing_angle_deg'),
        'clutter_class': clutter_gate.get('clutter_class', 'surface_bragg'),
    }

    updated_modes = []
    for m in modes:
        if 'label' not in m:
            m['label'] = classify_mode_3d(m)
        m['delta_tau_ms'] = float(m.get('tau_spread_ms', 0.0))
        updated = apply_clutter_to_mode(
            m, freq_MHz=cfg.FREQ_MHZ, params=clutter_params,
            Pt_W=cfg.PT_W, Gt=cfg.GT, Gr=cfg.GR, apply_to_power=True)
        updated_modes.append(updated)
    modes = _tag_clutter_radar_modes(updated_modes)

    print("  Clutter tube modes: {}".format(len(modes)))
    print_mode_table(modes)
    if modes:
        best = modes[0]
        print("  Best: tau2={:.3f}ms  Doppler={:.3f}Hz  sigma_eff={:.3e}m2".format(
            best.get('tau_2way_ms', 0.0), best.get('doppler_hz', 0.0),
            best.get('sigma_eff_m2', 0.0)))

    tau_ax, pd, clutter_pd_meta = build_surface_clutter_pd_spectrum(
        modes, freq_MHz=cfg.FREQ_MHZ, params=clutter_params)
    for m in modes:
        m.update(clutter_pd_meta)
    main, _, ms = _identify_main(tau_ax, pd, modes, delay_key='tau_2way_ms')
    _print_main(main)
    if ms:
        _print_summary(ms)

    csv_path = _save_modes_csv('radar_clutter', modes, tau_ax, pd,
                               pd_delay_key='tau_2way_ms')
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _plot_pd(
        tau_ax, pd, modes,
        title='P-D - Single-Station OTH Surface Clutter (Bragg range-Doppler projection)',
        filename='pd_radar_clutter.png', delay_key='tau_2way_ms')


def run_radar_baseline():
    """Backward-compatible alias for the point-target radar scenario."""
    return run_radar_point_target()


# -- Scenario SF: Spread-F ------------------------------------------------------

def run_with_spreadf():
    """IRI + spread-F only (Rino 1979 Ne perturbation), 3-D RT.
    Saves ray_fan_spreadf.png (mode paths), pd_spreadf.png."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    print("=" * 55)
    print("  Scenario SF: Spread-F only (3-D RT)")
    print("=" * 55)
    sf = cfg.SPREAD_F
    print("  Cs={} p={} h_screen={} km L0={} km".format(
        sf['Cs'], sf['p'], sf['h_screen_km'], sf['L0_km']))

    iono  = IonosphereModel(spread_f_params={**cfg.SPREAD_F, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d          = _make_n3d(Ne_2d)
    _rays_seed   = _shoot_rays_fan_reference_3d(n3d)
    fig_fan, _ = plot_tx_ray_fan(_rays3d_to_fan(_rays_seed), Ne_2d, cfg.BG_X, cfg.BG_Z,
                                 title='TX Ray Fan - Spread-F',
                                 link_info=_make_link_info(),
                                 pert_info=_make_spreadf_pert_info())
    _save(fig_fan, 'tx_fan_spreadf.png')
    modes  = _find_all_rays_reference_3d(n3d)
    _add_power_3d(modes, n3d=n3d)
    _add_ray_pe_reference_3d(modes, Ne_2d, {**cfg.BUBBLE_3D, 'enable': True, 'x0_km': cfg.RX_RANGE / 2.0, 'z0_km': cfg.SPREAD_F['h_screen_km'], 'sigma_perp_km': cfg.SPREAD_F.get('L0_km', 50.0), 'sigma_para_km': 300.0})
    _add_scatter_3d(modes, n3d)
    _annotate_p2p_causes(modes, 'spreadf')
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
    """OTH radar 3-D flux-tube tracer (TubeRayTracer3D, Coleman 1997/1998).
    Saves pd_tube_radar.png, modes_tube_radar.csv."""
    from models.ionosphere_model import IonosphereModel
    from models.tube_tracer_3d import TubeRayTracer3D
    from models.point_to_point_3d import classify_mode_3d
    from models.radar_clutter import apply_clutter_to_mode, build_surface_clutter_pd_spectrum

    print("=" * 55)
    print("  Scenario TR: OTH Radar 3-D Tube Tracer")
    print("=" * 55)
    print("  target_range={} km  sigma0={} dBsm".format(
        cfg.RADAR['target_range_km'], cfg.RADAR['sigma0_ground_dB']))

    iono   = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)

    tracer = TubeRayTracer3D(n3d, cfg.FREQ_MHZ,
                              tube_params=cfg.TUBE_3D, rt_params=cfg.RT_3D)
    sigma0 = 10.0**(float(cfg.RADAR['sigma0_ground_dB']) / 10.0)
    modes  = tracer.compute(cfg.TX_POS,
                             float(cfg.RADAR['target_range_km']),
                             Pt_W=cfg.PT_W, Gt=cfg.GT, Gr=cfg.GR, sigma0=sigma0)

    clutter_gate = cfg.UPDATE2.get('RADAR_CLUTTER', {})
    clutter_params = {
        'surface_type': clutter_gate.get('surface_type', 'sea'),
        'sigma0_dB': (clutter_gate.get('sigma0_dB')
                      if clutter_gate.get('sigma0_dB') is not None
                      else cfg.RADAR.get('sigma0_ground_dB', -20.0)),
        'radial_velocity_ms': clutter_gate.get('radial_velocity_ms', 0.0),
        'grazing_angle_deg': clutter_gate.get('grazing_angle_deg'),
        'clutter_class': clutter_gate.get('clutter_class', 'surface_bragg'),
    }

    # Add label and delta_tau_ms for table / CSV compatibility
    updated_modes = []
    for m in modes:
        if 'label' not in m:
            m['label'] = classify_mode_3d(m)
        m['delta_tau_ms'] = float(m.get('tau_spread_ms', 0.0))
        updated_modes.append(apply_clutter_to_mode(
            m, freq_MHz=cfg.FREQ_MHZ, params=clutter_params,
            Pt_W=cfg.PT_W, Gt=cfg.GT, Gr=cfg.GR,
            apply_to_power=True))
    modes = _tag_clutter_radar_modes(updated_modes)

    print("  Tube modes: {}".format(len(modes)))
    print_mode_table(modes)
    if modes:
        best = modes[0]
        print("  Best: tau={:.3f}ms  F_focus_3d={:.3f}  A_tube={:.1f} km^2".format(
            best['tau_ms'], best.get('F_focus_3d', 0.0), best.get('A_tube_km2', 0.0)))

    tau_ax, pd, clutter_pd_meta = build_surface_clutter_pd_spectrum(
        modes, freq_MHz=cfg.FREQ_MHZ, params=clutter_params)
    for m in modes:
        m.update(clutter_pd_meta)
    main, _, ms = _identify_main(tau_ax, pd, modes, delay_key='tau_2way_ms')
    _print_main(main)
    if ms:
        _print_summary(ms)

    csv_path = _save_modes_csv('tube_radar', modes, tau_ax, pd,
                               pd_delay_key='tau_2way_ms')
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    _plot_pd(tau_ax, pd, modes,
             title='P-D - OTH Radar 3-D Tube Clutter (Bragg range-Doppler projection)',
             filename='pd_tube_radar.png', delay_key='tau_2way_ms')


# -- Entry point ---------------------------------------------------------------

if __name__ == '__main__':
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
    run_radar_clutter()
    run_with_spreadf()
    run_tube_radar()




