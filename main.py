"""
main.py - Entry point for the HF hybrid propagation model.
All scenarios use 3-D Haselgrove ray tracing (Jones & Stephenson 1975).
TID irregularities: 3-component multi-direction superposition (config.TID_3D).
Plasma bubble: 3-D field-aligned Gaussian depletion (config.BUBBLE_3D).

Usage:
    conda run -n pytorch_cpu python main.py

Scenarios (each saves ray_fan_*.png and pd_*.png to output/):
    run_baseline()        IRI only, 3-D RT
    run_with_tid()        IRI + 3-dir TID (3-D), lateral deflection plot
    run_with_es()         IRI + Es, 3-D RT + Hao 2017 power
    run_with_bubble()     IRI + 3-D field-aligned bubble
    run_full()            3-dir TID + Es + 3-D bubble (mode paths overlay)
    run_with_OX()         3-dir TID, 3-D fan; OX power via 2-D AH (pending 3-D port)
    run_radar_baseline()  OTH radar, 3-D RT
    run_with_spreadf()    3-dir TID + spread-F, 3-D RT
    run_tube_radar()      OTH radar, 3-D flux-tube tracer (Coleman 1997/1998)

Test link: TX @ 30N 120E, RX @ 1169 km, f = 10 MHz.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import config as cfg
from viz.plot_utils import (print_mode_table, plot_ray_fan, plot_pd_spectrum,
                             save_modes_csv, build_ray_info_text)
from models.hybrid_model import build_pd_spectrum, identify_main_mode

_OUT = os.path.join(os.path.dirname(__file__), 'output')
os.makedirs(_OUT, exist_ok=True)


# -- Shared helpers -------------------------------------------------------------

def _save(fig, name):
    """Save figure to output/ and close it."""
    fig.savefig(os.path.join(_OUT, name), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: output/{}".format(name))


# -- 3-D model helpers ----------------------------------------------------------

def _make_n3d(Ne_2d, tid_params=None, bubble_3d_params=None):
    """Build RefractiveIndex3D with optional TID and 3-D field-aligned bubble."""
    from models.ray_tracer_3d import RefractiveIndex3D
    return RefractiveIndex3D(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                              tid_params=tid_params,
                              bubble_3d_params=bubble_3d_params,
                              earth_flat=True)


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


def _add_power_3d(modes, radar_mode=False, sigma_rcs=1.0):
    """Fill Pr_W, Pr_dBW, delta_tau_ms in 3-D P2P mode dicts."""
    from utils import free_space_loss_dB, to_dBW, radar_equation_W
    for m in modes:
        gp = float(m['group_path_km'])
        if radar_mode:
            m['Pr_W'] = radar_equation_W(cfg.PT_W, cfg.GT, cfg.GR,
                                          cfg.FREQ_MHZ, gp, sigma_rcs)
        else:
            L_dB      = free_space_loss_dB(gp, cfg.FREQ_MHZ)
            m['Pr_W'] = cfg.PT_W * cfg.GT * cfg.GR * 10.0**(-L_dB / 10.0)
        m['Pr_dBW']       = to_dBW(m['Pr_W'])
        m['delta_tau_ms'] = 0.0


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


def _print_main(main):
    if main:
        pr_str = ("{:.1f} dBW".format(main['Pr_dBW'])
                  if 'Pr_dBW' in main else "n/a")
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
                        cmap='Blues', shading='auto', alpha=0.6)
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


# -- Scenario 0: IRI baseline ---------------------------------------------------

def run_baseline():
    """IRI background only, 3-D ray tracer.
    Saves ray_fan_baseline.png, pd_baseline.png."""
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
    rays3d = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI baseline (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info())
    _save(fig, 'ray_fan_baseline.png')

    modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                             rays_fan=rays3d)
    _add_power_3d(modes)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('baseline', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D - IRI baseline (3-D)')
    _save(fig2, 'pd_baseline.png')


# -- Scenario 1: IRI + multi-direction TID -------------------------------------

def run_with_tid():
    """IRI + 3-component multi-direction TID, 3-D ray tracer.
    Saves ray_fan_tid.png, pd_tid.png, ray_deflect_tid.png."""
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
    rays3d = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + TID ({}-dir, 3-D)'.format(n_comp),
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_tid=True))
    _save(fig, 'ray_fan_tid.png')

    # Mini-fan with 3-D TID active to show lateral y-deflection (cross-path TID)
    _rt_mini  = {**cfg.RT_3D, 'n_fan': 7, 'max_steps': 3000}
    n3d_tid   = _make_n3d(Ne_2d, tid_params=tid_p)
    rays_mini = shoot_rays_fan_3d(cfg.TX_POS, n3d_tid, cfg.FREQ_MHZ, _rt_mini)
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

    modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                             rays_fan=rays3d)
    _add_power_3d(modes)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('tid', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - IRI + TID ({}-dir, 3-D)'.format(n_comp))
    _save(fig2, 'pd_tid.png')


# -- Scenario 2: IRI + Es -------------------------------------------------------

def run_with_es():
    """IRI + Es layer. 3-D RT + Es three-segment power correction (Hao 2017).
    Saves ray_fan_es.png, pd_es.png."""
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
    rays3d = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + Es (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_es=True))
    _save(fig, 'ray_fan_es.png')

    modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                             rays_fan=rays3d)
    _add_power_3d(modes)

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

    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('es', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes, title='P-D - IRI + Es (3-D)')
    _save(fig2, 'pd_es.png')


# -- Scenario 3: IRI + 3-D field-aligned bubble --------------------------------

def run_with_bubble():
    """IRI + 3-D field-aligned plasma bubble (Ma et al. 2023, Martinis 2023).
    Saves ray_fan_bubble.png, pd_bubble.png."""
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
    rays3d  = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    print("  Ray fan : {} rays".format(len(rays3d)))

    fig, _ = plot_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - IRI + 3-D Field-Aligned Bubble',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_bubble=True))
    _save(fig, 'ray_fan_bubble.png')

    modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                             rays_fan=rays3d)
    _add_power_3d(modes)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('bubble', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - IRI + 3-D Bubble')
    _save(fig2, 'pd_bubble.png')


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
    _rays_seed   = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    modes  = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                               rays_fan=_rays_seed)
    _add_power_3d(modes)

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

    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)
    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('full', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))

    # Mode-path overlay on full Ne background
    fig, _ = _plot_mode_paths(
        modes, Ne_2d, cfg.BG_X, cfg.BG_Z,
        title='3-D Mode Paths - {}-dir TID + Es + Bubble'.format(n_comp),
        link_info=_make_link_info(),
        pert_info=_make_pert_info(enable_tid=True, enable_es=True, enable_bubble=True))
    _save(fig, 'ray_fan_full.png')

    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - Full 3-D ({}-dir TID + Es + Bubble)'.format(n_comp))
    _save(fig2, 'pd_full.png')


# -- Scenario OX: O/X splitting ------------------------------------------------

def run_with_OX():
    """
    IRI + TID with Appleton-Hartree O/X mode splitting.
    Ray fan: 3-D Haselgrove (isotropic n^2).
    P-D / modes: 2-D AH via HybridPropagationModel.
    NOTE: Full AH magneto-ionic n^2 is not yet ported to RefractiveIndex3D.
    Saves ray_fan_OX.png, pd_OX.png.
    """
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.hybrid_model import HybridPropagationModel

    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    print("  Scenario OX: {}-dir TID + O/X (3-D fan / 2-D AH modes)".format(n_comp))
    print("=" * 55)
    print("  fH = {} MHz   dip = {} deg   decl = {} deg".format(
        cfg.GEOMAG['fH_MHz'], cfg.GEOMAG['dip_deg'], cfg.GEOMAG['decl_deg']))

    # 3-D ray fan (isotropic geometry)
    tid_p = {**cfg.TID_3D, 'enable': True}
    iono  = IonosphereModel(tid_params=tid_p)
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)
    rays3d = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    print("  Ray fan : {} rays".format(len(rays3d)))
    fig, _ = plot_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - O/X (3-D isotropic + {}-dir TID)'.format(n_comp),
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info(enable_tid=True))
    _save(fig, 'ray_fan_OX.png')

    # OX mode splitting via 2-D Appleton-Hartree (AH not yet in RefractiveIndex3D)
    model = HybridPropagationModel(
        _make_iono_params(enable_tid=True),
        _make_radar_params(),
        geomag_params={**cfg.GEOMAG, 'enable_OX': True})
    modes, tau_ax, pd, main = model.compute(cfg.TX_POS, cfg.RX_POS)
    print("  P2P modes: {} (expect ~2x isotropic count)".format(len(modes)))
    print_mode_table(modes)
    _print_main(main)
    csv_path = save_modes_csv('OX', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - O/X splitting ({}-dir TID, 2-D AH)'.format(n_comp))
    _save(fig2, 'pd_OX.png')


# -- Scenario R: OTH radar baseline -------------------------------------------

def run_radar_baseline():
    """OTH radar monostatic, IRI background, 3-D RT.
    Saves ray_fan_radar.png, pd_radar.png."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    print("=" * 55)
    print("  Scenario R: OTH Radar baseline (3-D RT)")
    print("=" * 55)
    print("  sigma_rcs = {} m^2   bearing = {} deg".format(
        cfg.RADAR['sigma_rcs_m2'], cfg.LINK_BEARING_DEG))

    iono   = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)
    rays3d = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    print("  Ray fan : {} rays".format(len(rays3d)))
    fig, _ = plot_ray_fan(_rays3d_to_fan(rays3d), Ne_2d, cfg.BG_X, cfg.BG_Z,
                          title='Ray Fan - OTH Radar (3-D)',
                          link_info=_make_link_info(),
                          pert_info=_make_pert_info())
    _save(fig, 'ray_fan_radar.png')

    modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                             rays_fan=rays3d)
    _add_power_3d(modes, radar_mode=True,
                  sigma_rcs=float(cfg.RADAR['sigma_rcs_m2']))
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('radar', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - OTH Radar (3-D)')
    _save(fig2, 'pd_radar.png')


# -- Scenario SF: Spread-F ------------------------------------------------------

def run_with_spreadf():
    """IRI + 3-dir TID + spread-F (Rino 1979 Ne perturbation), 3-D RT.
    Saves ray_fan_spreadf.png (mode paths), pd_spreadf.png."""
    from models.ionosphere_model import IonosphereModel
    from models.ray_tracer_3d import shoot_rays_fan_3d
    from models.point_to_point_3d import find_all_rays_3d

    n_comp = cfg.TID_3D['n_components']
    print("=" * 55)
    print("  Scenario SF: {}-dir TID + Spread-F (3-D RT)".format(n_comp))
    print("=" * 55)
    sf = cfg.SPREAD_F
    print("  Cs={} p={} h_screen={} km L0={} km".format(
        sf['Cs'], sf['p'], sf['h_screen_km'], sf['L0_km']))

    tid_p = {**cfg.TID_3D, 'enable': True}
    iono  = IonosphereModel(tid_params=tid_p,
                             spread_f_params={**cfg.SPREAD_F, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d          = _make_n3d(Ne_2d)
    _rays_seed   = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
    modes  = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                               rays_fan=_rays_seed)
    _add_power_3d(modes)
    print("  P2P modes: {}".format(len(modes)))
    print_mode_table(modes)

    tau_ax, pd = build_pd_spectrum(modes)
    main, _, _ = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    csv_path = save_modes_csv('spreadf', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))

    fig, _ = _plot_mode_paths(modes, Ne_2d, cfg.BG_X, cfg.BG_Z,
                               title='3-D Mode Paths - {}-dir TID + Spread-F'.format(n_comp),
                               link_info=_make_link_info(),
                               pert_info=_make_pert_info(enable_tid=True))
    _save(fig, 'ray_fan_spreadf.png')

    fig2, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - {}-dir TID + Spread-F (3-D)'.format(n_comp))
    _save(fig2, 'pd_spreadf.png')


# -- Scenario TR: 3-D flux-tube OTH radar --------------------------------------

def run_tube_radar():
    """OTH radar 3-D flux-tube tracer (TubeRayTracer3D, Coleman 1997/1998).
    Saves pd_tube_radar.png, modes_tube_radar.csv."""
    from models.ionosphere_model import IonosphereModel
    from models.tube_tracer_3d import TubeRayTracer3D
    from models.point_to_point_3d import classify_mode_3d

    print("=" * 55)
    print("  Scenario TR: OTH Radar 3-D Tube Tracer")
    print("=" * 55)
    print("  target_range={} km  sigma0={} dBsm".format(
        cfg.RADAR['target_range_km'], cfg.RADAR['sigma0_ground_dB']))

    iono   = IonosphereModel()
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n3d    = _make_n3d(Ne_2d)

    sigma0 = 10.0**(float(cfg.RADAR['sigma0_ground_dB']) / 10.0)
    tracer = TubeRayTracer3D(n3d, cfg.FREQ_MHZ,
                              tube_params=cfg.TUBE_3D, rt_params=cfg.RT_3D)
    modes  = tracer.compute(cfg.TX_POS,
                             float(cfg.RADAR['target_range_km']),
                             Pt_W=cfg.PT_W, Gt=cfg.GT, Gr=cfg.GR, sigma0=sigma0)

    # Add label and delta_tau_ms for table / CSV compatibility
    for m in modes:
        if 'label' not in m:
            m['label'] = classify_mode_3d(m)
        m['delta_tau_ms'] = float(m.get('tau_spread_ms', 0.0))

    print("  Tube modes: {}".format(len(modes)))
    print_mode_table(modes)
    if modes:
        best = modes[0]
        print("  Best: tau={:.3f}ms  F_focus_3d={:.3f}  A_tube={:.1f} km^2".format(
            best['tau_ms'], best.get('F_focus_3d', 0.0), best.get('A_tube_km2', 0.0)))

    tau_ax, pd = build_pd_spectrum(modes, method='tube')
    main, _, ms = identify_main_mode(tau_ax, pd, modes)
    _print_main(main)
    if ms:
        _print_summary(ms)

    csv_path = save_modes_csv('tube_radar', modes, tau_ax, pd, cfg.FREQ_MHZ, _OUT)
    print("  Saved: output/{}".format(os.path.basename(csv_path)))
    fig, _ = plot_pd_spectrum(tau_ax, pd, modes,
                               title='P-D - OTH Radar 3-D Tube Tracer')
    _save(fig, 'pd_tube_radar.png')


# -- Entry point ---------------------------------------------------------------

if __name__ == '__main__':
    if not check_env():
        sys.exit(1)

    run_baseline()
    run_with_tid()
    run_with_es()
    run_with_bubble()
    run_full()
    run_with_OX()
    run_radar_baseline()
    run_with_spreadf()
    run_tube_radar()
