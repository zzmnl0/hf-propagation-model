"""
Shared visualization utilities.
Imported by main.py and viz scripts.
"""
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from utils import to_dBW


# ── Annotation helper ─────────────────────────────────────────────────────────

def build_ray_info_text(link_info: dict,
                        ne_2d: np.ndarray | None,
                        z_km: np.ndarray | None,
                        pert_info: dict | None = None) -> str:
    """
    Build the parameter annotation string for ray-fan / mode-path plots.

    Parameters
    ----------
    link_info : dict with keys TX_LAT, TX_LON, freq_MHz, Pt_W, RX_km,
                IRI_DT (datetime), IRI_LAT, IRI_LON
    ne_2d     : (Nx, Nz) electron density [m^-3]; used to compute Ne_max
    z_km      : (Nz,) height array [km]
    pert_info : dict with optional keys 'TID', 'Es', 'Bubble', each a sub-dict
    """
    ne_line = ""
    if ne_2d is not None and z_km is not None:
        ne_col = ne_2d.max(axis=0)            # max over x -> (Nz,)
        iz     = int(np.argmax(ne_col))
        ne_max = float(ne_col[iz])
        z_ne   = float(z_km[iz])
        ne_line = "  Ne_max: {:.2e} m^-3 @ {:.0f} km".format(ne_max, z_ne)

    dt = link_info.get('IRI_DT', '')
    if hasattr(dt, 'strftime'):
        dt = dt.strftime('%Y-%m-%d %H:%M UTC')

    lines = [
        "Test Link",
        "  TX: {:.1f}N {:.1f}E   f = {:.1f} MHz".format(
            link_info.get('TX_LAT', 0.0),
            link_info.get('TX_LON', 0.0),
            link_info.get('freq_MHz', 0.0)),
        "  RX: {:.0f} km   Pt = {:.0f} W".format(
            link_info.get('RX_km', 0.0),
            link_info.get('Pt_W', 0.0)),
        "IRI Background",
        "  {}".format(dt),
        "  Lat: {:.1f}N   Lon: {:.1f}E".format(
            link_info.get('IRI_LAT', 0.0),
            link_info.get('IRI_LON', 0.0)),
    ]
    if link_info.get('coordinate_model'):
        lines.extend([
            "Coordinates",
            "  earth={} R={:.1f} km".format(
                link_info.get('earth_model', 'sphere'),
                float(link_info.get('earth_radius_km', 0.0))),
            "  {}".format(link_info.get('coordinate_model')),
            "  {} | AH3D={} | Ray-PE={} | scatter={}".format(
                link_info.get('oneway_or_twoway', 'one_way'),
                link_info.get('reflection_coupling', 'on'),
                link_info.get('ray_pe_coupling', 'off'),
                link_info.get('scatter_coupling', 'diagnostic_only')),
        ])
    if ne_line:
        lines.append(ne_line)

    pert = pert_info or {}
    if 'TID' in pert:
        p = pert['TID']
        amps = p.get('amplitude_list', [p.get('amplitude', 0)])
        lambdas = p.get('lambda_h_km_list', [p.get('lambda_h_km', 0)])
        periods = p.get('period_s_list', [p.get('T_s', 0)])
        lines.append("[TID] amp={}  lambda={} km  T={} s  z0={:.0f} km".format(
            "/".join("{:.2f}".format(float(v)) for v in amps),
            "/".join("{:.0f}".format(float(v)) for v in lambdas),
            "/".join("{:.0f}".format(float(v)) for v in periods),
            float(p.get('z0_km', 280.0))))
    if 'Es' in pert:
        p = pert['Es']
        lines.append("[Es]  foEs={:.1f} MHz  h_Es={:.0f} km".format(
            p.get('foEs_MHz', 0), p.get('h_Es_km', 0)))
    if 'Bubble' in pert:
        p = pert['Bubble']
        lines.append("[Bub] z0={:.0f} km  sperp={:.0f} km  spar={:.0f} km  dep={:.2f}  zbase={:.0f} km".format(
            p.get('z0_km', 0),
            p.get('sigma_perp_km', p.get('Lx_km', 0)),
            p.get('sigma_para_km', p.get('Lz_km', 0)),
            p.get('delta_max', 0),
            p.get('bottomside_height_km', 200.0)))
    if 'SpreadF' in pert:
        p = pert['SpreadF']
        lines.append("[Spread-F] dNe_rms={:.2f} Cs={:.1e} p={:.2f} h={:.0f} km".format(
            p.get('delta_Ne_rel_rms', 0), p.get('Cs', 0),
            p.get('p', 0), p.get('h_screen_km', 0)))

    return "\n".join(lines)


# ── CSV export ────────────────────────────────────────────────────────────────

def save_modes_csv(scenario: str,
                   mode_results: list,
                   tau_axis: np.ndarray,
                   pd_W: np.ndarray,
                   freq_MHz: float,
                   out_dir: str = 'output',
                   pd_delay_key: str = 'tau_ms') -> str:
    """
    Save per-mode propagation parameters to a CSV file.

    Columns
    -------
    scenario, label, freq_MHz,
    tau_ms, delta_tau_ms,
    Pr_W, Pr_dBW, pd_power_W,
    h_reflect_km, group_path_km, beta_deg, phi_deg

    pd_power_W  P-D spectrum value at each mode's tau_ms (option A):
                pd(tau_i) = Pr_i + sum_{j!=i} Pr_j * gauss_tail_at_tau_i.
                Equals Pr_W when all modes are well-separated in delay.

    Returns the path of the written file.
    """
    os.makedirs(out_dir, exist_ok=True)
    fields = [
        'scenario', 'label', 'freq_MHz',
        'provider_name', 'provider_components', 'solver_chain',
        'earth_model', 'earth_radius_km', 'coordinate_model',
        'coordinate_semantics', 'state_semantics', 'plot_coordinate_semantics',
        'feature_gates', 'pd_kernel', 'observable_type',
        'delay_key', 'oneway_or_twoway', 'ray_pe_dimensionality',
        'scatter_status', 'calibration_status',
        'tau_ms', 'tau_oneway_ms', 'tau_2way_ms', 'pd_delay_key',
        'pd_delay_ms', 'delta_tau_ms',
        'Pr_W', 'Pr_dBW', 'pd_power_W',
        'h_reflect_km', 'group_path_km',
        'range_km', 'radar_observable', 'radar_equation',
        'beta_deg', 'phi_deg',
        'geom_path_km', 'x_land_km', 'y_land_km', 'z_land_km',
        'az_deg_tx', 'az_deg_rx', 'azimuth_deflect_deg', 'az_deflect_deg',
        'abs_dB', 'abs_total_dB', 'abs_model', 'abs_complex_dB',
        'abs_status', 'collision_model', 'abs_index_model',
        'wave_mode', 'cause_layer', 'cause_launch', 'cause_mechanism',
        'cause_summary', 'p2p_seed_source',
        'OX_gp_ratio', 'ox_solver', 'branch_status', 'group_index_model',
        'es_mode', 'theta_Es_deg', 'aspect_Es_deg', 'es_loss_dB',
        'scatter_method', 'sigma_phi_sq', 'sigma_phi_rad',
        'S4', 'sigma_tau_ms', 'MPS_required',
        'scatter_solver', 'scatter_solver_dimensionality', 'scatter_validity',
        'scatter_power_applied', 'Pr_before_scatter', 'Pr_after_scatter',
        'scatter_loss_dB', 'energy_conservation_error',
        'kernel_normalization', 'pd_kernel_source',
        'pe_coupled', 'pe_apply_to_power', 'pe_status', 'pe_frame_status',
        'pe_solver_dimensionality', 'pe_validity',
        'pe_Pr_correction_dB', 'pe_power_ratio', 'pe_delta_tau_ms',
        'pe_mean_aoa_deg', 'pe_angle_spread_deg',
        'pe_x_entry_km', 'pe_y_entry_km', 'pe_z_entry_km', 'pe_entry_q',
        'F_focus', 'F_focus_el', 'F_focus_az', 'F_focus_3d',
        'A_tube_km2', 'A_tube_3d_km2', 'tau_spread_ms',
        'surface_type', 'clutter_class', 'sigma0_linear', 'sigma0_dB',
        'sigma_eff_m2', 'grazing_angle_deg', 'bragg_wavenumber_rad_m',
        'bragg_wavelength_m', 'radial_velocity_ms', 'doppler_hz',
        'clutter_Pr_W', 'clutter_Pr_dBW', 'mode_confidence',
        'clutter_pd_model', 'bragg_doppler_minus_hz',
        'bragg_doppler_plus_hz', 'range_gate_span_km',
        'n_range_bins', 'n_doppler_bins',
        'physics_basis',
    ]
    path = os.path.join(out_dir, 'modes_{}.csv'.format(scenario))
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for m in mode_results:
            row = {k: m.get(k, '') for k in fields}
            row['scenario']    = scenario
            row['freq_MHz']    = freq_MHz
            row['phi_deg']     = float(m.get('phi_deg', 0.0))
            row['Pr_dBW']      = m.get('Pr_dBW',
                                        to_dBW(m.get('Pr_W', 1e-30)))
            pd_delay = float(m.get(pd_delay_key, m.get('tau_ms', 0.0)))
            row['pd_delay_key'] = pd_delay_key
            row['pd_delay_ms'] = pd_delay
            row['pd_power_W'] = (float(np.interp(pd_delay, tau_axis, pd_W))
                                 if len(tau_axis) and len(pd_W) else '')
            if row.get('az_deflect_deg', '') == '' and 'azimuth_deflect_deg' in m:
                row['az_deflect_deg'] = m.get('azimuth_deflect_deg', '')
            tau2 = m.get('tau_2way_ms')
            row['tau_2way_ms'] = '' if tau2 is None else float(tau2)
            writer.writerow(row)
    return path


# ── Density field plot ────────────────────────────────────────────────────────

def plot_ne_field(Ne_2d: np.ndarray, x_km: np.ndarray, z_km: np.ndarray,
                  title: str = 'Electron Density Field',
                  save_path: str | None = None):
    """2-D colour map of electron density."""
    fig, ax = plt.subplots(figsize=(11, 5))
    pcm = ax.pcolormesh(x_km, z_km, Ne_2d.T / 1e10,
                        cmap='jet', shading='auto', alpha=0.5)
    plt.colorbar(pcm, ax=ax, label='Ne  [x10^10 m^-3]')
    ax.set_xlabel('Distance  (km)')
    ax.set_ylabel('Height  (km)')
    ax.set_title(title)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig, ax


# ── Ray fan plot ──────────────────────────────────────────────────────────────

def plot_ray_fan(rays: list,
                 ne_2d: np.ndarray | None = None,
                 x_km: np.ndarray | None = None,
                 z_km: np.ndarray | None = None,
                 tx_pos=(0.0, 0.0), rx_pos=(1169.0, 0.0),
                 title: str = 'Ray Fan',
                 link_info: dict | None = None,
                 pert_info: dict | None = None,
                 save_path: str | None = None,
                 xlim: tuple[float, float] | None = None):
    """
    Overlay ray paths on density background.

    Additions vs. prior version
    ---------------------------
    - Electron-density colorbar (Ne [x10^10 m^-3])
    - Parameter annotation box (upper-right): test-link params, IRI background,
      Ne_max, and any active perturbation params (TID / Es / Bubble)

    New parameters
    --------------
    link_info : dict — see build_ray_info_text for required keys
    pert_info : dict — keys 'TID', 'Es', 'Bubble' for active perturbations
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    if ne_2d is not None and x_km is not None and z_km is not None:
        pcm = ax.pcolormesh(x_km, z_km, ne_2d.T / 1e10,
                            cmap='jet', shading='auto', alpha=0.5)
        plt.colorbar(pcm, ax=ax, label='Ne  [x10^10 m^-3]',
                     fraction=0.046, pad=0.04)

    seen_branches = set()
    styles = {'O': ('crimson', '-'), 'X': ('deepskyblue', '--'),
              'iso': ('crimson', '-')}
    for r in rays:
        if 'trajectory_geocentric_lat_lon_height' in r:
            xs, zs = _mode_xz(r)
        else:
            xs = [s[0] for s in r['trajectory']]
            zs = [s[1] for s in r['trajectory']]
        branch = str(r.get('wave_mode', 'iso'))
        color, linestyle = styles.get(branch, ('crimson', '-'))
        label = '{} ray fan'.format(branch) if branch not in seen_branches else None
        ax.plot(xs, zs, color=color, ls=linestyle, lw=0.8,
                alpha=0.65, label=label)
        seen_branches.add(branch)

    ax.axvline(tx_pos[0], color='g', lw=1.2, ls='--', label='TX')
    ax.axvline(rx_pos[0], color='b', lw=1.2, ls='--', label='RX')
    ax.set_xlabel('Distance  (km)')
    ax.set_ylabel('Height  (km)')
    ax.set_ylim([0, 500])
    if xlim is not None:
        ax.set_xlim([float(xlim[0]), float(xlim[1])])
    ax.set_title(title)
    ax.legend(fontsize=8, loc='upper left')

    if link_info is not None:
        info_text = build_ray_info_text(link_info, ne_2d, z_km, pert_info)
        ax.text(0.99, 0.97, info_text,
                transform=ax.transAxes,
                fontsize=6.5, va='top', ha='right',
                family='monospace',
                bbox=dict(boxstyle='round,pad=0.4',
                          facecolor='white', alpha=0.88,
                          edgecolor='gray', lw=0.8))

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig, ax


def plot_tx_ray_fan(rays: list,
                    ne_2d: np.ndarray | None = None,
                    x_km: np.ndarray | None = None,
                    z_km: np.ndarray | None = None,
                    tx_pos=(0.0, 0.0), rx_pos=(1169.0, 0.0),
                    title: str = 'TX Ray Fan',
                    link_info: dict | None = None,
                    pert_info: dict | None = None,
                    save_path: str | None = None,
                    xlim: tuple[float, float] | None = None):
    """Plot only the TX-launched ray fan; it is not a P2P mode plot."""
    return plot_ray_fan(rays, ne_2d, x_km, z_km, tx_pos, rx_pos,
                        title, link_info, pert_info, save_path, xlim=xlim)


def _mode_xz(mode: dict):
    if 'trajectory_geocentric_lat_lon_height' in mode:
        points = np.asarray(
            mode['trajectory_geocentric_lat_lon_height'], dtype=float)
        if points.ndim == 2 and points.shape[0] >= 1 and points.shape[1] == 3:
            lat = np.radians(points[:, 0])
            lon = np.radians(points[:, 1])
            unit = np.column_stack((
                np.cos(lat) * np.cos(lon),
                np.cos(lat) * np.sin(lon),
                np.sin(lat)))
            angles = np.arccos(np.clip(
                np.sum(unit[1:] * unit[:-1], axis=1), -1.0, 1.0))
            ground = np.r_[0.0, np.cumsum(
                float(mode['earth_radius_km']) * angles)]
            return ground.tolist(), points[:, 2].tolist()
    if 'points' in mode:
        pts = np.asarray(mode['points'], dtype=float)
        if pts.ndim == 2 and pts.shape[1] >= 2:
            return pts[:, 0].tolist(), pts[:, 1].tolist()
    if 'trajectory' in mode:
        pts = np.asarray(mode['trajectory'], dtype=float)
        if pts.ndim == 2 and pts.shape[1] >= 2:
            return pts[:, 0].tolist(), pts[:, 1].tolist()
    return [], []


def plot_p2p_modes(modes: list,
                   ne_2d: np.ndarray | None = None,
                   x_km: np.ndarray | None = None,
                   z_km: np.ndarray | None = None,
                   tx_pos=(0.0, 0.0), rx_pos=(1169.0, 0.0),
                   title: str = 'P2P Modes',
                   link_info: dict | None = None,
                   pert_info: dict | None = None,
                   save_path: str | None = None,
                   xlim: tuple[float, float] | None = None):
    """Plot only TX/RX-hitting P2P modes with mode labels and cause fields."""
    fig, ax = plt.subplots(figsize=(12, 6))
    if ne_2d is not None and x_km is not None and z_km is not None:
        pcm = ax.pcolormesh(x_km, z_km, ne_2d.T / 1e10,
                            cmap='jet', shading='auto', alpha=0.5)
        plt.colorbar(pcm, ax=ax, label='Ne  [x10^10 m^-3]',
                     fraction=0.046, pad=0.04)
    styles = {'O': ('crimson', '-'), 'X': ('deepskyblue', '--'),
              'iso': ('purple', '-')}
    for m in modes:
        xs, zs = _mode_xz(m)
        if not xs:
            continue
        branch = str(m.get('wave_mode', 'iso'))
        color, linestyle = styles.get(branch, ('purple', '-'))
        beta = float(m.get('beta_deg_tx', m.get('elevation_tx_deg', m.get('beta_deg', 0.0))))
        landing = float(m.get('p2p_landing_error_km', float('nan')))
        label = '{} tau={:.3f}ms h={:.0f}km beta={:.1f}deg err={:.3f}km'.format(
            m.get('label', '?'), float(m.get('tau_ms', 0.0)),
            float(m.get('h_reflect_km', 0.0)), beta, landing)
        ax.plot(xs, zs, lw=1.8, color=color, ls=linestyle, label=label)
    ax.axvline(tx_pos[0], color='g', lw=1.2, ls='--', label='TX')
    ax.axvline(rx_pos[0], color='b', lw=1.2, ls='--', label='RX')
    ax.set_xlabel('Distance  (km)')
    ax.set_ylabel('Height  (km)')
    ax.set_ylim([0, 500])
    if xlim is not None:
        ax.set_xlim([float(xlim[0]), float(xlim[1])])
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7, loc='upper left', ncol=1)
    if not modes:
        ax.text(
            0.5, 0.48,
            'No converged receiver-hitting root from the local fan seeds.\n'
            'This is not evidence that the scene has no propagation mode.',
            transform=ax.transAxes, ha='center', va='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                      alpha=0.9, edgecolor='crimson', lw=1.0))
    if link_info is not None:
        info_text = build_ray_info_text(link_info, ne_2d, z_km, pert_info)
        ax.text(0.99, 0.97, info_text, transform=ax.transAxes,
                fontsize=6.5, va='top', ha='right', family='monospace',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          alpha=0.88, edgecolor='gray', lw=0.8))
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig, ax


# ── P-D spectrum plot ─────────────────────────────────────────────────────────

def plot_pd_spectrum(tau_axis: np.ndarray, pd_W: np.ndarray,
                     mode_results: list | None = None,
                     title: str = 'P-D Spectrum',
                     save_path: str | None = None,
                     delay_key: str = 'tau_ms',
                     xlabel: str = 'Group Delay  (ms)',
                     info_text: str | None = None):
    """Power-delay spectrum with optional mode markers."""
    fig, ax = plt.subplots(figsize=(10, 4))
    pd_dBW = 10.0 * np.log10(np.maximum(pd_W, 1e-30))
    ax.plot(tau_axis, pd_dBW, 'b-', lw=1.5)
    if mode_results:
        ymin = pd_dBW.min()
        for m in mode_results:
            delay = float(m.get(delay_key, m.get('tau_ms', 0.0)))
            ax.axvline(delay, color='r', ls='--', alpha=0.55, lw=0.9)
            ax.text(delay + 0.03, ymin + 1,
                    m.get('label', ''), rotation=90, fontsize=7, color='r')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Power  (dBW)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if info_text:
        ax.text(0.99, 0.96, info_text, transform=ax.transAxes,
                fontsize=6.8, va='top', ha='right', family='monospace',
                bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                          alpha=0.88, edgecolor='gray', lw=0.8))
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig, ax


# ── Mode table ────────────────────────────────────────────────────────────────

def print_mode_table(mode_results: list) -> None:
    """Pretty-print a table of propagation mode results."""
    header = (f"{'Label':<16} {'tau [ms]':>8} {'dtau [ms]':>9} "
              f"{'Pr [dBW]':>10} {'h_r [km]':>10} {'Path [km]':>11}")
    print(header)
    print('-' * len(header))
    for m in mode_results:
        label  = m.get('label', '?')
        tau    = m.get('tau_ms', float('nan'))
        dtau   = m.get('delta_tau_ms', 0.0)
        Pr_dBW = to_dBW(m.get('Pr_W', 1e-30))
        h_r    = m.get('h_reflect_km', float('nan'))
        path   = m.get('group_path_km', float('nan'))
        print(f"{label:<16} {tau:>8.3f} {dtau:>9.4f} "
              f"{Pr_dBW:>10.1f} {h_r:>10.1f} {path:>11.1f}")


