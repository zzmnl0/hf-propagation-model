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
    if ne_line:
        lines.append(ne_line)

    pert = pert_info or {}
    if 'TID' in pert:
        p = pert['TID']
        lines.append("[TID] amp={:.2f}  lambda={:.0f} km  T={:.0f} s".format(
            p.get('amplitude', 0), p.get('lambda_h_km', 0), p.get('T_s', 0)))
    if 'Es' in pert:
        p = pert['Es']
        lines.append("[Es]  foEs={:.1f} MHz  h_Es={:.0f} km".format(
            p.get('foEs_MHz', 0), p.get('h_Es_km', 0)))
    if 'Bubble' in pert:
        p = pert['Bubble']
        lines.append("[Bub] z0={:.0f} km  Lx={:.0f} km  dep={:.2f}".format(
            p.get('z0_km', 0), p.get('Lx_km', 0), p.get('delta_max', 0)))

    return "\n".join(lines)


# ── CSV export ────────────────────────────────────────────────────────────────

def save_modes_csv(scenario: str,
                   mode_results: list,
                   tau_axis: np.ndarray,
                   pd_W: np.ndarray,
                   freq_MHz: float,
                   out_dir: str = 'output') -> str:
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
        'tau_ms', 'tau_2way_ms', 'delta_tau_ms',
        'Pr_W', 'Pr_dBW', 'pd_power_W',
        'h_reflect_km', 'group_path_km',
        'beta_deg', 'phi_deg',
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
            row['pd_power_W']  = float(np.interp(m['tau_ms'], tau_axis, pd_W))
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
                        cmap='viridis', shading='auto')
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
                 save_path: str | None = None):
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
                            cmap='Blues', shading='auto', alpha=0.6)
        plt.colorbar(pcm, ax=ax, label='Ne  [x10^10 m^-3]',
                     fraction=0.046, pad=0.04)

    for r in rays:
        xs = [s[0] for s in r['trajectory']]
        zs = [s[1] for s in r['trajectory']]
        ax.plot(xs, zs, 'r-', lw=0.8, alpha=0.65)

    ax.axvline(tx_pos[0], color='g', lw=1.2, ls='--', label='TX')
    ax.axvline(rx_pos[0], color='b', lw=1.2, ls='--', label='RX')
    ax.set_xlabel('Distance  (km)')
    ax.set_ylabel('Height  (km)')
    ax.set_ylim([0, 500])
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


# ── P-D spectrum plot ─────────────────────────────────────────────────────────

def plot_pd_spectrum(tau_axis: np.ndarray, pd_W: np.ndarray,
                     mode_results: list | None = None,
                     title: str = 'P-D Spectrum',
                     save_path: str | None = None):
    """Power-delay spectrum with optional mode markers."""
    fig, ax = plt.subplots(figsize=(10, 4))
    pd_dBW = 10.0 * np.log10(np.maximum(pd_W, 1e-30))
    ax.plot(tau_axis, pd_dBW, 'b-', lw=1.5)
    if mode_results:
        ymin = pd_dBW.min()
        for m in mode_results:
            ax.axvline(m['tau_ms'], color='r', ls='--', alpha=0.55, lw=0.9)
            ax.text(m['tau_ms'] + 0.03, ymin + 1,
                    m.get('label', ''), rotation=90, fontsize=7, color='r')
    ax.set_xlabel('Group Delay  (ms)')
    ax.set_ylabel('Power  (dBW)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
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
