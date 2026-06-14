"""
Shared visualization utilities.
Imported by main.py and viz scripts.
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import to_dBW


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


def plot_ray_fan(rays: list,
                 ne_2d: np.ndarray | None = None,
                 x_km: np.ndarray | None = None,
                 z_km: np.ndarray | None = None,
                 tx_pos=(0.0, 0.0), rx_pos=(1169.0, 0.0),
                 title: str = 'Ray Fan',
                 save_path: str | None = None):
    """Overlay ray paths on optional density background."""
    fig, ax = plt.subplots(figsize=(12, 6))
    if ne_2d is not None and x_km is not None and z_km is not None:
        ax.pcolormesh(x_km, z_km, ne_2d.T / 1e10,
                      cmap='Blues', shading='auto', alpha=0.35)
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
    ax.legend(fontsize=8)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig, ax


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
