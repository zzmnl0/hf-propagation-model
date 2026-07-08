"""Scatter-coupled P-D spectrum kernels.

The kernels broaden delay distributions from mode diagnostics without treating
sigma_phi_sq or S4 as direct power multipliers. Total mode power is conserved
for implemented Born/Rytov kernels. MPS_required is left unchanged unless a real
MPS result is available in the mode diagnostics.
"""
from __future__ import annotations

import numpy as np


def _normalized_gaussian(tau_axis: np.ndarray, tau0: float, sigma_ms: float) -> np.ndarray:
    sigma = max(float(sigma_ms), 1e-9)
    k = np.exp(-0.5 * ((tau_axis - float(tau0)) / sigma) ** 2)
    s = float(np.sum(k))
    if s <= 0.0:
        out = np.zeros_like(tau_axis, dtype=float)
        out[int(np.argmin(np.abs(tau_axis - tau0)))] = 1.0
        return out
    return k / s


def make_delay_kernel(mode: dict, tau_axis: np.ndarray, tau_res_ms: float | None = None) -> np.ndarray:
    method = str(mode.get('scatter_method', '') or '').lower()
    if tau_res_ms is None:
        tau_res_ms = float(np.median(np.diff(tau_axis))) if len(tau_axis) > 1 else 0.02
    tau0 = float(mode.get('tau_ms', 0.0))
    base_sigma = max(float(mode.get('delta_tau_ms', 0.0)), float(tau_res_ms))
    if method == 'mps_required':
        return _normalized_gaussian(tau_axis, tau0, base_sigma)
    if method == 'born':
        sigma = max(base_sigma, float(mode.get('sigma_tau_ms', 0.0)), float(tau_res_ms))
        return _normalized_gaussian(tau_axis, tau0, sigma)
    if method == 'rytov':
        sigma = max(base_sigma, float(mode.get('sigma_tau_ms', 0.0)), float(tau_res_ms))
        return _normalized_gaussian(tau_axis, tau0, sigma)
    if method == 'mps' and 'mps_delay_kernel' in mode:
        vals = np.asarray(mode['mps_delay_kernel'], dtype=float)
        if vals.shape == tau_axis.shape and float(np.sum(vals)) > 0.0:
            return vals / float(np.sum(vals))
    return _normalized_gaussian(tau_axis, tau0, base_sigma)


def apply_scatter_kernel(mode_results: list,
                         tau_axis: np.ndarray | None = None,
                         tau_res_ms: float = 0.02) -> tuple[np.ndarray, np.ndarray]:
    if not mode_results:
        return np.array([0.0]), np.array([0.0])
    taus = [float(m['tau_ms']) for m in mode_results]
    spreads = [max(float(m.get('delta_tau_ms', 0.0)), float(m.get('sigma_tau_ms', 0.0)), tau_res_ms)
               for m in mode_results]
    if tau_axis is None:
        tau_min = max(0.0, min(t - 4.0 * s for t, s in zip(taus, spreads)) - 0.5)
        tau_max = max(t + 4.0 * s for t, s in zip(taus, spreads)) + 0.5
        tau_axis = np.arange(tau_min, tau_max, tau_res_ms)
    pd = np.zeros(len(tau_axis), dtype=float)
    for mode in mode_results:
        Pr = float(mode.get('Pr_W', 0.0))
        kernel = make_delay_kernel(mode, tau_axis, tau_res_ms)
        pd += Pr * kernel
    return tau_axis, pd
