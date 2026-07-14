"""Diagnostics for unknown 3-D ionospheric irregularities.

This module does not choose propagation algorithms.  It only summarizes the
input Ne field so an external Ne cube can be propagated without a required
TID/Spread-F/Es/bubble label.
"""

from __future__ import annotations

import numpy as np


def _spectral_slope_1d(vals: np.ndarray, dx_km: float) -> float:
    arr = np.asarray(vals, dtype=float)
    arr = arr - np.mean(arr)
    if arr.size < 8 or np.std(arr) <= 0.0:
        return 0.0
    psd = np.abs(np.fft.rfft(arr)) ** 2
    k = np.fft.rfftfreq(arr.size, d=max(float(dx_km), 1e-12))
    mask = (k > 0.0) & (psd > 0.0)
    if np.count_nonzero(mask) < 3:
        return 0.0
    coeff = np.polyfit(np.log(k[mask]), np.log(psd[mask]), 1)
    return float(-coeff[0])


def diagnose_irregularity_field(provider, x_km, y_km, z_km, t: float = 0.0) -> dict:
    """Return gradient, anisotropy, spectrum, and candidate-label diagnostics."""
    x = np.asarray(x_km, dtype=float)
    y = np.asarray(y_km, dtype=float)
    z = np.asarray(z_km, dtype=float)
    cube = provider.sample_grid3d(x, y, z, t=t)
    dx = float(np.mean(np.diff(x))) if x.size > 1 else 1.0
    dy = float(np.mean(np.diff(y))) if y.size > 1 else 1.0
    dz = float(np.mean(np.diff(z))) if z.size > 1 else 1.0
    gx, gy, gz = np.gradient(cube, dx, dy, dz, edge_order=1)
    grad_mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2)
    mid_y = y.size // 2
    mid_z = z.size // 2
    slope_x = _spectral_slope_1d(cube[:, mid_y, mid_z], dx)
    var_x = float(np.var(np.diff(cube, axis=0))) if x.size > 1 else 0.0
    var_y = float(np.var(np.diff(cube, axis=1))) if y.size > 1 else 0.0
    var_z = float(np.var(np.diff(cube, axis=2))) if z.size > 1 else 0.0
    anisotropy = max(var_x, var_y, var_z) / max(min(v for v in (var_x, var_y, var_z) if v >= 0.0), 1e-30)
    meta = provider.get_metadata()
    components = meta.get("components", [])
    candidates = []
    for comp in components:
        if comp.get("enabled", False):
            candidates.append(comp.get("component_type", "unknown"))
    if not candidates:
        rel_std = float(np.std(cube) / max(np.mean(cube), 1e-30))
        if slope_x >= 2.0 and rel_std > 1.0e-3:
            candidates.append("SpreadF_like_power_law")
        if np.max(cube[:, :, 0]) > np.mean(cube) * 1.5:
            candidates.append("Es_like_low_altitude_enhancement")
    return {
        "diagnostic_only": True,
        "Ne_min_m3": float(np.min(cube)),
        "Ne_max_m3": float(np.max(cube)),
        "Ne_mean_m3": float(np.mean(cube)),
        "grad_Ne_rms_m3_per_km": float(np.sqrt(np.mean(grad_mag ** 2))),
        "spectral_slope_x": slope_x,
        "anisotropy_indicator": float(anisotropy),
        "candidate_components": candidates,
        "provider_name": meta.get("name", ""),
        "provider_components": [c.get("component_type", "") for c in components],
    }
