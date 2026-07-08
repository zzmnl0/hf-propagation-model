"""Appleton-Hartree collisional complex absorption path integral.

The absorption model evaluates the collisional Appleton-Hartree complex
phase index along an already traced ray and integrates intensity attenuation,

    A_dB = 8.686 * integral(2 k0 Im(n) ds).

This couples absorption to electron density, collision frequency, geomagnetic
field, propagation direction, and O/X branch without changing the real-valued
ray trajectory solver.
"""

from __future__ import annotations

import numpy as np

from config import GEOMAG, K_FP, LINK_BEARING_DEG
from utils import freq_to_k0


def _b_hat_from_geomag(geomag: dict | None = None) -> np.ndarray:
    g = GEOMAG if geomag is None else geomag
    dip = np.radians(float(g.get("dip_deg", 48.7)))
    decl = np.radians(float(g.get("decl_deg", -5.5)))
    link = np.radians(float(g.get("link_bearing_deg", LINK_BEARING_DEG)))
    b = np.array([
        np.cos(dip) * np.cos(decl - link),
        np.cos(dip) * np.sin(decl - link),
        -np.sin(dip),
    ], dtype=float)
    return b / max(float(np.linalg.norm(b)), 1e-30)


def collision_frequency_s(freq_alt_km: float, params: dict | None = None) -> float:
    """Return electron-neutral collision frequency [s^-1]."""
    p = params or {}
    model = str(p.get("collision_model", "d_region_exponential"))
    if model == "constant":
        return max(float(p.get("nu0_s", 0.0)), 0.0)

    nu_ref = max(float(p.get("nu_ref_s", 5.0e6)), 0.0)
    h_ref = float(p.get("h_ref_km", 80.0))
    scale = max(float(p.get("scale_height_km", 8.0)), 1e-6)
    h_min = float(p.get("h_min_km", 50.0))
    h_max = float(p.get("h_max_km", 140.0))
    h = float(freq_alt_km)
    if h < h_min or h > h_max:
        return 0.0
    return float(nu_ref * np.exp(-(h - h_ref) / scale))


def complex_refractive_index(freq_MHz: float, Ne_m3: float, nu_s: float) -> complex:
    """Cold unmagnetized collisional plasma index; retained for Y=0 checks."""
    return complex_refractive_index_ah(
        freq_MHz, Ne_m3, nu_s, k_hat=np.array([1.0, 0.0, 0.0]),
        wave_mode="O", geomag={**GEOMAG, "fH_MHz": 0.0})


def complex_refractive_index_ah(
    freq_MHz: float,
    Ne_m3: float,
    nu_s: float,
    k_hat,
    wave_mode: str = "O",
    geomag: dict | None = None,
) -> complex:
    """Collisional Appleton-Hartree complex phase index."""
    f_hz = float(freq_MHz) * 1e6
    X = (K_FP ** 2) * max(float(Ne_m3), 0.0) / (f_hz * f_hz)
    Z = max(float(nu_s), 0.0) / (2.0 * np.pi * f_hz)
    g = GEOMAG if geomag is None else geomag
    Y = float(g.get("fH_MHz", 1.197)) / float(freq_MHz)
    U = complex(1.0, -Z)

    if abs(Y) < 1e-14:
        n2 = 1.0 - X / U
    else:
        kh = np.asarray(k_hat, dtype=float)
        kh = kh / max(float(np.linalg.norm(kh)), 1e-30)
        bh = _b_hat_from_geomag(g)
        cos_alpha = float(np.clip(np.dot(kh, bh), -1.0, 1.0))
        sin_alpha_sq = max(1.0 - cos_alpha * cos_alpha, 0.0)
        YL = Y * cos_alpha
        YT_sq = (Y * Y) * sin_alpha_sq
        UX = U - X
        if abs(UX) < 1e-12:
            UX = complex(1e-12, UX.imag)
        half = YT_sq / (2.0 * UX)
        disc = np.sqrt(YT_sq * YT_sq / (4.0 * UX * UX) + YL * YL + 0j)
        if str(wave_mode).upper() == "O":
            denom = U - half + disc
        else:
            denom = U - half - disc
        if abs(denom) < 1e-12:
            denom = complex(1e-12, denom.imag)
        n2 = 1.0 - X / denom

    n = np.sqrt(n2 + 0j)
    if n.imag < 0.0:
        n = complex(n.real, -n.imag)
    return n


def _path_xyz(path: list | tuple | np.ndarray) -> np.ndarray:
    arr = np.asarray(path, dtype=float)
    if arr.ndim != 2 or arr.shape[0] < 2:
        return np.empty((0, 3), dtype=float)
    if arr.shape[1] >= 3:
        return arr[:, :3]
    out = np.zeros((arr.shape[0], 3), dtype=float)
    out[:, 0] = arr[:, 0]
    out[:, 2] = arr[:, 1]
    return out


def complex_absorption_path_dB(
    path: list | tuple | np.ndarray,
    freq_MHz: float,
    Ne_sampler,
    params: dict | None = None,
    two_way: bool = False,
    wave_mode: str = "O",
    geomag: dict | None = None,
) -> dict:
    """Integrate one-way or two-way AH-collisional absorption along a ray path."""
    p = params or {}
    xyz = _path_xyz(path)
    model = p.get("collision_model", "d_region_exponential")
    if xyz.shape[0] < 2:
        return {"abs_complex_dB": 0.0, "abs_total_dB": 0.0,
                "abs_status": "empty_path", "collision_model": model,
                "abs_index_model": "appleton_hartree_collision"}
    if Ne_sampler is None:
        return {"abs_complex_dB": 0.0, "abs_total_dB": 0.0,
                "abs_status": "no_provider", "collision_model": model,
                "abs_index_model": "appleton_hartree_collision"}

    k0_m = freq_to_k0(freq_MHz) / 1000.0
    one_way = 0.0
    for a, b in zip(xyz[:-1], xyz[1:]):
        seg = b - a
        ds_m = float(np.linalg.norm(seg)) * 1000.0
        if ds_m <= 0.0:
            continue
        mid = 0.5 * (a + b)
        Ne = max(float(Ne_sampler(float(mid[0]), float(mid[1]), float(mid[2]))), 0.0)
        nu = collision_frequency_s(float(mid[2]), p)
        if nu <= 0.0 or Ne <= 0.0:
            continue
        n_c = complex_refractive_index_ah(
            freq_MHz, Ne, nu, seg, wave_mode=wave_mode, geomag=geomag)
        one_way += 8.686 * 2.0 * k0_m * max(float(n_c.imag), 0.0) * ds_m

    one_way = float(max(one_way, 0.0))
    total = one_way * (2.0 if two_way else 1.0)
    return {
        "abs_complex_dB": one_way,
        "abs_total_dB": float(total),
        "abs_status": "ok",
        "collision_model": model,
        "abs_index_model": "appleton_hartree_collision",
    }
