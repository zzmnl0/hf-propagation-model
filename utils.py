"""
Shared physics utility functions used across all modules.
Units convention (unless otherwise stated):
    distances  -> km
    heights    -> km
    frequency  -> MHz
    time       -> s (or ms where noted)
    power      -> W (or dBW where noted)

Visualization utilities are in viz/plot_utils.py.
"""
import numpy as np
from config import C_MS, C_KMS, K_FP


# ── Electromagnetic helpers ───────────────────────────────────────────────────

def ne_to_n2(Ne: np.ndarray, freq_MHz: float) -> np.ndarray:
    """Electron density [m^-3] -> refractive index squared n^2  (isotropic)."""
    fp2 = K_FP**2 * Ne
    n2  = 1.0 - fp2 / (freq_MHz * 1e6)**2
    return np.maximum(n2, 1e-6)

def n2_to_n(n2: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(n2, 1e-6))

def freq_to_k0(freq_MHz: float) -> float:
    """Free-space wavenumber k0 = 2pi/lambda  [km^-1]."""
    lam_km = C_KMS / (freq_MHz * 1e6)
    return 2.0 * np.pi / lam_km

def wavelength_m(freq_MHz: float) -> float:
    """Free-space wavelength [m]."""
    return C_MS / (freq_MHz * 1e6)

def free_space_loss_dB(D_km: float, freq_MHz: float) -> float:
    """Free-space path loss  L = 20 log10(4*pi*D/lambda)  [dB]."""
    lam_m = wavelength_m(freq_MHz)
    return 20.0 * np.log10(4.0 * np.pi * D_km * 1e3 / lam_m)

def radar_equation_W(Pt_W: float, Gt: float, Gr: float,
                     freq_MHz: float, group_path_km: float,
                     sigma_m2: float) -> float:
    """
    Monostatic radar received power (one-way group path length).

    Pr = (Pt * Gt * Gr * lambda^2 * sigma) / ((4*pi)^3 * R^4)

    group_path_km : one-way group path [km] (ray-traced, not slant range)
    sigma_m2      : target radar cross-section [m^2]
    """
    lam_m = C_MS / (freq_MHz * 1e6)
    R_m   = group_path_km * 1e3
    return (Pt_W * Gt * Gr * lam_m**2 * sigma_m2
            / ((4.0 * np.pi)**3 * R_m**4))

def group_delay_ms(group_path_km: float) -> float:
    """Group path [km] -> group delay [ms]."""
    return group_path_km / C_KMS * 1e3


# ── Geomagnetic helper ────────────────────────────────────────────────────────

def get_geomag(lat: float, lon: float,
               alt_km: float = 300.0,
               dt=None) -> dict:
    """
    Compute IGRF geomagnetic parameters at a given location via ppigrf.

    Returns dict with fH_MHz, dip_deg, decl_deg.
    Falls back to hard-coded mid-latitude defaults if ppigrf is unavailable.
    """
    try:
        import ppigrf
        from datetime import datetime as _dt
        dt = dt or _dt(2020, 1, 1)
        Be, Bn, Bu = ppigrf.igrf(lon, lat, alt_km, dt)
        Be = float(np.ravel(Be)[0])
        Bn = float(np.ravel(Bn)[0])
        Bu = float(np.ravel(Bu)[0])
        F  = np.sqrt(Be**2 + Bn**2 + Bu**2)
        I  = float(np.degrees(np.arctan2(-Bu, np.sqrt(Be**2 + Bn**2))))
        D  = float(np.degrees(np.arctan2(Be, Bn)))
        fH = 2.7994e-5 * F      # fH [MHz] = 2.7994e10 [Hz/T] * F[nT]*1e-9 / 1e6
        return {'fH_MHz': round(fH, 3), 'dip_deg': round(I, 1), 'decl_deg': round(D, 1)}
    except Exception:
        return {'fH_MHz': 1.197, 'dip_deg': 48.7, 'decl_deg': -5.5}


# ── Power unit conversions ────────────────────────────────────────────────────

def to_dBW(P_W: float) -> float:
    return 10.0 * np.log10(max(P_W, 1e-30))

def from_dBW(P_dBW: float) -> float:
    return 10.0 ** (P_dBW / 10.0)


# ── Coordinate helpers ────────────────────────────────────────────────────────

def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance between two geographic points [km]."""
    from config import RE_KM
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    return 2.0 * RE_KM * np.arcsin(np.sqrt(a))

def bearing_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Initial bearing from point 1 to point 2 [deg, clockwise from North]."""
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlon = lon2 - lon1
    x = np.cos(lat2) * np.sin(dlon)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    return (np.degrees(np.arctan2(x, y)) + 360.0) % 360.0
