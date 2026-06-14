"""
Module M2 (basic) – 2-D ray tracer.
Integrates the Haselgrove equations (isotropic, 2-D Cartesian)
with an RK4 stepper and adaptive step size.

Each ray dict returned by trace_single_ray / shoot_rays_fan:
    {
      'beta_deg'       : launch elevation [deg]
      'trajectory'     : list of np.array([x, z, kx, kz]) along path
      'group_path_km'  : accumulated group path [km]
      'tau_ms'         : group delay [ms]
      'h_reflect_km'   : maximum height reached [km]
      'beta_recv_deg'  : arrival elevation at ground [deg]
      'at_Es'          : dict at Es upward crossing, or None
      'at_bubble'      : dict at bubble-bottom upward crossing, or None
      'L_bg_dB'        : background loss (spreading + D-layer abs) [dB]
    }
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from config import C_KMS, FREQ_MHZ, RT, IRI_DT, IRI_LAT
from utils import freq_to_k0, group_delay_ms, free_space_loss_dB


# ── Refractive-index model ────────────────────────────────────────────────────

class RefractiveIndex:
    """
    Interpolable refractive-index model backed by a 2-D Ne grid.

    Parameters
    ----------
    Ne_2d   : (Nx, Nz) electron density [m^-3]
    x_array : (Nx,) horizontal distance [km]
    z_array : (Nz,) height [km]
    freq_MHz: float  operating frequency
    """

    def __init__(self,
                 Ne_2d:   np.ndarray,
                 x_array: np.ndarray,
                 z_array: np.ndarray,
                 freq_MHz: float = FREQ_MHZ):
        self.freq_MHz = freq_MHz
        freq_Hz = freq_MHz * 1e6
        fp2_field = (8.98 ** 2) * Ne_2d          # plasma freq^2 [Hz^2]
        n2_field  = np.maximum(1.0 - fp2_field / freq_Hz**2, 1e-6)

        self._n2_interp = RegularGridInterpolator(
            (x_array, z_array), n2_field,
            method='linear', bounds_error=False, fill_value=1.0
        )
        self._Ne_interp = RegularGridInterpolator(
            (x_array, z_array), Ne_2d,
            method='linear', bounds_error=False, fill_value=0.0
        )

    def n2(self, x: float, z: float) -> float:
        """n^2(x, z)  – interpolated, clamped to [1e-6, inf)."""
        return float(self._n2_interp([[x, z]]))

    def n(self, x: float, z: float) -> float:
        """n(x, z)."""
        return np.sqrt(max(self.n2(x, z), 1e-6))

    def grad_n2(self, x: float, z: float,
                dx: float = 0.5, dz: float = 0.5
                ) -> tuple[float, float]:
        """Central-difference gradient (dn^2/dx, dn^2/dz)  [km^-1]."""
        dn2_dx = (self.n2(x + dx, z) - self.n2(x - dx, z)) / (2.0 * dx)
        dn2_dz = (self.n2(x, z + dz) - self.n2(x, z - dz)) / (2.0 * dz)
        return dn2_dx, dn2_dz


# ── Haselgrove equations (2-D Cartesian, isotropic) ──────────────────────────

def haselgrove_deriv(state: np.ndarray,
                     n_model: RefractiveIndex,
                     k0: float) -> np.ndarray:
    """
    Right-hand side of 2-D Haselgrove equations.

    State : [x, z, kx, kz]  (km, km, km^-1, km^-1)
    k0    : free-space wavenumber [km^-1]
    Returns d(state)/dP  (P = group path [km]).

        dx/dP  = kx / |k|
        dz/dP  = kz / |k|
        dkx/dP = (k0^2 / 2|k|) * dn^2/dx
        dkz/dP = (k0^2 / 2|k|) * dn^2/dz
    """
    x, z, kx, kz = state
    dn2_dx, dn2_dz = n_model.grad_n2(x, z)
    kmag = max(np.hypot(kx, kz), 1e-30)
    half_k02_over_kmag = 0.5 * k0 ** 2 / kmag
    return np.array([
        kx / kmag,
        kz / kmag,
        half_k02_over_kmag * dn2_dx,
        half_k02_over_kmag * dn2_dz,
    ])


def rk4_step(state: np.ndarray,
             n_model: RefractiveIndex,
             k0: float,
             ds: float) -> np.ndarray:
    """Single RK4 step of length ds [km]. Returns new state."""
    k1 = haselgrove_deriv(state,             n_model, k0)
    k2 = haselgrove_deriv(state + ds/2 * k1, n_model, k0)
    k3 = haselgrove_deriv(state + ds/2 * k2, n_model, k0)
    k4 = haselgrove_deriv(state + ds   * k3, n_model, k0)
    return state + (ds / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


# ── Single-ray tracer ─────────────────────────────────────────────────────────

def _martyn_absorption_dB(group_path_km: float,
                          freq_MHz: float) -> float:
    """
    D-layer absorption via Martyn empirical formula.
    L_abs [dB] ~ 677.2 * sec(chi) / (f + fH)^2
    Uses config IRI_DT / IRI_LAT to estimate solar zenith angle chi.
    Returns 0 at nighttime (chi > 90 deg).
    """
    import math
    from datetime import datetime as _dt

    dt  = IRI_DT
    lat = IRI_LAT

    doy = dt.timetuple().tm_yday
    B   = 2.0 * math.pi * (doy - 1) / 365.0
    # solar declination (Spencer formula) [deg]
    delta = math.degrees(
        0.006918
        - 0.399912 * math.cos(B)   + 0.070257 * math.sin(B)
        - 0.006758 * math.cos(2*B) + 0.000907 * math.sin(2*B)
    )
    # local solar hour angle – treat IRI_DT.hour as local solar time
    ha_deg = (dt.hour + dt.minute / 60.0 - 12.0) * 15.0
    lat_r  = math.radians(lat)
    dec_r  = math.radians(delta)
    ha_r   = math.radians(ha_deg)
    cos_chi = (math.sin(lat_r) * math.sin(dec_r)
               + math.cos(lat_r) * math.cos(dec_r) * math.cos(ha_r))
    if cos_chi <= 0.0:          # nighttime: no D-layer
        return 0.0

    fH = 1.4                    # longitudinal gyrofrequency [MHz] (mid-lat)
    return 677.2 / (cos_chi * (freq_MHz + fH) ** 2)


def trace_single_ray(tx_pos:          tuple,
                     beta_deg:         float,
                     n_model:          RefractiveIndex,
                     freq_MHz:         float = FREQ_MHZ,
                     rt_params:        dict  = RT,
                     h_Es_km:          float | None = None,
                     h_bubble_bot_km:  float | None = None
                     ) -> dict:
    """
    Trace one ray from tx_pos at launch elevation beta_deg.

    Adaptive step size: reduced near the reflection layer (n^2 < 0.25)
    to maintain accuracy at the turning point.

    Parameters
    ----------
    tx_pos         : (x0, z0) [km, km]
    beta_deg       : launch elevation above horizon [deg]
    n_model        : RefractiveIndex instance
    freq_MHz       : operating frequency [MHz]
    rt_params      : RT config dict (ds_km, ds_min_km, ds_max_km, ...)
    h_Es_km        : record first upward crossing of this height
    h_bubble_bot_km: record first upward crossing of this height

    Returns ray dict (see module docstring).
    """
    k0       = freq_to_k0(freq_MHz)          # [km^-1]
    ds_nom   = rt_params['ds_km']
    ds_min   = rt_params['ds_min_km']
    ds_max   = rt_params['ds_max_km']
    z_stop   = rt_params['z_stop_km']
    max_steps = rt_params['max_steps']

    beta_rad = np.radians(beta_deg)
    x0, z0   = float(tx_pos[0]), float(tx_pos[1])
    n0       = n_model.n(x0, z0)
    k_mag    = n0 * k0

    state = np.array([x0, z0,
                      k_mag * np.cos(beta_rad),
                      k_mag * np.sin(beta_rad)], dtype=float)

    trajectory  = [state.copy()]
    group_path  = 0.0
    h_reflect   = z0
    at_es       = None
    at_bubble   = None
    prev_z      = z0
    has_risen   = False          # guard: don't terminate until ray has left ground

    for _ in range(max_steps):
        x, z, kx, kz = state

        # ── Adaptive step: shrink near reflection (n^2 -> 0) ─────────────────
        n2_here = n_model.n2(x, z)
        if n2_here < 0.04:
            ds = ds_min
        elif n2_here < 0.20:
            ds = max(ds_min, ds_nom * 0.4)
        else:
            ds = ds_nom

        state      = rk4_step(state, n_model, k0, ds)
        group_path += ds

        x, z, kx, kz = state
        trajectory.append(state.copy())

        if z > 50.0:
            has_risen = True
        h_reflect = max(h_reflect, z)

        # ── Es upward crossing ────────────────────────────────────────────────
        if h_Es_km is not None and prev_z < h_Es_km <= z and at_es is None:
            grazing = np.degrees(np.arctan2(abs(kz), abs(kx)))
            at_es = {
                'x': x, 'z': z, 'kx': kx, 'kz': kz,
                'theta_Es_deg': grazing,
                'group_path_km': group_path,
            }

        # ── Bubble-bottom upward crossing ─────────────────────────────────────
        if (h_bubble_bot_km is not None
                and prev_z < h_bubble_bot_km <= z
                and at_bubble is None):
            at_bubble = {
                'x': x, 'z': z, 'kx': kx, 'kz': kz,
                'A_inc': None,           # filled by hybrid model (Part 6)
                'group_path_km': group_path,
            }

        prev_z = z

        # ── Termination ───────────────────────────────────────────────────────
        if z <= 0.0 and has_risen:
            break
        if z >= z_stop:
            break

    # ── Arrival angle at ground ───────────────────────────────────────────────
    x_f, z_f, kx_f, kz_f = state
    beta_recv_deg = float(np.degrees(np.arctan2(abs(kz_f), abs(kx_f))))

    # ── Background loss ───────────────────────────────────────────────────────
    L_spread = free_space_loss_dB(group_path, freq_MHz)
    L_abs    = _martyn_absorption_dB(group_path, freq_MHz)
    L_bg_dB  = L_spread + L_abs

    return {
        'beta_deg'      : float(beta_deg),
        'trajectory'    : trajectory,
        'group_path_km' : group_path,
        'tau_ms'        : group_delay_ms(group_path),
        'h_reflect_km'  : h_reflect,
        'beta_recv_deg' : beta_recv_deg,
        'at_Es'         : at_es,
        'at_bubble'     : at_bubble,
        'L_bg_dB'       : L_bg_dB,
    }


# ── Ray fan ───────────────────────────────────────────────────────────────────

def shoot_rays_fan(tx_pos:         tuple,
                   n_model:         RefractiveIndex,
                   freq_MHz:        float = FREQ_MHZ,
                   rt_params:       dict  = RT,
                   h_Es_km:         float | None = None,
                   h_bubble_bot_km: float | None = None
                   ) -> list[dict]:
    """
    Launch n_fan rays spanning beta_min..beta_max degrees.
    Returns list of ray dicts (see module docstring).
    """
    betas = np.linspace(rt_params['beta_min'],
                        rt_params['beta_max'],
                        rt_params['n_fan'])
    return [
        trace_single_ray(tx_pos, float(b), n_model, freq_MHz,
                         rt_params, h_Es_km, h_bubble_bot_km)
        for b in betas
    ]
