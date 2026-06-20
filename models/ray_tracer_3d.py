"""
Module: 3-D HF ray tracer (Cartesian-on-sphere).

Implements the 3-D Haselgrove (1955) equations in Cartesian coordinates
with Earth-curvature correction (Earth-flattening transform) and optional
3-D TID perturbation including cross-path (y) variation.

State vector:   [x, y, z, kx, ky, kz]   (km, km, km, km^-1, km^-1, km^-1)
Integration parameter: group path P [km]

Haselgrove / Jones & Stephenson (1975) Cartesian form:
  dx/dP = kx / |k|,  dy/dP = ky / |k|,  dz/dP = kz / |k|
  dkx/dP = (k0^2 / 2|k|) * dn^2_eff / dx
  dky/dP = (k0^2 / 2|k|) * dn^2_eff / dy
  dkz/dP = (k0^2 / 2|k|) * dn^2_eff / dz

where |k| = n_eff * k0 (on-shell condition H = |k|^2 - n^2 k0^2 = 0).

Earth-flattening (Jones 1975, PHaRLAP; first-order parabolic, valid <1500 km):
  z_phys = z_cart - (x^2 + y^2) / (2 * R_E)
  n^2_eff(x, y, z_cart) = n^2_plasma(x, z_phys) * (1 + z_phys/R_E)^2

Coordinate system:
  x : along great-circle from TX to RX  [km]
  y : cross-path (perpendicular, horizontal)  [km]
  z : "flattened" height (approximation to height above ellipsoid)  [km]

Reference: Jones & Stephenson (1975); Cervera & Harris (2014 JGR).
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from config import (C_KMS, FREQ_MHZ, RT, RE_KM, K_FP,
                    IRI_DT, IRI_LAT, IRI_LON, LINK_BEARING_DEG,
                    TX_LAT, TX_LON)
from utils import freq_to_k0, group_delay_ms, destination_point


# ── 3-D refractive-index model ────────────────────────────────────────────────

class RefractiveIndex3D:
    """
    Position-dependent refractive index n^2(x, y, z) for 3-D ray tracing.

    Background (no y-variation): 2-D IRI Ne(x, z) grid broadcast in y.
    TID perturbation (Phase 5/6): 3-D δNe(x, y, z, t) with full wave vector
      including cross-path component ky = k_h * sin(az_TID - az_link).
    Earth-flattening: z_phys = z_cart - (x^2+y^2)/(2*R_E)  (Jones 1975).

    Parameters
    ----------
    Ne_2d      : (Nx, Nz) IRI background electron density [m^-3]
    x_km       : (Nx,)  horizontal positions along link [km]
    z_km       : (Nz,)  heights [km]
    freq_MHz   : float  operating frequency [MHz]
    tid_params : dict   from config.TID (optional; None = no TID)
    earth_flat : bool   apply Earth-flattening correction (default True)
    """

    def __init__(self,
                 Ne_2d:      np.ndarray,
                 x_km:       np.ndarray,
                 z_km:       np.ndarray,
                 freq_MHz:   float = FREQ_MHZ,
                 tid_params: dict | None = None,
                 earth_flat: bool = True):
        self.freq_MHz   = freq_MHz
        self.Ne_2d      = Ne_2d
        self.x_km       = x_km
        self.z_km       = z_km
        self.earth_flat = earth_flat
        self.tid        = tid_params

        freq_Hz  = freq_MHz * 1e6
        fp2_2d   = (K_FP ** 2) * Ne_2d                        # [Hz^2]
        n2_2d    = np.maximum(1.0 - fp2_2d / freq_Hz ** 2, 1e-6)

        self._n2_interp = RegularGridInterpolator(
            (x_km, z_km), n2_2d,
            method='linear', bounds_error=False, fill_value=1.0)
        self._Ne_interp = RegularGridInterpolator(
            (x_km, z_km), Ne_2d,
            method='linear', bounds_error=False, fill_value=0.0)

        # Pre-compute TID normalization factors (from 2D cross-section at y=0)
        self._tid_norms = []
        if self.tid is not None and self.tid.get('enable', False):
            self._tid_norms = self._precompute_tid_norms()

    def _precompute_tid_norms(self) -> list:
        """
        Compute per-component normalization factors for the Hooke (1968) TID.
        Returns list of (norm_factor, kx_i, ky_i, kz_i, omega_i, k_para_i) tuples.
        """
        p       = self.tid
        n_comp  = int(p.get('n_components', 1))
        omega_b = float(p.get('omega_b_rad_s', 2.0 * np.pi / 1200.0))
        I_rad   = np.deg2rad(float(p.get('I_dip_deg', 50.0)))
        H_m     = float(p.get('H_km', 60.0)) * 1e3
        link_az = float(p.get('link_bearing_deg', LINK_BEARING_DEG))

        az_list  = p.get('az_deg_list',      [link_az])
        amp_list = p.get('amplitude_list',   [float(p.get('amplitude', 0.1))])
        T_list   = p.get('period_s_list',    [float(p.get('T_s', 2400.0))])
        lam_list = p.get('lambda_h_km_list', [float(p.get('lambda_h_km', 300.0))])

        norms = []
        for i in range(n_comp):
            lam_h_m = float(lam_list[i]) * 1e3
            k_h     = 2.0 * np.pi / lam_h_m               # full horizontal |k| [rad/m]
            T_s     = float(T_list[i])
            omega   = 2.0 * np.pi / T_s                    # [rad/s]
            amp     = float(amp_list[i])

            kz2 = k_h**2 * (omega_b**2 / omega**2 - 1.0) - 1.0/(4.0*H_m**2)
            if kz2 <= 0.0:
                norms.append(None)
                continue
            kz = np.sqrt(kz2)                              # vertical wavenumber [rad/m]

            # 3-D wave vector components
            az_rel  = np.radians(float(az_list[i]) - link_az)
            kx_i    = k_h * np.cos(az_rel)                # along-link [rad/m]
            ky_i    = k_h * np.sin(az_rel)                # cross-link [rad/m]

            # k_para: projection onto geomagnetic B direction
            # B_hat ~ (cos(I)*cos(0), cos(I)*sin(0), sin(I)) in (x, y, z) frame
            # For the Hooke formula, the dominant term uses |k_h| and kz:
            # k_para = k_h * cos(I) + kz * sin(I)  (AGW dispersion projection)
            k_para  = k_h * np.cos(I_rad) + kz * np.sin(I_rad)
            sinI    = np.sin(I_rad)

            # Compute normalization on 2-D slice (y=0) at t=0
            dz_m      = (self.z_km[1] - self.z_km[0]) * 1e3
            Ne_mid    = self.Ne_2d[len(self.x_km)//2, :]       # mid-path profile
            dNe_dz    = np.gradient(Ne_mid, dz_m)
            X_m = self.x_km * 1e3
            Z_m = self.z_km * 1e3
            phase_2d  = kx_i * X_m[:, None] + kz * Z_m[None, :] # (Nx, Nz) at y=0,t=0
            dNe_raw   = (-dNe_dz[None, :] * sinI * np.sin(phase_2d)
                         - k_para * Ne_mid[None, :] * np.cos(phase_2d))

            Ne_peak = Ne_mid.max()
            f_mask  = Ne_mid > Ne_peak * 0.01
            if not f_mask.any() or Ne_peak < 1.0:
                norms.append(None)
                continue
            max_ratio = np.max(np.abs(dNe_raw[:, f_mask]) / Ne_mid[None, f_mask])
            if max_ratio < 1e-30:
                norms.append(None)
                continue

            norms.append((amp / max_ratio, kx_i, ky_i, kz, omega, k_para, sinI))

        return norms

    def _tid_delta_Ne(self, x_km: float, y_km: float, z_phys_km: float,
                      t: float = 0.0) -> float:
        """
        Evaluate total TID δNe(x, y, z, t) [m^-3] using Hooke (1968) formula
        with 3-D wave vector including cross-path ky component.
        """
        if not self._tid_norms:
            return 0.0

        I_rad  = np.deg2rad(float(self.tid.get('I_dip_deg', 50.0)))
        dz_m   = (self.z_km[1] - self.z_km[0]) * 1e3
        z_m    = z_phys_km * 1e3

        # Look up local Ne and its vertical gradient from 2-D IRI grid
        pt    = np.array([[x_km, z_phys_km]])
        Ne_bg = float(self._Ne_interp(pt)[0])
        # Finite-difference dNe/dz at this point
        pt_up  = np.array([[x_km, z_phys_km + dz_m/1e3]])
        pt_dn  = np.array([[x_km, z_phys_km - dz_m/1e3]])
        dNe_dz = (float(self._Ne_interp(pt_up)[0]) -
                  float(self._Ne_interp(pt_dn)[0])) / (2.0 * dz_m)

        dNe_total = 0.0
        x_m = x_km * 1e3
        y_m = y_km * 1e3

        for entry in self._tid_norms:
            if entry is None:
                continue
            norm, kx_i, ky_i, kz_i, omega_i, k_para_i, sinI = entry
            phase = kx_i * x_m + ky_i * y_m + kz_i * z_m - omega_i * t
            dNe_i = (-dNe_dz * sinI * np.sin(phase)
                     - k_para_i * Ne_bg * np.cos(phase))
            dNe_total += norm * dNe_i

        return dNe_total

    def n2(self, x_km: float, y_km: float, z_cart_km: float,
           t: float = 0.0) -> float:
        """
        n^2_eff(x, y, z_cart) with Earth-flattening and optional 3-D TID.

        Steps:
          1. z_phys = z_cart - (x^2+y^2)/(2*R_E)   [km]
          2. n^2_plasma from 2-D IRI at (x, z_phys)
          3. Add TID δNe contribution
          4. Apply Earth-flattening factor (1 + z_phys/R_E)^2
        """
        if self.earth_flat:
            z_phys = z_cart_km - (x_km**2 + y_km**2) / (2.0 * RE_KM)
        else:
            z_phys = z_cart_km
        z_phys = max(z_phys, float(self.z_km[0]))

        # Base n^2 from IRI
        pt      = np.array([[x_km, z_phys]])
        n2_base = float(self._n2_interp(pt)[0])

        # TID perturbation
        if self.tid is not None and self.tid.get('enable', False):
            freq_Hz = self.freq_MHz * 1e6
            Ne_bg   = float(self._Ne_interp(pt)[0])
            dNe     = self._tid_delta_Ne(x_km, y_km, z_phys, t)
            Ne_tot  = max(Ne_bg + dNe, 0.0)
            fp2_tot = (K_FP**2) * Ne_tot
            n2_base = max(1.0 - fp2_tot / freq_Hz**2, 1e-6)

        # Earth-flattening correction
        if self.earth_flat:
            ef = (1.0 + z_phys / RE_KM) ** 2
            n2_base *= ef

        return max(n2_base, 1e-6)

    def grad_n2(self, x_km: float, y_km: float, z_cart_km: float,
                dx: float = 0.5, dy: float = 0.5, dz: float = 0.5,
                t: float = 0.0) -> tuple[float, float, float]:
        """
        Central-difference gradient (dn^2/dx, dn^2/dy, dn^2/dz) [km^-1].
        """
        dn2_dx = (self.n2(x_km+dx, y_km,    z_cart_km,    t)
                - self.n2(x_km-dx, y_km,    z_cart_km,    t)) / (2.0*dx)
        dn2_dy = (self.n2(x_km,    y_km+dy, z_cart_km,    t)
                - self.n2(x_km,    y_km-dy, z_cart_km,    t)) / (2.0*dy)
        dn2_dz = (self.n2(x_km,    y_km,    z_cart_km+dz, t)
                - self.n2(x_km,    y_km,    z_cart_km-dz, t)) / (2.0*dz)
        return dn2_dx, dn2_dy, dn2_dz

    def Ne_at(self, x_km: float, y_km: float, z_cart_km: float,
              t: float = 0.0) -> float:
        """Background Ne [m^-3] at (x, y, z_cart) with Earth-flattening."""
        if self.earth_flat:
            z_phys = z_cart_km - (x_km**2 + y_km**2) / (2.0 * RE_KM)
        else:
            z_phys = z_cart_km
        z_phys = max(z_phys, float(self.z_km[0]))
        pt = np.array([[x_km, z_phys]])
        return max(float(self._Ne_interp(pt)[0]), 0.0)


# ── 3-D Haselgrove ODE ───────────────────────────────────────────────────────

def haselgrove_3d_deriv(state: np.ndarray,
                         n3d:   RefractiveIndex3D,
                         k0:    float,
                         t:     float = 0.0) -> np.ndarray:
    """
    RHS of 3-D Haselgrove equations.
    State: [x, y, z, kx, ky, kz]  (km, km, km, km^-1, km^-1, km^-1)
    Returns d(state)/dP  (P = group path [km]).
    """
    x, y, z, kx, ky, kz = state
    dn2_dx, dn2_dy, dn2_dz = n3d.grad_n2(x, y, z, t=t)
    kmag = max(np.sqrt(kx**2 + ky**2 + kz**2), 1e-30)
    half_k02_km = 0.5 * k0**2 / kmag
    return np.array([
        kx / kmag,
        ky / kmag,
        kz / kmag,
        half_k02_km * dn2_dx,
        half_k02_km * dn2_dy,
        half_k02_km * dn2_dz,
    ])


def rk4_step_3d(state: np.ndarray,
                n3d:   RefractiveIndex3D,
                k0:    float,
                ds:    float,
                t:     float = 0.0) -> np.ndarray:
    """Single RK4 step of length ds [km]. Returns new state."""
    k1 = haselgrove_3d_deriv(state,             n3d, k0, t)
    k2 = haselgrove_3d_deriv(state + ds/2 * k1, n3d, k0, t)
    k3 = haselgrove_3d_deriv(state + ds/2 * k2, n3d, k0, t)
    k4 = haselgrove_3d_deriv(state + ds   * k3, n3d, k0, t)
    return state + (ds / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


def _renorm_k(state: np.ndarray, n3d: RefractiveIndex3D, k0: float) -> np.ndarray:
    """
    Renormalize wave-vector components to maintain on-shell condition |k| = n*k0.
    PHaRLAP applies this after each RK4 step to prevent drift.
    """
    x, y, z, kx, ky, kz = state
    n_eff  = np.sqrt(n3d.n2(x, y, z))
    kmag   = np.sqrt(kx**2 + ky**2 + kz**2)
    target = n_eff * k0
    if kmag > 1e-30:
        scale = target / kmag
        state = state.copy()
        state[3] *= scale
        state[4] *= scale
        state[5] *= scale
    return state


# ── Single 3-D ray tracer ─────────────────────────────────────────────────────

def trace_single_ray_3d(tx_pos:    tuple,
                         beta_deg:  float,
                         az_deg:    float,
                         n3d:       RefractiveIndex3D,
                         freq_MHz:  float = FREQ_MHZ,
                         rt_params: dict  = RT,
                         t:         float = 0.0
                         ) -> dict:
    """
    Trace one 3-D ray from tx_pos at (elevation, azimuth).

    Parameters
    ----------
    tx_pos    : (x0, z0) in km  (y=0 assumed; TX on link axis)
    beta_deg  : launch elevation above horizon [deg]
    az_deg    : azimuth offset from link direction (+y direction) [deg]
                az_deg=0 means ray stays in x-z plane (2-D equivalent)
    n3d       : RefractiveIndex3D instance
    freq_MHz  : operating frequency [MHz]
    rt_params : RT config dict
    t         : snapshot time for TID [s]

    Returns dict with fields:
      x_land_km, y_land_km, z_land_km : landing position
      group_path_km, tau_ms           : propagation delay
      h_reflect_km                    : peak height above Earth (physical)
      reflect_x_km, reflect_y_km     : horizontal position at reflection
      beta_deg_tx                     : TX elevation (input)
      az_deg_tx                       : TX azimuth (input)
      beta_deg_rx                     : arrival elevation at ground
      az_deg_rx                       : arrival azimuth at ground
      azimuth_deflect_deg             : lateral deflection y_land/x_land [deg]
      trajectory_3d                   : list of [x,y,z,kx,ky,kz] states
    """
    k0      = freq_to_k0(freq_MHz)
    ds_nom  = rt_params['ds_km']
    ds_min  = rt_params['ds_min_km']
    ds_max  = rt_params['ds_max_km']
    z_stop  = rt_params['z_stop_km']
    max_st  = rt_params['max_steps']

    x0, z0  = float(tx_pos[0]), float(tx_pos[1])
    y0      = 0.0
    beta_r  = np.radians(beta_deg)
    az_r    = np.radians(az_deg)

    n0   = np.sqrt(n3d.n2(x0, y0, z0, t=t))
    kmag = n0 * k0
    # Initial wave vector: elevation beta from horizontal, azimuth az from x-axis
    kx0 = kmag * np.cos(beta_r) * np.cos(az_r)
    ky0 = kmag * np.cos(beta_r) * np.sin(az_r)
    kz0 = kmag * np.sin(beta_r)
    state = np.array([x0, y0, z0, kx0, ky0, kz0], dtype=float)

    traj       = [state.copy()]
    group_path = 0.0
    h_reflect  = z0     # maximum z (Cartesian), used for reflection height
    x_reflect  = x0
    y_reflect  = y0
    has_risen  = False

    for _ in range(max_st):
        x, y, z, kx, ky, kz = state
        n2_here = n3d.n2(x, y, z, t=t)

        # Adaptive step: shrink near reflection
        if n2_here < 0.04:
            ds = ds_min
        elif n2_here < 0.20:
            ds = max(ds_min, ds_nom * 0.4)
        else:
            ds = ds_nom

        state      = rk4_step_3d(state, n3d, k0, ds, t=t)
        state      = _renorm_k(state, n3d, k0)        # PHaRLAP on-shell renorm
        group_path += ds

        x, y, z, kx, ky, kz = state
        traj.append(state.copy())

        if z > 50.0:
            has_risen = True
        if z > h_reflect:
            h_reflect = z
            x_reflect = x
            y_reflect = y

        if z <= 0.0 and has_risen:
            break
        if z >= z_stop:
            break

    x_f, y_f, z_f, kx_f, ky_f, kz_f = state
    kmag_f = max(np.sqrt(kx_f**2 + ky_f**2 + kz_f**2), 1e-30)

    # Arrival angles
    beta_rx = float(np.degrees(np.arcsin(np.clip(abs(kz_f)/kmag_f, 0, 1))))
    az_rx   = float(np.degrees(np.arctan2(ky_f, kx_f)))

    # Physical reflection height (correct for Earth curvature)
    if n3d.earth_flat:
        h_phys = h_reflect - (x_reflect**2 + y_reflect**2) / (2.0 * RE_KM)
    else:
        h_phys = h_reflect

    # Azimuth deflection: angle of landing point off link axis
    az_deflect = float(np.degrees(np.arctan2(y_f, max(abs(x_f), 1e-6))))

    return {
        'x_land_km'        : float(x_f),
        'y_land_km'        : float(y_f),
        'z_land_km'        : float(z_f),
        'group_path_km'    : group_path,
        'tau_ms'           : group_delay_ms(group_path),
        'h_reflect_km'     : float(h_phys),
        'reflect_x_km'     : float(x_reflect),
        'reflect_y_km'     : float(y_reflect),
        'beta_deg_tx'      : float(beta_deg),
        'az_deg_tx'        : float(az_deg),
        'beta_deg_rx'      : beta_rx,
        'az_deg_rx'        : az_rx,
        'azimuth_deflect_deg': az_deflect,
        'trajectory_3d'    : traj,
    }


# ── 3-D ray fan ───────────────────────────────────────────────────────────────

def shoot_rays_fan_3d(tx_pos:    tuple,
                       n3d:       RefractiveIndex3D,
                       freq_MHz:  float = FREQ_MHZ,
                       rt_params: dict  = RT,
                       az_deg:    float = 0.0,
                       t:         float = 0.0
                       ) -> list[dict]:
    """
    Launch a fan of 3-D rays at fixed azimuth, scanning elevation.
    az_deg=0 matches the 2-D fan for backward compatibility.
    """
    betas = np.linspace(rt_params['beta_min'],
                        rt_params['beta_max'],
                        rt_params['n_fan'])
    return [
        trace_single_ray_3d(tx_pos, float(b), az_deg, n3d,
                            freq_MHz, rt_params, t)
        for b in betas
    ]
