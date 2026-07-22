"""Legacy local-coordinate Appleton-Hartree O/X ray tracing.

This first update2.1 implementation uses a numerical Hamiltonian finite
difference formulation.  It prioritizes explicit formulas and regression to the
isotropic 3-D tracer when ``fH_MHz = 0``.  It is retained only while C0R-5
migrates existing callers.  New geocentric/ECEF work must use
``ray_tracer_ah_3d_ecef``; this module is legacy_coordinate_reference_only.
"""

from __future__ import annotations

import numpy as np

from config import FREQ_MHZ, GEOMAG, LINK_BEARING_DEG, RE_KM, RT, UPDATE2
from utils import freq_to_k0, group_delay_ms
from models.dispersion_ah import (
    C_KM_S, appleton_hartree_roots, dispersion_D, group_velocity, select_root,
)
from .ray_tracer_3d import RefractiveIndex3D


def _b_hat_from_geomag(geomag: dict | None = None) -> np.ndarray:
    g = GEOMAG if geomag is None else geomag
    dip = np.radians(float(g.get("dip_deg", 48.7)))
    decl = np.radians(float(g.get("decl_deg", -5.5)))
    link = np.radians(float(g.get("link_bearing_deg", LINK_BEARING_DEG)))
    bx = np.cos(dip) * np.cos(decl - link)
    by = np.cos(dip) * np.sin(decl - link)
    bz = -np.sin(dip)
    b = np.array([bx, by, bz], dtype=float)
    return b / max(float(np.linalg.norm(b)), 1e-30)


class RefractiveIndexAH3D:
    """Appleton-Hartree 3-D refractive-index model."""

    def __init__(
        self,
        n3d: RefractiveIndex3D,
        freq_MHz: float = FREQ_MHZ,
        wave_mode: str = "O",
        geomag: dict | None = None,
        n2_floor: float = 1e-6,
        use_dispersion_v2: bool | None = None,
        integration_profile: str | None = None,
    ):
        self.n3d = n3d
        self.freq_MHz = float(freq_MHz)
        self.wave_mode = str(wave_mode).upper()
        self.geomag = GEOMAG if geomag is None else geomag
        self.fH_MHz = float(self.geomag.get("fH_MHz", 1.197))
        self.b_hat = _b_hat_from_geomag(self.geomag)
        self.n2_floor = float(n2_floor)
        solver = str(UPDATE2.get("AH3D_SOLVER", "v2")).lower()
        if solver not in ("v2", "legacy"):
            raise ValueError("AH3D_SOLVER must be 'v2' or 'legacy'")
        self.use_dispersion_v2 = (solver == "v2" if use_dispersion_v2 is None
                                  else bool(use_dispersion_v2))
        self.integration_profile = integration_profile
        self.earth_flat = bool(getattr(n3d, "earth_flat", True))
        self.z_km = n3d.z_km
        self.near_cutoff_n2 = max(20.0 * self.n2_floor, 2.0e-2)

    def _z_phys(self, x_km: float, y_km: float, z_cart_km: float) -> float:
        if self.earth_flat:
            z_phys = z_cart_km - (x_km * x_km + y_km * y_km) / (2.0 * RE_KM)
        else:
            z_phys = z_cart_km
        return max(float(z_phys), float(self.z_km[0]))

    def root(self, x_km, y_km, z_cart_km, k_vec=None, t=0.0, previous_n2=None):
        """Return an unclamped, classified magnetoionic root."""
        if k_vec is None:
            k_vec = np.array([1.0, 0.0, 0.0])
        Ne = max(float(self.n3d.Ne_at(x_km, y_km, z_cart_km, t=t)), 0.0)
        omega = 2.0 * np.pi * self.freq_MHz * 1e6
        roots = appleton_hartree_roots(omega, Ne, k_vec, self.b_hat, self.fH_MHz)
        root = select_root(roots, self.wave_mode, previous_n2, self.near_cutoff_n2)
        if self.earth_flat and np.isfinite(root.n2.real):
            scale = (1.0 + self._z_phys(x_km, y_km, z_cart_km) / RE_KM) ** 2
            roots = {key: value * scale for key, value in roots.items()}
            root = select_root(roots, self.wave_mode,
                               None if previous_n2 is None else previous_n2,
                               self.near_cutoff_n2)
        return root

    def n2_raw(self, x_km: float, y_km: float, z_cart_km: float, k_vec=None, t: float = 0.0) -> float:
        """Return unclamped AH phase-index squared for physical status tests.

        The numerical Hamiltonian still needs a small positive floor in n2(),
        but propagation validity is decided from this raw value following the
        Appleton-Hartree magnetoionic dispersion relation (Budden, 1985;
        Jones and Stephenson, 1975).
        """
        return float(self.root(x_km, y_km, z_cart_km, k_vec, t).n2.real)

    def n2(self, x_km: float, y_km: float, z_cart_km: float, k_vec=None, t: float = 0.0) -> float:
        """Return the physical, unclamped AH phase-index squared."""
        value = self.n2_raw(x_km, y_km, z_cart_km, k_vec=k_vec, t=t)
        if self.use_dispersion_v2:
            return value
        return max(value, self.n2_floor) if np.isfinite(value) else self.n2_floor

    def branch_status(self, x_km: float, y_km: float, z_cart_km: float, k_vec=None, t: float = 0.0) -> str:
        return self.root(x_km, y_km, z_cart_km, k_vec, t).status


def hamiltonian_ah3d(state: np.ndarray, n_ah: RefractiveIndexAH3D, k0: float, t: float = 0.0) -> float:
    x, y, z, kx, ky, kz = state
    k_vec = np.array([kx, ky, kz], dtype=float)
    return float(np.dot(k_vec, k_vec) - k0 * k0 * n_ah.n2_raw(x, y, z, k_vec=k_vec, t=t))


def hamiltonian_residual_ah3d(state: np.ndarray, n_ah: RefractiveIndexAH3D, k0: float, t: float = 0.0) -> float:
    denom = max(k0 * k0, 1e-30)
    return abs(hamiltonian_ah3d(state, n_ah, k0, t=t)) / denom


def group_index_ah3d(x_km: float, y_km: float, z_km: float, k_vec, n_ah: RefractiveIndexAH3D, t: float = 0.0) -> float:
    """Group index from the dispersion gradient v_g=-D_k/D_omega."""
    if not n_ah.use_dispersion_v2:
        f0 = n_ah.freq_MHz
        df = max(1e-3, 1e-4 * f0)
        n0 = np.sqrt(max(n_ah.n2_raw(x_km, y_km, z_km, k_vec, t), n_ah.n2_floor))
        values = []
        for freq in (f0 + df, max(f0 - df, df)):
            model = RefractiveIndexAH3D(n_ah.n3d, freq, n_ah.wave_mode,
                                        n_ah.geomag, n_ah.n2_floor, False)
            values.append(np.sqrt(max(model.n2_raw(x_km, y_km, z_km, k_vec, t), model.n2_floor)))
        return max(float(n0 + f0 * (values[0] - values[1]) / (2.0 * df)), 1e-4)
    vg = group_velocity_ah3d(x_km, y_km, z_km, k_vec, n_ah, t=t)
    speed = float(np.linalg.norm(vg))
    return float(C_KM_S / speed) if np.isfinite(speed) and speed > 0.0 else np.nan


def group_velocity_ah3d(x_km, y_km, z_km, k_vec, n_ah, t=0.0):
    Ne = max(float(n_ah.n3d.Ne_at(x_km, y_km, z_km, t=t)), 0.0)
    scale = ((1.0 + n_ah._z_phys(x_km, y_km, z_km) / RE_KM) ** 2
             if n_ah.earth_flat else 1.0)

    def d_fn(k, omega):
        roots = appleton_hartree_roots(omega, Ne, k, n_ah.b_hat, n_ah.fH_MHz)
        return dispersion_D(k, omega, roots[n_ah.wave_mode] * scale)

    return group_velocity(d_fn, k_vec, 2.0 * np.pi * n_ah.freq_MHz * 1e6)


def ah3d_deriv(state: np.ndarray, n_ah: RefractiveIndexAH3D, k0: float, t: float = 0.0) -> np.ndarray:
    """Hamiltonian finite-difference RHS normalized to geometric path length."""
    h_r = 0.5
    k_mag = max(float(np.linalg.norm(state[3:6])), 1e-30)
    h_k = max(1e-6, 1e-5 * k_mag)

    grad_r = np.empty(3, dtype=float)
    for i in range(3):
        sp = state.copy()
        sm = state.copy()
        sp[i] += h_r
        sm[i] -= h_r
        grad_r[i] = (hamiltonian_ah3d(sp, n_ah, k0, t=t) - hamiltonian_ah3d(sm, n_ah, k0, t=t)) / (2.0 * h_r)

    grad_k = np.empty(3, dtype=float)
    for i in range(3):
        sp = state.copy()
        sm = state.copy()
        sp[3 + i] += h_k
        sm[3 + i] -= h_k
        grad_k[i] = (hamiltonian_ah3d(sp, n_ah, k0, t=t) - hamiltonian_ah3d(sm, n_ah, k0, t=t)) / (2.0 * h_k)

    speed = max(float(np.linalg.norm(grad_k)), 1e-30)
    dr_ds = grad_k / speed
    dk_ds = -grad_r / speed
    return np.array([dr_ds[0], dr_ds[1], dr_ds[2], dk_ds[0], dk_ds[1], dk_ds[2]], dtype=float)


def rk4_step_ah3d(state: np.ndarray, n_ah: RefractiveIndexAH3D, k0: float, ds: float, t: float = 0.0) -> np.ndarray:
    k1 = ah3d_deriv(state, n_ah, k0, t=t)
    k2 = ah3d_deriv(state + 0.5 * ds * k1, n_ah, k0, t=t)
    k3 = ah3d_deriv(state + 0.5 * ds * k2, n_ah, k0, t=t)
    k4 = ah3d_deriv(state + ds * k3, n_ah, k0, t=t)
    return state + (ds / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _renorm_k_ah3d(state: np.ndarray, n_ah: RefractiveIndexAH3D, k0: float, t: float = 0.0) -> np.ndarray:
    k_vec = state[3:6]
    k_mag = max(float(np.linalg.norm(k_vec)), 1e-30)
    n2 = n_ah.n2_raw(float(state[0]), float(state[1]), float(state[2]), k_vec=k_vec, t=t)
    if not np.isfinite(n2) or n2 <= 0.0:
        return state.copy()
    n_eff = np.sqrt(n2)
    out = state.copy()
    out[3:6] *= (n_eff * k0) / k_mag
    return out


def _locate_cutoff(prev, trial, n_ah, t=0.0, iterations=30):
    """Locate n2=0 on one accepted-to-evanescent step by bisection."""
    lo, hi = prev.copy(), trial.copy()
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        value = n_ah.n2_raw(*mid[:3], k_vec=mid[3:6], t=t)
        if np.isfinite(value) and value > 0.0:
            lo = mid
        else:
            hi = mid
    return lo


def trace_single_ray_ah_3d(
    tx_pos: tuple,
    beta_deg: float,
    az_deg: float,
    n_ah: RefractiveIndexAH3D,
    freq_MHz: float = FREQ_MHZ,
    rt_params: dict = RT,
    t: float = 0.0,
    with_sensitivity: bool = False,
) -> dict:
    """Trace one AH3D O/X ray."""
    if n_ah.use_dispersion_v2:
        from .ray_tracer_ah_3d_v2 import trace_single_ray_ah_3d_v2
        return trace_single_ray_ah_3d_v2(
            tx_pos, beta_deg, az_deg, n_ah, freq_MHz, rt_params, t=t,
            profile=n_ah.integration_profile,
            with_sensitivity=with_sensitivity,
        )
    if with_sensitivity:
        raise ValueError("landing sensitivity requires AH3D_SOLVER='v2'")
    k0 = freq_to_k0(freq_MHz)
    ds_min = float(rt_params["ds_min_km"])
    ds_max = float(rt_params["ds_max_km"])
    z_stop = float(rt_params["z_stop_km"])
    max_st = int(rt_params["max_steps"])

    x0, z0 = float(tx_pos[0]), float(tx_pos[1])
    y0 = 0.0
    beta_r = np.radians(beta_deg)
    az_r = np.radians(az_deg)
    launch_hat = np.array([
        np.cos(beta_r) * np.cos(az_r),
        np.cos(beta_r) * np.sin(az_r),
        np.sin(beta_r),
    ])
    n0 = np.sqrt(n_ah.n2(x0, y0, z0, k_vec=launch_hat, t=t))
    state = np.array([x0, y0, z0, *(n0 * k0 * launch_hat)], dtype=float)

    traj = [state.copy()]
    geom_path = 0.0
    group_path = 0.0
    h_reflect = z0
    x_reflect = x0
    y_reflect = y0
    has_risen = False
    statuses = set()
    status_history = []
    branch_gap_min = float("inf")
    branch_continuity_max = 0.0
    previous_root = None
    n2_min_raw = float("inf")
    max_hamiltonian_residual = hamiltonian_residual_ah3d(state, n_ah, k0, t=t)

    step_count = 0
    for step_count in range(max_st):
        x, y, z, kx, ky, kz = state
        root = n_ah.root(
            x, y, z, k_vec=state[3:6], t=t,
            previous_n2=previous_root if n_ah.use_dispersion_v2 else None,
        )
        step_status = root.status
        previous_root = root.n2
        branch_gap_min = min(branch_gap_min, root.branch_gap)
        branch_continuity_max = max(branch_continuity_max, root.branch_continuity)
        statuses.add(step_status)
        if not status_history or status_history[-1] != step_status:
            status_history.append(step_status)
        n2_raw_here = n_ah.n2_raw(x, y, z, k_vec=state[3:6], t=t)
        if np.isfinite(n2_raw_here):
            n2_min_raw = min(n2_min_raw, float(n2_raw_here))
        if step_status in ("invalid_dispersion", "resonance_excluded", "branch_lost", "evanescent"):
            break
        n2_here = n_ah.n2(x, y, z, k_vec=state[3:6], t=t)
        ds = ds_max if n2_here > 0.20 else (max(ds_min, 0.4) if n2_here > 0.04 else ds_min)
        prev = state.copy()
        trial = rk4_step_ah3d(state, n_ah, k0, ds, t=t)
        if not np.all(np.isfinite(trial)):
            statuses.add("step_failed")
            status_history.append("step_failed")
            state = prev
            break
        trial_status = n_ah.branch_status(*trial[:3], k_vec=trial[3:6], t=t)
        if (n_ah.use_dispersion_v2
                and trial_status in ("evanescent", "invalid_dispersion", "resonance_excluded")):
            state = _locate_cutoff(prev, trial, n_ah, t=t)
            statuses.add("turning_point")
            status_history.append("turning_point")
            traj.append(state.copy())
            break
        state = trial
        state = _renorm_k_ah3d(state, n_ah, k0, t=t)

        dx = float(state[0] - prev[0])
        dy = float(state[1] - prev[1])
        dz = float(state[2] - prev[2])
        ds_geom = (dx * dx + dy * dy + dz * dz) ** 0.5
        x_mid = float(0.5 * (state[0] + prev[0]))
        y_mid = float(0.5 * (state[1] + prev[1]))
        z_mid = float(0.5 * (state[2] + prev[2]))
        k_mid = 0.5 * (state[3:6] + prev[3:6])
        ng_mid = group_index_ah3d(x_mid, y_mid, z_mid, k_mid, n_ah, t=t)
        if not np.isfinite(ng_mid) or ng_mid <= 0.0:
            statuses.add("invalid_dispersion")
            status_history.append("invalid_dispersion")
            state = prev
            break
        geom_path += ds_geom
        group_path += ds_geom * ng_mid
        max_hamiltonian_residual = max(max_hamiltonian_residual, hamiltonian_residual_ah3d(state, n_ah, k0, t=t))

        x, y, z, kx, ky, kz = state
        traj.append(state.copy())
        if z > 50.0:
            has_risen = True
        if z > h_reflect:
            h_reflect = float(z)
            x_reflect = float(x)
            y_reflect = float(y)
        if z <= 0.0 and has_risen:
            break
        if z >= z_stop:
            break

    x_f, y_f, z_f, kx_f, ky_f, kz_f = state
    kmag_f = max(float(np.sqrt(kx_f * kx_f + ky_f * ky_f + kz_f * kz_f)), 1e-30)
    beta_rx = float(np.degrees(np.arcsin(np.clip(abs(kz_f) / kmag_f, 0.0, 1.0))))
    az_rx = float(np.degrees(np.arctan2(ky_f, kx_f)))
    if n_ah.earth_flat:
        h_phys = h_reflect - (x_reflect * x_reflect + y_reflect * y_reflect) / (2.0 * RE_KM)
    else:
        h_phys = h_reflect
    az_deflect = float(np.degrees(np.arctan2(y_f, max(abs(x_f), 1e-6))))

    if not np.isfinite(n2_min_raw):
        n2_min_raw = float("nan")

    if z_f > 0.0 and not (z_f <= 0.0 and has_risen):
        statuses.add("mode_lost")
        if not status_history or status_history[-1] != "mode_lost":
            status_history.append("mode_lost")

    branch_status = "ok"
    for status in ("step_failed", "invalid_dispersion", "resonance_excluded",
                   "branch_lost", "evanescent", "mode_lost", "turning_point",
                   "near_cutoff"):
        if status in statuses:
            branch_status = status
            break

    return {
        "x_land_km": float(x_f),
        "y_land_km": float(y_f),
        "z_land_km": float(z_f),
        "geom_path_km": float(geom_path),
        "group_path_km": float(group_path),
        "tau_ms": group_delay_ms(group_path),
        "h_reflect_km": float(h_phys),
        "reflect_x_km": float(x_reflect),
        "reflect_y_km": float(y_reflect),
        "beta_deg_tx": float(beta_deg),
        "az_deg_tx": float(az_deg),
        "elevation_tx_deg": float(beta_deg),
        "beta_deg_rx": beta_rx,
        "elevation_rx_deg": beta_rx,
        "az_deg_rx": az_rx,
        "azimuth_deflect_deg": az_deflect,
        "trajectory_3d": traj,
        "wave_mode": n_ah.wave_mode,
        "branch_status": branch_status,
        "branch_status_history": ",".join(status_history),
        "n2_min_raw": float(n2_min_raw),
        "group_index_model": ("dispersion_gradient_-Dk_Domega" if n_ah.use_dispersion_v2
                              else "legacy_fd_dn_df"),
        "group_path_definition": ("c_integral_dt" if n_ah.use_dispersion_v2
                                  else "integral_n_g_ds"),
        "ray_attenuation_semantics": "real_ray_complex_attenuation",
        "branch_id": n_ah.wave_mode,
        "branch_gap_min": float(branch_gap_min),
        "branch_continuity_max": float(branch_continuity_max),
        "n2_floor_propagation": not n_ah.use_dispersion_v2,
        "hamiltonian_residual_max": float(max_hamiltonian_residual),
        "ox_solver": "ah3d",
        "freq_MHz": float(freq_MHz),
        "theory_reference": "Appleton-Hartree; Haselgrove Hamiltonian; Jones-Stephenson 3D ray tracing",
    }


def shoot_rays_fan_ah_3d(
    tx_pos: tuple,
    n_ah: RefractiveIndexAH3D,
    freq_MHz: float = FREQ_MHZ,
    rt_params: dict = RT,
    az_center_deg: float = 0.0,
    t: float = 0.0,
) -> list:
    betas = np.linspace(rt_params["beta_min"], rt_params["beta_max"], rt_params["n_fan"])
    return [
        trace_single_ray_ah_3d(tx_pos, float(beta), az_center_deg, n_ah, freq_MHz, rt_params, t=t)
        for beta in betas
    ]
