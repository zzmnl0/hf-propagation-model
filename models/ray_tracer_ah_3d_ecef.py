"""ECEF Appleton-Hartree Hamiltonian core for spherical and WGS84 Earth.

This module deliberately has no local-Cartesian or Earth-flattening entry.
Ground/upper-boundary events and downstream migration belong to C0R-4/C0R-5.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.integrate import DOP853

from config import K_FP, UPDATE2
from utils import (
    ecef_to_geocentric,
    ecef_to_geodetic,
    enu_basis_wgs84,
    geocentric_to_ecef,
    geodetic_to_ecef,
    launch_direction_ecef,
    launch_direction_ecef_wgs84,
)
from .dispersion_ah import (
    C_KM_S,
    AHRoot,
    appleton_hartree_roots,
    dispersion_D,
    dispersion_partials,
    select_root,
)
from .ne_provider_geocentric import COORDINATE_SYSTEM, GeocentricNeProvider3D
from .geomagnetic_field_3d import (
    GYRO_MHZ_PER_T,
    ConstantGeomagneticField3D,
)


COORDINATE_MODEL = "geocentric_spherical_ecef"
WGS84_COORDINATE_MODEL = "wgs84_geodetic_ecef"
STATE_SEMANTICS = "r_ecef_km,k_ecef_km^-1,t_group_s"
_COORDINATE_MAPPING_HASH = hashlib.sha256(
    b"C0R-3:geocentric_spherical<->ECEF:R_explicit:v1"
).hexdigest()


class SpatialDerivativeError(RuntimeError):
    """Raised when an ECEF spatial derivative is non-finite or unresolved."""


def _metadata_hash(metadata: dict) -> str:
    return hashlib.sha256(json.dumps(
        metadata, sort_keys=True, separators=(",", ":"), default=str,
    ).encode("utf-8")).hexdigest()


class RefractiveIndexAH3DECEF:
    """Real cold-plasma AH dispersion model queried only in physical ECEF."""

    def __init__(
        self,
        provider: GeocentricNeProvider3D,
        *,
        earth_radius_km: float,
        freq_MHz: float,
        wave_mode: str = "O",
        gyro_MHz: float = 0.0,
        b_ecef_hat=(0.0, 0.0, 1.0),
        geomagnetic_provider=None,
        coordinate_model: str = COORDINATE_MODEL,
        collision_s: float = 0.0,
        near_cutoff_n2: float = 2.0e-2,
        spatial_fd_base_km: float | None = None,
        spatial_fd_scales=(0.5, 1.0, 2.0),
        spatial_fd_stability_rtol: float = 2.0e-2,
        spatial_fd_stability_atol_n2_per_km: float = 2.0e-5,
        integration_profile: str = "reference",
        plasma_coefficient: float | None = None,
        compatibility_mode: str | None = None,
    ):
        if not (hasattr(provider, "get_Ne_ecef")
                and hasattr(provider, "domain_margin_ecef")
                and hasattr(provider, "get_metadata")):
            raise TypeError("provider must implement the physical ECEF Ne contract")
        self.provider = provider
        self.earth_radius_km = float(earth_radius_km)
        self.coordinate_model = str(coordinate_model)
        if self.coordinate_model not in {
                COORDINATE_MODEL, WGS84_COORDINATE_MODEL}:
            raise ValueError("coordinate_model must be spherical ECEF or WGS84 ECEF")
        self.freq_MHz = float(freq_MHz)
        self.wave_mode = str(wave_mode).upper()
        self.gyro_MHz = float(gyro_MHz)
        self.collision_s = float(collision_s)
        self.near_cutoff_n2 = float(near_cutoff_n2)
        if compatibility_mode not in {None, "pharlap_radial_b0"}:
            raise ValueError("unsupported compatibility_mode")
        self.compatibility_mode = compatibility_mode
        compat = UPDATE2["PHARLAP_RADIAL_B0"] if compatibility_mode else {}
        self.plasma_coefficient = float(
            compat.get("plasma_coefficient", K_FP**2)
            if plasma_coefficient is None else plasma_coefficient)
        if not np.isfinite(self.plasma_coefficient) or self.plasma_coefficient <= 0.0:
            raise ValueError("plasma_coefficient must be finite and positive")
        self.spatial_fd_base_km = float(
            compat.get("spatial_fd_base_km", 0.5)
            if spatial_fd_base_km is None else spatial_fd_base_km)
        self.spatial_fd_scales = tuple(float(value) for value in spatial_fd_scales)
        self.spatial_fd_stability_rtol = float(spatial_fd_stability_rtol)
        self.spatial_fd_stability_atol_n2_per_km = float(
            spatial_fd_stability_atol_n2_per_km)
        self.integration_profile = str(integration_profile)
        self.omega = 2.0 * np.pi * self.freq_MHz * 1.0e6
        # Use the same exact c as dispersion_D and v_g=-D_k/D_omega.
        self.k0 = self.omega / C_KM_S

        metadata = provider.get_metadata()
        provider_radius = metadata.get("earth_radius_km")
        expected_provider_coordinates = (
            COORDINATE_SYSTEM if self.coordinate_model == COORDINATE_MODEL
            else "wgs84_geodetic")
        if metadata.get("coordinate_system") != expected_provider_coordinates:
            raise ValueError(
                "Provider coordinates do not match the AH3D coordinate model")
        if (self.coordinate_model == COORDINATE_MODEL
                and (provider_radius is None
                     or float(provider_radius) != self.earth_radius_km)):
            raise ValueError("Provider and AH3D Earth radii must match exactly")
        if not np.isfinite(self.freq_MHz) or self.freq_MHz <= 0.0:
            raise ValueError("freq_MHz must be finite and positive")
        if self.wave_mode not in {"O", "X"}:
            raise ValueError("wave_mode must be 'O' or 'X'")
        if not np.isfinite(self.gyro_MHz) or self.gyro_MHz < 0.0:
            raise ValueError("gyro_MHz must be finite and nonnegative")
        if self.collision_s != 0.0:
            raise ValueError(
                "C0R-3 traces the real collisionless dispersion surface; "
                "collisional attenuation is a separate path integral")
        b_hat = np.asarray(b_ecef_hat, dtype=float)
        if b_hat.shape != (3,) or not np.all(np.isfinite(b_hat)):
            raise ValueError("b_ecef_hat must be a finite length-3 ECEF vector")
        b_norm = float(np.linalg.norm(b_hat))
        if b_norm <= 0.0:
            raise ValueError("b_ecef_hat must be nonzero")
        self.b_ecef_hat = b_hat / b_norm
        if geomagnetic_provider is None:
            constant_magnitude_T = self.gyro_MHz / GYRO_MHZ_PER_T
            geomagnetic_provider = ConstantGeomagneticField3D(
                constant_magnitude_T * self.b_ecef_hat)
        if not (hasattr(geomagnetic_provider, "get_B_ecef_T")
                and hasattr(geomagnetic_provider, "get_metadata")):
            raise TypeError(
                "geomagnetic_provider must implement get_B_ecef_T and get_metadata")
        self.geomagnetic_provider = geomagnetic_provider
        if compatibility_mode:
            if (self.coordinate_model != COORDINATE_MODEL or self.earth_radius_km != 6371.0
                    or metadata.get("interpolation") != "radial_hermite"):
                raise ValueError("radial compatibility requires a radial Hermite Provider on the 6371 km sphere")
            if (self.gyro_MHz != 0.0 or self.collision_s != 0.0
                    or not isinstance(geomagnetic_provider, ConstantGeomagneticField3D)
                    or np.any(geomagnetic_provider.b_ecef_T != 0.0)):
                raise ValueError("radial compatibility requires zero constant magnetic field and collisions")
        self.geomagnetic_metadata = geomagnetic_provider.get_metadata()
        self.geomagnetic_provider_hash = str(
            self.geomagnetic_metadata.get(
                "field_hash", _metadata_hash(self.geomagnetic_metadata)))
        if (not np.isfinite(self.spatial_fd_base_km)
                or self.spatial_fd_base_km <= 0.0):
            raise ValueError("spatial_fd_base_km must be finite and positive")
        if (len(self.spatial_fd_scales) < 3
                or any(not np.isfinite(v) or v <= 0.0
                       for v in self.spatial_fd_scales)):
            raise ValueError("spatial_fd_scales requires at least three positive scales")
        if (tuple(sorted(self.spatial_fd_scales)) != self.spatial_fd_scales
                or any(b <= a for a, b in zip(
                    self.spatial_fd_scales, self.spatial_fd_scales[1:]))):
            raise ValueError("spatial_fd_scales must be strictly increasing")
        if (not np.isfinite(self.spatial_fd_stability_rtol)
                or self.spatial_fd_stability_rtol <= 0.0):
            raise ValueError("spatial_fd_stability_rtol must be positive")
        if (not np.isfinite(self.spatial_fd_stability_atol_n2_per_km)
                or self.spatial_fd_stability_atol_n2_per_km <= 0.0):
            raise ValueError("spatial_fd_stability_atol_n2_per_km must be positive")

        self.provider_metadata = metadata
        self.provider_hash = str(metadata.get(
            "provider_hash", metadata.get("grid_hash", _metadata_hash(metadata))))
        self.fd_evaluation_count = 0
        self.fd_rejection_count = 0
        self.fd_relative_spread_max = 0.0
        self.fd_absolute_spread_max = 0.0
        self.fd_stability_ratio_max = 0.0
        self.fd_selected_platform_counts = {}

    def _Ne(self, ecef_km, t: float = 0.0) -> float:
        point = np.asarray(ecef_km, dtype=float)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("ECEF position must be a finite length-3 vector")
        value = float(self.provider.get_Ne_ecef(
            point, self.earth_radius_km, t=float(t)))
        if not np.isfinite(value) or value < 0.0:
            raise ValueError("Provider returned non-finite or negative Ne")
        return value

    def _magnetic_field(self, ecef_km, t: float = 0.0):
        field = np.asarray(
            self.geomagnetic_provider.get_B_ecef_T(ecef_km, t=t), dtype=float)
        if field.shape != (3,) or not np.all(np.isfinite(field)):
            raise ValueError("Geomagnetic Provider returned an invalid ECEF field")
        magnitude = float(np.linalg.norm(field))
        if magnitude <= 1.0e-18:
            return 0.0, self.b_ecef_hat
        return GYRO_MHZ_PER_T * magnitude, field / magnitude

    def coordinates_from_ecef(self, ecef_km) -> np.ndarray:
        if self.coordinate_model == WGS84_COORDINATE_MODEL:
            return ecef_to_geodetic(ecef_km)
        return ecef_to_geocentric(ecef_km, self.earth_radius_km)

    def position_to_ecef(self, lat_deg, lon_deg, height_km) -> np.ndarray:
        if self.coordinate_model == WGS84_COORDINATE_MODEL:
            return geodetic_to_ecef(lat_deg, lon_deg, height_km)
        return geocentric_to_ecef(
            lat_deg, lon_deg, height_km, self.earth_radius_km)

    def launch_direction(self, lat_deg, lon_deg, elevation_deg,
                         bearing_deg) -> np.ndarray:
        if self.coordinate_model == WGS84_COORDINATE_MODEL:
            return launch_direction_ecef_wgs84(
                lat_deg, lon_deg, elevation_deg, bearing_deg)
        return launch_direction_ecef(
            lat_deg, lon_deg, elevation_deg, bearing_deg,
            self.earth_radius_km)

    def physical_height_km(self, ecef_km) -> float:
        return float(self.coordinates_from_ecef(ecef_km)[2])

    def provider_domain_margin(self, ecef_km) -> float:
        return float(self.provider.domain_margin_ecef(
            ecef_km, self.earth_radius_km))

    def roots(self, ecef_km, k_ecef, *, omega=None, t: float = 0.0):
        frequency = self.omega if omega is None else omega
        gyro_MHz, b_ecef_hat = self._magnetic_field(ecef_km, t=t)
        return appleton_hartree_roots(
            frequency, self._Ne(ecef_km, t=t), k_ecef,
            b_ecef_hat, gyro_MHz, self.collision_s,
            plasma_coefficient=self.plasma_coefficient,
        )

    def root(self, ecef_km, k_ecef, *, omega=None, t: float = 0.0) -> AHRoot:
        return select_root(
            self.roots(ecef_km, k_ecef, omega=omega, t=t),
            self.wave_mode, near_cutoff_n2=self.near_cutoff_n2,
        )

    def n2_raw(self, ecef_km, k_ecef, *, omega=None, t: float = 0.0) -> float:
        root = self.root(ecef_km, k_ecef, omega=omega, t=t)
        if not np.isfinite(root.n2.real) or abs(root.n2.imag) > 1.0e-10:
            raise ValueError("AH root is not a finite real propagation root")
        return float(root.n2.real)

    def metadata(self) -> dict:
        return {
            "plasma_coefficient": self.plasma_coefficient,
            "compatibility_mode": self.compatibility_mode,
            "compatibility_scope": "continuous_radial_B0_collisionless" if self.compatibility_mode else None,
            "coordinate_model": self.coordinate_model,
            "state_semantics": STATE_SEMANTICS,
            "earth_model": (
                "WGS84" if self.coordinate_model == WGS84_COORDINATE_MODEL
                else "sphere"),
            "earth_radius_km": self.earth_radius_km,
            "frequency_MHz": self.freq_MHz,
            "branch_id": self.wave_mode,
            "gyrofrequency_MHz": self.gyro_MHz,
            "b_field_coordinates": self.geomagnetic_metadata.get(
                "coordinate_system", "ecef"),
            "b_ecef_hat": self.b_ecef_hat.tolist(),
            "geomagnetic_provider": self.geomagnetic_metadata,
            "geomagnetic_provider_hash": self.geomagnetic_provider_hash,
            "collision_frequency_s^-1": self.collision_s,
            "hamiltonian": "D=|k_ecef|^2-(omega/c)^2*n_AH^2",
            "spatial_derivative": "multiscale_centered_ECEF",
            "spatial_fd_base_km": self.spatial_fd_base_km,
            "spatial_fd_scales": list(self.spatial_fd_scales),
            "spatial_fd_stability_rtol": self.spatial_fd_stability_rtol,
            "spatial_fd_stability_atol_n2_per_km": (
                self.spatial_fd_stability_atol_n2_per_km),
            "provider_hash": self.provider_hash,
            "coordinate_mapping_hash": (
                _metadata_hash({
                    "coordinate_model": self.coordinate_model,
                    "mapping": "analytic_ECEF",
                }) if self.coordinate_model == WGS84_COORDINATE_MODEL
                else _COORDINATE_MAPPING_HASH),
        }


def hamiltonian_ah3d_ecef(state, model: RefractiveIndexAH3DECEF,
                          *, omega=None, t: float = 0.0) -> float:
    """Return D=|k|^2-(omega/c)^2*n_AH^2 in one ECEF basis."""
    value = np.asarray(state, dtype=float)
    if value.ndim != 1 or value.size not in (6, 7) or not np.all(np.isfinite(value)):
        raise ValueError("state must be a finite ECEF vector of length 6 or 7")
    angular_frequency = model.omega if omega is None else omega
    n2 = model.n2_raw(value[:3], value[3:6], omega=angular_frequency, t=t)
    result = dispersion_D(value[3:6], angular_frequency, n2)
    if not np.isfinite(result.real) or abs(result.imag) > 1.0e-10:
        raise ValueError("Hamiltonian is not finite and real")
    return float(result.real)


def hamiltonian_residual_ah3d_ecef(state, model, *, t: float = 0.0) -> float:
    return abs(hamiltonian_ah3d_ecef(state, model, t=t)) / max(
        model.k0 * model.k0, 1.0e-30)


def spatial_gradient_ah3d_ecef(state, model: RefractiveIndexAH3DECEF,
                               *, t: float = 0.0):
    """Resolve D_r on a 0.5/1/2-scale centered-difference platform."""
    value = np.asarray(state, dtype=float)
    estimates = []
    try:
        for scale in model.spatial_fd_scales:
            step = model.spatial_fd_base_km * scale
            gradient = np.empty(3, dtype=float)
            for axis in range(3):
                plus = value.copy()
                minus = value.copy()
                plus[axis] += step
                minus[axis] -= step
                gradient[axis] = (
                    hamiltonian_ah3d_ecef(plus, model, t=t)
                    - hamiltonian_ah3d_ecef(minus, model, t=t)
                ) / (2.0 * step)
            estimates.append(gradient)
    except (ValueError, FloatingPointError) as exc:
        model.fd_rejection_count += 1
        raise SpatialDerivativeError(
            f"ECEF spatial difference crossed an invalid Provider/dispersion domain: {exc}"
        ) from exc
    model.fd_evaluation_count += 1
    estimates = np.asarray(estimates)
    if not np.all(np.isfinite(estimates)):
        model.fd_rejection_count += 1
        raise SpatialDerivativeError("non-finite ECEF spatial derivative")
    platforms = []
    for start in range(len(estimates) - 2):
        candidate = estimates[start:start + 3]
        reference = np.median(candidate, axis=0)
        spread = float(np.max(np.linalg.norm(candidate - reference, axis=1)))
        magnitude = max(float(np.linalg.norm(reference)), model.k0**2 * 1.0e-12)
        platforms.append((spread / magnitude, start, reference))
    relative_spread, selected_start, reference = min(
        platforms, key=lambda item: item[0])
    selected_steps = tuple(model.spatial_fd_base_km * s for s in
                           model.spatial_fd_scales[selected_start:selected_start + 3])
    model.fd_relative_spread_max = max(
        model.fd_relative_spread_max, float(relative_spread))
    selected_spread = float(platforms[selected_start][0]
                            * max(np.linalg.norm(reference), model.k0**2 * 1.0e-12))
    model.fd_absolute_spread_max = max(
        model.fd_absolute_spread_max, selected_spread)
    stability_limit = (
        model.spatial_fd_stability_rtol * float(np.linalg.norm(reference))
        + model.k0**2 * model.spatial_fd_stability_atol_n2_per_km)
    stability_ratio = selected_spread / stability_limit
    model.fd_stability_ratio_max = max(
        model.fd_stability_ratio_max, float(stability_ratio))
    model.fd_selected_platform_counts[selected_steps] = (
        model.fd_selected_platform_counts.get(selected_steps, 0) + 1)
    if stability_ratio > 1.0:
        model.fd_rejection_count += 1
        raise SpatialDerivativeError(
            "no stable multiscale ECEF spatial-difference platform: "
            f"relative_spread={relative_spread:.6g}, "
            f"absolute_spread={selected_spread:.6g}, "
            f"combined_limit={stability_limit:.6g}")
    return reference, {
        "steps_km": [model.spatial_fd_base_km * s
                     for s in model.spatial_fd_scales],
        "selected_steps_km": list(selected_steps),
        "relative_spread": relative_spread,
        "absolute_spread": selected_spread,
        "stability_ratio": stability_ratio,
        "estimates": estimates,
    }


def dispersion_gradients_ah3d_ecef(state, model: RefractiveIndexAH3DECEF,
                                    *, t: float = 0.0):
    value = np.asarray(state, dtype=float)
    grad_r, diagnostics = spatial_gradient_ah3d_ecef(value, model, t=t)
    position = value[:3].copy()

    def root_fn(k_ecef, omega):
        n2 = model.roots(position, k_ecef, omega=omega, t=t)[model.wave_mode]
        return dispersion_D(k_ecef, omega, n2)

    grad_k, grad_omega = dispersion_partials(root_fn, value[3:6], model.omega)
    if (not np.all(np.isfinite(grad_k)) or not np.isfinite(grad_omega)
            or abs(grad_omega) < 1.0e-20):
        raise ValueError("non-finite or singular AH dispersion derivative")
    return grad_r, grad_k, grad_omega, diagnostics


def group_velocity_ah3d_ecef(state, model: RefractiveIndexAH3DECEF,
                              *, t: float = 0.0) -> np.ndarray:
    value = np.asarray(state, dtype=float)
    position = value[:3].copy()

    def root_fn(k_ecef, omega):
        n2 = model.roots(position, k_ecef, omega=omega, t=t)[model.wave_mode]
        return dispersion_D(k_ecef, omega, n2)

    grad_k, grad_omega = dispersion_partials(root_fn, value[3:6], model.omega)
    if not np.isfinite(grad_omega) or abs(grad_omega) < 1.0e-20:
        raise ValueError("singular D_omega prevents group-velocity evaluation")
    velocity = -grad_k / grad_omega
    if not np.all(np.isfinite(velocity)):
        raise ValueError("non-finite group velocity")
    return velocity


def ah3d_ecef_rhs(path_km: float, state, model: RefractiveIndexAH3DECEF,
                   *, medium_time_s: float = 0.0) -> np.ndarray:
    """Hamilton equations normalized to ECEF physical geometric arclength."""
    del path_km
    root = model.root(state[:3], state[3:6], t=medium_time_s)
    if root.status not in {"ok", "near_cutoff"}:
        raise ValueError(f"non-propagating AH branch: {root.status}")
    grad_r, grad_k, grad_omega, _ = dispersion_gradients_ah3d_ecef(
        state, model, t=medium_time_s)
    norm_grad_k = float(np.linalg.norm(grad_k))
    if not np.isfinite(norm_grad_k) or norm_grad_k <= 1.0e-20:
        raise ValueError("singular D_k prevents geometric-path normalization")
    velocity = -grad_k / grad_omega
    speed = float(np.linalg.norm(velocity))
    if not np.isfinite(speed) or speed <= 0.0:
        raise ValueError("non-positive AH group speed")
    return np.r_[grad_k / norm_grad_k, -grad_r / norm_grad_k, 1.0 / speed]


def project_dispersion_ah3d_ecef(state, model: RefractiveIndexAH3DECEF,
                                  *, t: float = 0.0, tolerance: float = 1.0e-10,
                                  max_iterations: int = 6):
    """Apply the minimum k-space correction along D_k to restore D=0."""
    projected = np.asarray(state, dtype=float).copy()
    correction_total = 0.0
    for _ in range(int(max_iterations)):
        residual = hamiltonian_ah3d_ecef(projected, model, t=t)
        if abs(residual) / max(model.k0**2, 1.0e-30) <= tolerance:
            return projected, correction_total, True
        position = projected[:3].copy()

        def root_fn(k_ecef, omega):
            n2 = model.roots(position, k_ecef, omega=omega, t=t)[model.wave_mode]
            return dispersion_D(k_ecef, omega, n2)

        grad_k, _ = dispersion_partials(root_fn, projected[3:6], model.omega)
        norm2 = float(np.dot(grad_k, grad_k))
        if not np.isfinite(norm2) or norm2 <= 1.0e-24:
            break
        correction = -residual * grad_k / norm2
        projected[3:6] += correction
        correction_total += float(np.linalg.norm(correction))
    return projected, correction_total, False


def _integration_profile(name: str) -> dict:
    profiles = UPDATE2.get("AH3D_INTEGRATION", {})
    if name not in profiles:
        raise ValueError(f"unknown AH3D integration profile: {name}")
    return dict(profiles[name])


def launch_state_ah3d_ecef(
    model: RefractiveIndexAH3DECEF,
    *,
    tx_lat_deg: float,
    tx_lon_deg: float,
    tx_height_km: float,
    elevation_deg: float,
    bearing_deg: float,
    medium_time_s: float = 0.0,
) -> np.ndarray:
    position = model.position_to_ecef(
        tx_lat_deg, tx_lon_deg, tx_height_km)
    direction = model.launch_direction(
        tx_lat_deg, tx_lon_deg, elevation_deg, bearing_deg)
    root = model.root(position, direction, t=medium_time_s)
    if root.status not in {"ok", "near_cutoff"} or root.n2.real <= 0.0:
        raise ValueError(f"invalid launch dispersion: {root.status}")
    state = np.r_[position, model.k0 * np.sqrt(root.n2.real) * direction, 0.0]
    if model.compatibility_mode:
        residual = abs(np.dot(state[3:6], state[3:6]) / model.k0**2
                       - model.n2_raw(position, state[3:6], t=medium_time_s))
        if not np.isfinite(residual) or residual > 1.0e-10:
            raise ValueError("initial normalized dispersion residual exceeds 1e-10")
    return state


def _physical_height_km(state, model: RefractiveIndexAH3DECEF) -> float:
    return model.physical_height_km(np.asarray(state)[:3])


def _radial_group_velocity(state, model, medium_time_s: float) -> float:
    position = np.asarray(state)[:3]
    if model.coordinate_model == WGS84_COORDINATE_MODEL:
        lat, lon, _ = ecef_to_geodetic(position)
        vertical = enu_basis_wgs84(float(lat), float(lon))[:, 2]
    else:
        vertical = position / np.linalg.norm(position)
    return float(np.dot(
        group_velocity_ah3d_ecef(state, model, t=medium_time_s),
        vertical,
    ))


def _locate_dense_event(dense, lo: float, hi: float, event_value,
                        *, tolerance_km: float = 1.0e-9):
    flo = float(event_value(dense(lo)))
    fhi = float(event_value(dense(hi)))
    if not (np.isfinite(flo) or np.isinf(flo)) or not np.isfinite(fhi):
        raise ValueError("event bracket is non-finite")
    if flo < 0.0 or fhi > 0.0:
        raise ValueError("event bracket does not cross from inside/above to outside/below")
    while hi - lo > tolerance_km:
        mid = 0.5 * (lo + hi)
        fmid = float(event_value(dense(mid)))
        if not (np.isfinite(fmid) or np.isinf(fmid)):
            raise ValueError("event function became non-finite")
        if fmid >= 0.0:
            lo = mid
        else:
            hi = mid
    path = 0.5 * (lo + hi)
    return path, np.asarray(dense(path), dtype=float), hi - lo


def integrate_fixed_path_ah3d_ecef(
    model: RefractiveIndexAH3DECEF,
    *,
    tx_lat_deg: float,
    tx_lon_deg: float,
    tx_height_km: float,
    elevation_deg: float,
    bearing_deg: float,
    path_length_km: float,
    medium_time_s: float = 0.0,
    profile: str | None = None,
    h_stop_km: float | None = None,
    stop_on_physical_events: bool = False,
    apply_dispersion_projection: bool = True,
) -> dict:
    """Integrate physical arclength, optionally stopping on C0R-4 events."""
    length = float(path_length_km)
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("path_length_km must be finite and positive")
    profile_name = str(profile or model.integration_profile)
    settings = _integration_profile(profile_name)
    initial = launch_state_ah3d_ecef(
        model, tx_lat_deg=tx_lat_deg, tx_lon_deg=tx_lon_deg,
        tx_height_km=tx_height_km, elevation_deg=elevation_deg,
        bearing_deg=bearing_deg, medium_time_s=medium_time_s)
    departure_tolerance_km = 1.0e-6
    if stop_on_physical_events:
        if h_stop_km is None or not np.isfinite(float(h_stop_km)):
            raise ValueError("h_stop_km must be finite for physical-event tracing")
        h_stop_km = float(h_stop_km)
        if h_stop_km <= float(tx_height_km):
            raise ValueError("h_stop_km must exceed TX physical height")
        if model.provider_domain_margin(initial[:3]) < 0.0:
            raise ValueError("TX is outside the Provider domain")
        initial_radial_velocity = _radial_group_velocity(
            initial, model, medium_time_s)
        if (float(tx_height_km) <= departure_tolerance_km
                and initial_radial_velocity <= 0.0):
            raise ValueError("surface launch must initially depart outward")
    atol = np.array(
        [settings["atol_r_km"]] * 3 + [settings["atol_k"]] * 3 + [1.0e-12])

    def rhs(path, state):
        return ah3d_ecef_rhs(path, state, model, medium_time_s=medium_time_s)

    first_step = min(float(settings["ds_min_km"]), length)
    solver = DOP853(
        rhs, 0.0, initial, length, rtol=float(settings["rtol"]), atol=atol,
        first_step=first_step, max_step=float(settings["ds_max_km"]))
    states = [initial.copy()]
    paths = [0.0]
    residual_max = hamiltonian_residual_ah3d_ecef(initial, model, t=medium_time_s)
    projection_max = 0.0
    accepted = 0
    termination = ("max_path" if stop_on_physical_events
                   else "fixed_geometric_path_completed")
    failure = None
    event_records = []
    event_error_max = 0.0
    turning_count = 0
    apex_state = None
    apex_path = None
    left_ground = float(tx_height_km) > departure_tolerance_km
    previous_height = _physical_height_km(initial, model)
    previous_margin = model.provider_domain_margin(initial[:3])
    previous_radial_velocity = _radial_group_velocity(
        initial, model, medium_time_s)
    try:
        while solver.status == "running":
            path_previous = float(solver.t)
            solver.step()
            if solver.status == "failed" or not np.all(np.isfinite(solver.y)):
                termination = "integration_failed"
                failure = "DOP853 failed or returned a non-finite state"
                break
            dense = solver.dense_output()
            if apply_dispersion_projection:
                projected, correction, ok = project_dispersion_ah3d_ecef(
                    solver.y, model, t=medium_time_s)
                projection_max = max(projection_max, correction)
                if not ok:
                    termination = "constraint_failed"
                    failure = "D=0 projection did not converge"
                    break
                solver.y[:] = projected
                solver.f = rhs(solver.t, solver.y)
            accepted += 1
            residual_max = max(
                residual_max,
                hamiltonian_residual_ah3d_ecef(
                    solver.y, model, t=medium_time_s),
            )
            if not stop_on_physical_events:
                states.append(solver.y.copy())
                paths.append(float(solver.t))
                continue

            current = solver.y.copy()
            current_path = float(solver.t)
            current_height = _physical_height_km(current, model)
            current_margin = model.provider_domain_margin(current[:3])
            current_radial_velocity = _radial_group_velocity(
                current, model, medium_time_s)
            terminal_candidates = []

            if (left_ground and previous_height > 0.0
                    and current_height <= 0.0
                    and current_radial_velocity < 0.0):
                terminal_candidates.append((
                    "ground",
                    _locate_dense_event(
                        dense, path_previous, current_path,
                        lambda state: _physical_height_km(state, model)),
                ))
            upper_previous = float(h_stop_km) - previous_height
            upper_current = float(h_stop_km) - current_height
            if upper_previous >= 0.0 and upper_current < 0.0:
                terminal_candidates.append((
                    "upper_boundary",
                    _locate_dense_event(
                        dense, path_previous, current_path,
                        lambda state: float(h_stop_km)
                        - _physical_height_km(state, model)),
                ))
            if (np.isfinite(previous_margin) and previous_margin >= 0.0
                    and current_margin < 0.0):
                terminal_candidates.append((
                    "grid_exit",
                    _locate_dense_event(
                        dense, path_previous, current_path,
                        lambda state: model.provider_domain_margin(state[:3])),
                ))

            turning_candidate = None
            if previous_radial_velocity > 0.0 and current_radial_velocity <= 0.0:
                turning_candidate = _locate_dense_event(
                    dense, path_previous, current_path,
                    lambda state: _radial_group_velocity(
                        state, model, medium_time_s))

            terminal = (min(terminal_candidates, key=lambda item: item[1][0])
                        if terminal_candidates else None)
            terminal_path = float("inf") if terminal is None else terminal[1][0]
            if turning_candidate is not None and turning_candidate[0] < terminal_path:
                turn_path, turn_state, turn_error = turning_candidate
                if apply_dispersion_projection:
                    turn_state, _, turn_ok = project_dispersion_ah3d_ecef(
                        turn_state, model, t=medium_time_s)
                    if not turn_ok:
                        termination = "constraint_failed"
                        failure = "D=0 projection failed at turning event"
                        break
                turning_count += 1
                apex_state = turn_state.copy()
                apex_path = float(turn_path)
                event_error_max = max(event_error_max, float(turn_error))
                event_records.append({
                    "event": "turning_apex",
                    "geom_path_km": float(turn_path),
                    "localization_bracket_km": float(turn_error),
                    "direction": "radial_group_velocity_positive_to_negative",
                })
                states.append(turn_state.copy())
                paths.append(float(turn_path))

            if terminal is not None:
                event_name, (event_path, event_state, event_error) = terminal
                if apply_dispersion_projection:
                    event_state, _, event_ok = project_dispersion_ah3d_ecef(
                        event_state, model, t=medium_time_s)
                    if not event_ok:
                        termination = "constraint_failed"
                        failure = f"D=0 projection failed at {event_name} event"
                        break
                event_error_max = max(event_error_max, float(event_error))
                event_records.append({
                    "event": event_name,
                    "geom_path_km": float(event_path),
                    "localization_bracket_km": float(event_error),
                    "direction": {
                        "ground": "descending",
                        "upper_boundary": "outward",
                        "grid_exit": "provider_domain_exit",
                    }[event_name],
                })
                states.append(event_state.copy())
                paths.append(float(event_path))
                residual_max = max(
                    residual_max,
                    hamiltonian_residual_ah3d_ecef(
                        event_state, model, t=medium_time_s),
                )
                termination = event_name
                break

            states.append(current)
            paths.append(current_path)
            left_ground = left_ground or current_height > departure_tolerance_km
            previous_height = current_height
            previous_margin = current_margin
            previous_radial_velocity = current_radial_velocity
    except (ValueError, SpatialDerivativeError) as exc:
        termination = "physics_domain_rejected"
        failure = str(exc)

    trajectory = np.asarray(states)
    coordinates = np.asarray([
        model.coordinates_from_ecef(row[:3]) for row in trajectory
    ])
    completed_path = float(paths[-1])
    group_path = float(C_KM_S * trajectory[-1, 6])
    positions = trajectory[:, :3]
    tangents = np.empty_like(positions)
    if len(positions) == 1:
        tangents[0] = model.launch_direction(
            tx_lat_deg, tx_lon_deg, elevation_deg, bearing_deg)
    else:
        tangents[0] = positions[1] - positions[0]
        tangents[-1] = positions[-1] - positions[-2]
        if len(positions) > 2:
            tangents[1:-1] = positions[2:] - positions[:-2]
        tangent_norms = np.linalg.norm(tangents, axis=1)
        tangents /= np.maximum(tangent_norms[:, None], 1.0e-30)
    metadata = model.metadata()
    landing = (coordinates[-1] if termination == "ground" else None)
    apex_coordinates = (
        None if apex_state is None
        else model.coordinates_from_ecef(apex_state[:3]))
    launch_gyro_MHz, launch_b_hat = model._magnetic_field(
        trajectory[0, :3], t=medium_time_s)
    magnetic_fields = []
    for point in trajectory[:, :3]:
        try:
            magnetic_fields.append(
                model.geomagnetic_provider.get_B_ecef_T(
                    point, t=medium_time_s))
        except ValueError:
            magnetic_fields.append(np.full(3, np.nan))
    magnetic_fields = np.asarray(magnetic_fields, dtype=float)
    coordinate_label = (
        "geodetic" if model.coordinate_model == WGS84_COORDINATE_MODEL
        else "geocentric")
    return {
        "trajectory_ecef_km": trajectory[:, :3],
        "plasma_coefficient": model.plasma_coefficient,
        "compatibility_mode": model.compatibility_mode,
        "compatibility_scope": metadata["compatibility_scope"],
        "trajectory_group_tangent_ecef": tangents,
        "trajectory_group_tangent_ecef_semantics": "normalized_position_chord_estimate",
        "trajectory_group_path_km": C_KM_S * trajectory[:, 6],
        "trajectory_group_time_s": trajectory[:, 6].copy(),
        "trajectory_group_path_semantics": "c_times_integrated_group_time_state",
        "wavevector_ecef_km^-1": trajectory[:, 3:6],
        "trajectory_B_ecef_T": magnetic_fields,
        "trajectory_gyrofrequency_MHz": (
            GYRO_MHZ_PER_T * np.linalg.norm(magnetic_fields, axis=1)),
        "trajectory_lat_lon_height": coordinates,
        "trajectory_geocentric_lat_lon_height": (
            coordinates if coordinate_label == "geocentric" else None),
        "trajectory_geodetic_lat_lon_height": (
            coordinates if coordinate_label == "geodetic" else None),
        "trajectory_physical_height_km": coordinates[:, 2],
        "trajectory_geom_path_km": np.asarray(paths),
        "geom_path_km": completed_path,
        "group_path_km": group_path,
        "tau_ms": float(1000.0 * trajectory[-1, 6]),
        "termination_status": termination,
        "event_history": event_records,
        "event_localization_error_max_km": float(event_error_max),
        "turning_point_count": turning_count,
        "landing_geocentric_lat_lon_height": (
            None if landing is None or coordinate_label != "geocentric"
            else landing.tolist()),
        "landing_geodetic_lat_lon_height": (
            None if landing is None or coordinate_label != "geodetic"
            else landing.tolist()),
        "landing_lat_lon_height": (
            None if landing is None else landing.tolist()),
        "landing_lat_deg": None if landing is None else float(landing[0]),
        "landing_lon_deg": None if landing is None else float(landing[1]),
        "landing_height_km": None if landing is None else float(landing[2]),
        "apex_ecef_km": None if apex_state is None else apex_state[:3].tolist(),
        "apex_wavevector_ecef_km^-1": (
            None if apex_state is None else apex_state[3:6].tolist()),
        "apex_geocentric_lat_lon_height": (
            None if apex_coordinates is None or coordinate_label != "geocentric"
            else apex_coordinates.tolist()),
        "apex_geodetic_lat_lon_height": (
            None if apex_coordinates is None or coordinate_label != "geodetic"
            else apex_coordinates.tolist()),
        "apex_lat_lon_height": (
            None if apex_coordinates is None else apex_coordinates.tolist()),
        "apex_height_km": (
            None if apex_coordinates is None else float(apex_coordinates[2])),
        "apex_geom_path_km": apex_path,
        "apex_radial_group_velocity_km_s": (
            None if apex_state is None else _radial_group_velocity(
                apex_state, model, medium_time_s)),
        "failure_reason": failure,
        "branch_id": model.wave_mode,
        "wave_mode": model.wave_mode,
        "frequency_MHz": model.freq_MHz,
        "gyrofrequency_MHz": float(launch_gyro_MHz),
        "b_ecef_hat": launch_b_hat.tolist(),
        "geomagnetic_field_semantics": "launch_value_only_local_field_used_along_path",
        "geomagnetic_provider_hash": model.geomagnetic_provider_hash,
        "tx_lat_lon_height": [
            float(tx_lat_deg), float(tx_lon_deg), float(tx_height_km)],
        "tx_geocentric_lat_lon_height": (
            [float(tx_lat_deg), float(tx_lon_deg), float(tx_height_km)]
            if coordinate_label == "geocentric" else None),
        "tx_geodetic_lat_lon_height": (
            [float(tx_lat_deg), float(tx_lon_deg), float(tx_height_km)]
            if coordinate_label == "geodetic" else None),
        "elevation_tx_deg": float(elevation_deg),
        "bearing_tx_deg": float(bearing_deg) % 360.0,
        "coordinate_model": model.coordinate_model,
        "state_semantics": STATE_SEMANTICS,
        "earth_model": (
            "WGS84" if model.coordinate_model == WGS84_COORDINATE_MODEL
            else "sphere"),
        "earth_radius_km": model.earth_radius_km,
        "integrator": "DOP853",
        "integration_profile": profile_name,
        "integration_scope": (
            ("wgs84_physical_events" if
             model.coordinate_model == WGS84_COORDINATE_MODEL
             else "spherical_physical_events_C0R4") if stop_on_physical_events
            else "fixed_physical_path_core_no_events_C0R3"),
        "accepted_steps": accepted,
        "nfev": int(solver.nfev),
        "hamiltonian_residual_max": float(residual_max),
        "dispersion_projection_correction_max": float(projection_max),
        "dispersion_projection_enabled": bool(apply_dispersion_projection),
        "group_velocity_definition": "v_g=-D_k/D_omega",
        "group_path_definition": "c*t_group",
        "provider_hash": model.provider_hash,
        "coordinate_mapping_hash": metadata["coordinate_mapping_hash"],
        "spatial_fd_evaluation_count": model.fd_evaluation_count,
        "spatial_fd_base_km": model.spatial_fd_base_km,
        "spatial_fd_scales": list(model.spatial_fd_scales),
        "spatial_fd_rejection_count": model.fd_rejection_count,
        "spatial_fd_relative_spread_max": float(model.fd_relative_spread_max),
        "spatial_fd_absolute_spread_max": float(model.fd_absolute_spread_max),
        "spatial_fd_stability_ratio_max": float(model.fd_stability_ratio_max),
        "spatial_fd_selected_platform_counts": {
            ",".join(f"{step:.12g}" for step in steps): int(count)
            for steps, count in model.fd_selected_platform_counts.items()
        },
        "theory_reference": (
            "Haselgrove (1955); Haselgrove and Haselgrove (1960); "
            "Jones and Stephenson (1975); Budden (1985)"),
    }


def trace_ray_ah3d_ecef(
    model: RefractiveIndexAH3DECEF,
    *,
    tx_lat_deg: float,
    tx_lon_deg: float,
    tx_height_km: float,
    elevation_deg: float,
    bearing_deg: float,
    h_stop_km: float,
    max_path_km: float,
    medium_time_s: float = 0.0,
    profile: str | None = None,
    apply_dispersion_projection: bool = True,
) -> dict:
    """Trace until ground, upper boundary, Provider grid exit, or max path."""
    return integrate_fixed_path_ah3d_ecef(
        model,
        tx_lat_deg=tx_lat_deg,
        tx_lon_deg=tx_lon_deg,
        tx_height_km=tx_height_km,
        elevation_deg=elevation_deg,
        bearing_deg=bearing_deg,
        path_length_km=max_path_km,
        medium_time_s=medium_time_s,
        profile=profile,
        h_stop_km=h_stop_km,
        stop_on_physical_events=True,
        apply_dispersion_projection=apply_dispersion_projection,
    )


__all__ = [
    "COORDINATE_MODEL",
    "STATE_SEMANTICS",
    "SpatialDerivativeError",
    "RefractiveIndexAH3DECEF",
    "hamiltonian_ah3d_ecef",
    "hamiltonian_residual_ah3d_ecef",
    "spatial_gradient_ah3d_ecef",
    "dispersion_gradients_ah3d_ecef",
    "group_velocity_ah3d_ecef",
    "ah3d_ecef_rhs",
    "project_dispersion_ah3d_ecef",
    "launch_state_ah3d_ecef",
    "integrate_fixed_path_ah3d_ecef",
    "trace_ray_ah3d_ecef",
]
