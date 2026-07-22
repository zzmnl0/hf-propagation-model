"""Strict spherical-geocentric electron-density providers for C0R and later."""

from __future__ import annotations

from abc import ABC, abstractmethod
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from utils import (ecef_to_geocentric, enu_basis, geocentric_to_ecef,
                   spherical_destination)
from .ionosphere_components_3d import (
    Bubble3DComponent,
    Es3DComponent,
    SpreadF3DComponent,
    TID3DComponent,
    geomag_unit_vector,
)


COORDINATE_SYSTEM = "geocentric_spherical"
AXIS_ORDER = "latitude,longitude,height"


def _hash_value(value) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("ascii")).hexdigest()


def _hash_arrays(metadata: dict, *arrays) -> str:
    digest = hashlib.sha256(json.dumps(
        metadata, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("ascii"))
    for value in arrays:
        arr = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest.update(str(arr.shape).encode("ascii"))
        digest.update(arr.tobytes())
    return digest.hexdigest()


def _scalar_text(value, name: str) -> str:
    arr = np.asarray(value)
    if arr.shape != ():
        raise ValueError(f"{name} must be a scalar string")
    return str(arr.item())


def _derive_link_coordinates(lat_deg, lon_deg, height_km, tx_lat_deg,
                             tx_lon_deg, link_bearing_deg, earth_radius_km):
    tx_up = geocentric_to_ecef(
        tx_lat_deg, tx_lon_deg, 0.0, earth_radius_km)
    tx_up /= np.linalg.norm(tx_up)
    point = geocentric_to_ecef(
        lat_deg, lon_deg, height_km, earth_radius_km)
    point_up = point / np.linalg.norm(point)
    dot = float(np.clip(np.dot(tx_up, point_up), -1.0, 1.0))
    central_angle = float(np.arctan2(np.linalg.norm(np.cross(tx_up, point_up)), dot))
    tangent = point_up - dot * tx_up
    tangent_norm = float(np.linalg.norm(tangent))
    relative_bearing = 0.0
    if tangent_norm > 1.0e-14:
        east, north = (tangent / tangent_norm) @ enu_basis(
            tx_lat_deg, tx_lon_deg, earth_radius_km)[:, :2]
        relative_bearing = (float(np.arctan2(east, north))
                            - np.radians(float(link_bearing_deg)))
    distance = float(earth_radius_km) * central_angle
    return np.array([distance * np.cos(relative_bearing),
                     distance * np.sin(relative_bearing), float(height_km)])


class GeocentricNeProvider3D(ABC):
    """Formal Ne interface: geocentric lat/lon/height or ECEF, never local xyz."""

    name = "geocentric_ne_provider_3d"

    @abstractmethod
    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        """Return Ne [m^-3] at geocentric latitude, east longitude, and height."""

    def get_Ne_ecef(self, ecef_km, earth_radius_km: float,
                    t: float = 0.0) -> float:
        provider_radius = getattr(self, "earth_radius_km", None)
        if (provider_radius is not None
                and float(earth_radius_km) != float(provider_radius)):
            raise ValueError("ECEF query Earth radius differs from Provider radius")
        lat, lon, height = ecef_to_geocentric(ecef_km, earth_radius_km)
        return self.get_Ne_geocentric(float(lat), float(lon), float(height), t=t)

    def domain_margin_geocentric(self, lat_deg: float, lon_deg: float,
                                  height_km: float) -> float:
        """Signed Provider-domain margin; positive is inside, inf is unbounded."""
        geocentric_to_ecef(
            lat_deg, lon_deg, height_km,
            float(getattr(self, "earth_radius_km")))
        return float("inf")

    def domain_margin_ecef(self, ecef_km, earth_radius_km: float) -> float:
        provider_radius = getattr(self, "earth_radius_km", None)
        if (provider_radius is None
                or float(earth_radius_km) != float(provider_radius)):
            raise ValueError("ECEF domain query Earth radius differs from Provider radius")
        lat, lon, height = ecef_to_geocentric(ecef_km, earth_radius_km)
        return self.domain_margin_geocentric(
            float(lat), float(lon), float(height))

    def get_grad_Ne_ecef(self, ecef_km, earth_radius_km: float,
                         t: float = 0.0, step_km: float = 1.0) -> np.ndarray:
        point = np.asarray(ecef_km, dtype=float)
        step = float(step_km)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("ecef_km must be a finite length-3 vector")
        if not np.isfinite(step) or step <= 0.0:
            raise ValueError("step_km must be finite and positive")
        gradient = np.empty(3, dtype=float)
        for axis in range(3):
            plus = point.copy()
            minus = point.copy()
            plus[axis] += step
            minus[axis] -= step
            gradient[axis] = (
                self.get_Ne_ecef(plus, earth_radius_km, t=t)
                - self.get_Ne_ecef(minus, earth_radius_km, t=t)
            ) / (2.0 * step)
        return gradient

    def sample_geocentric_grid(self, latitude_deg, longitude_deg, height_km,
                               t: float = 0.0) -> np.ndarray:
        lat = np.asarray(latitude_deg, dtype=float)
        lon = np.asarray(longitude_deg, dtype=float)
        height = np.asarray(height_km, dtype=float)
        if any(axis.ndim != 1 for axis in (lat, lon, height)):
            raise ValueError("latitude, longitude, and height axes must be 1-D")
        out = np.empty((lat.size, lon.size, height.size), dtype=float)
        for i, lat_value in enumerate(lat):
            for j, lon_value in enumerate(lon):
                for k, height_value in enumerate(height):
                    out[i, j, k] = self.get_Ne_geocentric(
                        float(lat_value), float(lon_value), float(height_value), t=t)
        return out

    def sample_link_plane(self, tx_geocentric, bearing_deg: float,
                          ground_range_km, height_km, t: float = 0.0) -> np.ndarray:
        """Sample a derived great-circle range/height plane for plotting only."""
        tx = np.asarray(tx_geocentric, dtype=float)
        ranges = np.asarray(ground_range_km, dtype=float)
        heights = np.asarray(height_km, dtype=float)
        if tx.shape != (3,) or ranges.ndim != 1 or heights.ndim != 1:
            raise ValueError("tx_geocentric must have 3 values and plot axes must be 1-D")
        out = np.empty((ranges.size, heights.size), dtype=float)
        for i, distance in enumerate(ranges):
            lat, lon = spherical_destination(
                float(tx[0]), float(tx[1]), float(bearing_deg),
                float(distance), float(self.earth_radius_km))
            for j, height in enumerate(heights):
                out[i, j] = self.get_Ne_geocentric(
                    lat, lon, float(height), t=t)
        return out

    def get_metadata(self) -> dict:
        return {
            "name": self.name,
            "dimensions": "3d",
            "units": "m^-3",
            "coordinate_system": COORDINATE_SYSTEM,
            "axis_order": AXIS_ORDER,
        }


class GeocentricGridNeProvider3D(GeocentricNeProvider3D):
    """Linear interpolation of Ne(latitude, longitude, height) on a sphere."""

    name = "geocentric_grid_ne_provider_3d"

    def __init__(self, Ne_m3, latitude_deg, longitude_deg, height_km,
                 earth_radius_km: float, *, fill_value_m3: float = 0.0,
                 longitude_periodic_complete: bool = False,
                 source_metadata: dict | None = None):
        self.latitude_deg = np.asarray(latitude_deg, dtype=float)
        self.longitude_deg = np.asarray(longitude_deg, dtype=float)
        self.height_km = np.asarray(height_km, dtype=float)
        self.earth_radius_km = float(earth_radius_km)
        self.fill_value_m3 = float(fill_value_m3)
        if not isinstance(longitude_periodic_complete, (bool, np.bool_)):
            raise ValueError("longitude_periodic_complete must be boolean")
        self.longitude_periodic_complete = bool(longitude_periodic_complete)
        self.Ne_m3 = np.asarray(Ne_m3, dtype=float)
        self.source_metadata = dict(source_metadata or {})
        self._validate()
        self._longitude_center_deg = 0.5 * (
            float(self.longitude_deg[0]) + float(self.longitude_deg[-1]))
        self.grid_hash = _hash_arrays({
            "coordinate_system": COORDINATE_SYSTEM,
            "axis_order": AXIS_ORDER,
            "units": "m^-3",
            "earth_radius_km": self.earth_radius_km,
            "fill_value_m3": self.fill_value_m3,
            "longitude_periodic_complete": self.longitude_periodic_complete,
        }, self.latitude_deg, self.longitude_deg, self.height_km, self.Ne_m3)
        interpolation_longitude = self.longitude_deg
        interpolation_Ne = self.Ne_m3
        if self.longitude_periodic_complete:
            interpolation_longitude = np.r_[
                self.longitude_deg[-1] - 360.0,
                self.longitude_deg,
                self.longitude_deg[0] + 360.0,
            ]
            interpolation_Ne = np.concatenate((
                self.Ne_m3[:, -1:, :], self.Ne_m3, self.Ne_m3[:, :1, :]),
                axis=1,
            )
        self._interp = RegularGridInterpolator(
            (self.latitude_deg, interpolation_longitude, self.height_km),
            interpolation_Ne,
            method="linear", bounds_error=False, fill_value=self.fill_value_m3)

    def _validate(self) -> None:
        if not np.isfinite(self.earth_radius_km) or self.earth_radius_km <= 0.0:
            raise ValueError("earth_radius_km must be finite and positive")
        if not np.isfinite(self.fill_value_m3) or self.fill_value_m3 < 0.0:
            raise ValueError("fill_value_m3 must be finite and nonnegative")
        axes = (
            ("latitude_deg", self.latitude_deg),
            ("longitude_deg", self.longitude_deg),
            ("height_km", self.height_km),
        )
        for name, axis in axes:
            if (axis.ndim != 1 or axis.size < 2 or not np.all(np.isfinite(axis))
                    or not np.all(np.diff(axis) > 0.0)):
                raise ValueError(f"{name} must be finite, 1-D, and strictly increasing")
        if self.latitude_deg[0] < -90.0 or self.latitude_deg[-1] > 90.0:
            raise ValueError("latitude_deg must be in [-90, 90]")
        if self.longitude_deg[-1] - self.longitude_deg[0] >= 360.0:
            raise ValueError("longitude axis must not duplicate a periodic endpoint")
        if self.longitude_periodic_complete:
            wrap_gap = self.longitude_deg[0] + 360.0 - self.longitude_deg[-1]
            if wrap_gap <= 0.0 or wrap_gap > 2.0 * np.max(np.diff(self.longitude_deg)):
                raise ValueError(
                    "complete periodic longitude grid has an unresolved wrap gap")
        if np.any(self.earth_radius_km + self.height_km <= 0.0):
            raise ValueError("height axis must remain above the Earth center")
        expected = (self.latitude_deg.size, self.longitude_deg.size,
                    self.height_km.size)
        if self.Ne_m3.shape != expected:
            raise ValueError(
                f"Ne shape {self.Ne_m3.shape} does not match axis order "
                f"{AXIS_ORDER} with shape {expected}")
        if not np.all(np.isfinite(self.Ne_m3)) or np.any(self.Ne_m3 < 0.0):
            raise ValueError("Ne_m3 must be finite and nonnegative")

    def _periodic_longitude(self, lon_deg: float) -> float:
        lon = float(lon_deg)
        if not np.isfinite(lon):
            raise ValueError("longitude must be finite")
        return self._longitude_center_deg + (
            (lon - self._longitude_center_deg + 180.0) % 360.0 - 180.0)

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        geocentric_to_ecef(lat_deg, lon_deg, height_km, self.earth_radius_km)
        point = [[float(lat_deg), self._periodic_longitude(lon_deg),
                  float(height_km)]]
        return max(float(self._interp(point)[0]), 0.0)

    def domain_margin_geocentric(self, lat_deg: float, lon_deg: float,
                                  height_km: float) -> float:
        geocentric_to_ecef(
            lat_deg, lon_deg, height_km, self.earth_radius_km)
        latitude = float(lat_deg)
        height = float(height_km)
        margins = [
            latitude - float(self.latitude_deg[0]),
            float(self.latitude_deg[-1]) - latitude,
            height - float(self.height_km[0]),
            float(self.height_km[-1]) - height,
        ]
        if not self.longitude_periodic_complete:
            longitude = self._periodic_longitude(lon_deg)
            margins.extend((
                longitude - float(self.longitude_deg[0]),
                float(self.longitude_deg[-1]) - longitude,
            ))
        return float(min(margins))

    def get_metadata(self) -> dict:
        return {
            **super().get_metadata(),
            "earth_radius_km": self.earth_radius_km,
            "latitude_range_deg": [float(self.latitude_deg[0]),
                                   float(self.latitude_deg[-1])],
            "longitude_range_deg_unwrapped": [float(self.longitude_deg[0]),
                                               float(self.longitude_deg[-1])],
            "height_range_km": [float(self.height_km[0]),
                                float(self.height_km[-1])],
            "grid_shape": list(self.Ne_m3.shape),
            "interpolation": "linear",
            "bounds_error": False,
            "fill_value_m^-3": self.fill_value_m3,
            "longitude_semantics": (
                "east_positive_periodic_full_interpolation"
                if self.longitude_periodic_complete
                else "east_positive_periodic_query_regional_domain"),
            "longitude_periodic_complete": self.longitude_periodic_complete,
            "grid_hash": self.grid_hash,
            "source_metadata": self.source_metadata,
        }


class GeocentricLinkGridNeProvider3D(GeocentricNeProvider3D):
    """Range-height background queried through formal spherical coordinates."""

    name = "geocentric_link_grid_ne_provider_3d"

    def __init__(self, Ne_m3, ground_range_km, height_km, *,
                 tx_lat_deg: float, tx_lon_deg: float,
                 link_bearing_deg: float, earth_radius_km: float,
                 fill_value_m3: float = 0.0, source_metadata: dict | None = None):
        self.Ne_m3 = np.asarray(Ne_m3, dtype=float)
        self.ground_range_km = np.asarray(ground_range_km, dtype=float)
        self.height_km = np.asarray(height_km, dtype=float)
        self.tx_lat_deg = float(tx_lat_deg)
        self.tx_lon_deg = float(tx_lon_deg)
        self.link_bearing_deg = float(link_bearing_deg) % 360.0
        self.earth_radius_km = float(earth_radius_km)
        self.fill_value_m3 = float(fill_value_m3)
        self.source_metadata = dict(source_metadata or {})
        if (self.ground_range_km.ndim != 1 or self.height_km.ndim != 1
                or self.Ne_m3.shape != (self.ground_range_km.size, self.height_km.size)
                or not np.all(np.diff(self.ground_range_km) > 0.0)
                or not np.all(np.diff(self.height_km) > 0.0)
                or not np.all(np.isfinite(self.Ne_m3)) or np.any(self.Ne_m3 < 0.0)):
            raise ValueError("invalid range-height background grid")
        geocentric_to_ecef(
            self.tx_lat_deg, self.tx_lon_deg, 0.0, self.earth_radius_km)
        self._interp = RegularGridInterpolator(
            (self.ground_range_km, self.height_km), self.Ne_m3,
            method="linear", bounds_error=False, fill_value=self.fill_value_m3)
        self.grid_hash = _hash_arrays({
            "coordinate_system": COORDINATE_SYSTEM,
            "derived_axes": "great_circle_ground_range,physical_height",
            "earth_radius_km": self.earth_radius_km,
        }, self.ground_range_km, self.height_km, self.Ne_m3)

    def derive_link_coordinates(self, lat_deg: float, lon_deg: float,
                                height_km: float) -> np.ndarray:
        return _derive_link_coordinates(
            lat_deg, lon_deg, height_km, self.tx_lat_deg, self.tx_lon_deg,
            self.link_bearing_deg, self.earth_radius_km)

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        along, _, height = self.derive_link_coordinates(
            lat_deg, lon_deg, height_km)
        return max(float(self._interp([[along, height]])[0]), 0.0)

    def domain_margin_geocentric(self, lat_deg: float, lon_deg: float,
                                  height_km: float) -> float:
        along, _, height = self.derive_link_coordinates(
            lat_deg, lon_deg, height_km)
        return float(min(
            along - self.ground_range_km[0], self.ground_range_km[-1] - along,
            height - self.height_km[0], self.height_km[-1] - height))

    def get_metadata(self) -> dict:
        frozen = {
            "earth_radius_km": self.earth_radius_km,
            "tx_lat_deg": self.tx_lat_deg,
            "tx_lon_deg": self.tx_lon_deg,
            "link_bearing_deg": self.link_bearing_deg,
            "background_semantics": "derived_link_range_height_y_invariant",
            "grid_hash": self.grid_hash,
            "source_metadata": self.source_metadata,
        }
        return {**super().get_metadata(), **frozen,
                "provider_hash": _hash_value(frozen)}


class ExternalGeocentricNeProvider3D(GeocentricGridNeProvider3D):
    """Strict NPZ reader; required metadata prevents axis/unit inference."""

    name = "external_geocentric_ne_provider_3d"
    REQUIRED = {
        "Ne", "latitude_deg", "longitude_deg", "height_km", "units",
        "coordinate_system", "axis_order", "earth_radius_km",
        "longitude_periodic_complete",
    }

    def __init__(self, path: str | Path, *, fill_value_m3: float = 0.0):
        source = Path(path)
        with np.load(source, allow_pickle=False) as data:
            missing = sorted(self.REQUIRED - set(data.files))
            if missing:
                raise ValueError("external Ne file missing required fields: "
                                 + ", ".join(missing))
            coordinate_system = _scalar_text(data["coordinate_system"],
                                             "coordinate_system")
            axis_order = _scalar_text(data["axis_order"], "axis_order")
            units = _scalar_text(data["units"], "units")
            if coordinate_system != COORDINATE_SYSTEM:
                raise ValueError(f"coordinate_system must be {COORDINATE_SYSTEM}")
            if axis_order != AXIS_ORDER:
                raise ValueError(f"axis_order must be {AXIS_ORDER}")
            Ne = np.asarray(data["Ne"], dtype=float)
            unit_key = units.strip().lower().replace(" ", "")
            if unit_key in {"m^-3", "m-3", "1/m3", "m**-3"}:
                Ne_m3 = Ne
                conversion = "none"
            elif unit_key in {"cm^-3", "cm-3", "1/cm3", "cm**-3"}:
                Ne_m3 = Ne * 1.0e6
                conversion = "cm^-3_to_m^-3_x1e6"
            else:
                raise ValueError(f"unsupported Ne units: {units}")
            latitude = np.asarray(data["latitude_deg"], dtype=float).copy()
            longitude = np.asarray(data["longitude_deg"], dtype=float).copy()
            height = np.asarray(data["height_km"], dtype=float).copy()
            earth_radius = float(np.asarray(data["earth_radius_km"]).item())
            periodic_complete = np.asarray(
                data["longitude_periodic_complete"]).item()
            if not isinstance(periodic_complete, (bool, np.bool_)):
                raise ValueError("longitude_periodic_complete must be boolean")
            if "mask" in data:
                mask = np.asarray(data["mask"], dtype=bool)
                if mask.shape != Ne_m3.shape:
                    raise ValueError("mask shape must match Ne")
                Ne_m3 = np.where(mask, Ne_m3, fill_value_m3)
        self.path = str(source.resolve())
        self.source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
        self.source_units = units
        super().__init__(
            Ne_m3, latitude, longitude, height, earth_radius,
            fill_value_m3=fill_value_m3,
            longitude_periodic_complete=periodic_complete,
            source_metadata={
                "source_path": self.path,
                "source_sha256": self.source_sha256,
                "source_units": units,
                "unit_conversion": conversion,
            },
        )

    def get_metadata(self) -> dict:
        return {**super().get_metadata(), "name": self.name}


class VacuumGeocentricNeProvider3D(GeocentricNeProvider3D):
    name = "vacuum_geocentric_ne_provider_3d"

    def __init__(self, earth_radius_km: float):
        geocentric_to_ecef(0.0, 0.0, 0.0, earth_radius_km)
        self.earth_radius_km = float(earth_radius_km)

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        geocentric_to_ecef(lat_deg, lon_deg, height_km, self.earth_radius_km)
        return 0.0

    def get_metadata(self) -> dict:
        frozen = {"model": "vacuum", "Ne_m^-3": 0.0,
                  "earth_radius_km": self.earth_radius_km}
        return {**super().get_metadata(), **frozen,
                "provider_hash": _hash_value(frozen)}


class ChapmanGeocentricNeProvider3D(GeocentricNeProvider3D):
    """Horizontally uniform Chapman F2 reference field on spherical height."""

    name = "chapman_geocentric_ne_provider_3d"

    def __init__(self, earth_radius_km: float, *, foF2_MHz: float = 7.0,
                 hmF2_km: float = 350.0, ymF2_km: float = 100.0):
        geocentric_to_ecef(0.0, 0.0, 0.0, earth_radius_km)
        self.earth_radius_km = float(earth_radius_km)
        self.foF2_MHz = float(foF2_MHz)
        self.hmF2_km = float(hmF2_km)
        self.ymF2_km = float(ymF2_km)
        if not all(np.isfinite(v) and v > 0.0 for v in (
                self.foF2_MHz, self.hmF2_km, self.ymF2_km)):
            raise ValueError("Chapman foF2, hmF2, and ymF2 must be positive finite")
        self.peak_Ne_m3 = (self.foF2_MHz * 1.0e6 / 8.98) ** 2

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        geocentric_to_ecef(lat_deg, lon_deg, height_km, self.earth_radius_km)
        u = (float(height_km) - self.hmF2_km) / self.ymF2_km
        exponent = 0.5 * (1.0 - u - np.exp(np.clip(-u, -700.0, 700.0)))
        return float(self.peak_Ne_m3 * np.exp(np.clip(exponent, -745.0, 700.0)))

    def get_metadata(self) -> dict:
        frozen = {
            "model": "Chapman_F2_horizontal_uniform",
            "foF2_MHz": self.foF2_MHz,
            "hmF2_km": self.hmF2_km,
            "ymF2_km": self.ymF2_km,
            "peak_Ne_m3": self.peak_Ne_m3,
            "earth_radius_km": self.earth_radius_km,
            "density_formula": "Ne=Nmax*exp(0.5*(1-u-exp(-u)))",
            "u": "(height_km-hmF2_km)/ymF2_km",
        }
        return {**super().get_metadata(), **frozen,
                "provider_hash": _hash_value(frozen)}


class GeocentricParameterizedIonosphere3DProvider(GeocentricNeProvider3D):
    """Existing component equations evaluated in a sphere-derived link frame."""

    name = "geocentric_parameterized_ionosphere_3d"

    def __init__(self, background_provider: GeocentricNeProvider3D, *,
                 tx_lat_deg: float, tx_lon_deg: float, link_bearing_deg: float,
                 earth_radius_km: float, tid_params: dict | None = None,
                 bubble_params: dict | None = None,
                 spreadf_params: dict | None = None,
                 es_params: dict | None = None, geomag: dict | None = None,
                 seed: int | None = None, traceability: dict | None = None):
        self.background = background_provider
        self.tx_lat_deg = float(tx_lat_deg)
        self.tx_lon_deg = float(tx_lon_deg)
        self.link_bearing_deg = float(link_bearing_deg) % 360.0
        self.earth_radius_km = float(earth_radius_km)
        geocentric_to_ecef(
            self.tx_lat_deg, self.tx_lon_deg, 0.0, self.earth_radius_km)
        background_metadata = self.background.get_metadata()
        if background_metadata.get("coordinate_system") != COORDINATE_SYSTEM:
            raise ValueError("background Provider must use geocentric_spherical")
        if ("earth_radius_km" not in background_metadata
                or float(background_metadata["earth_radius_km"])
                != self.earth_radius_km):
            raise ValueError("background and parameterized Provider radii must match")
        if not np.isfinite(float(link_bearing_deg)):
            raise ValueError("link_bearing_deg must be finite")
        self.geomag = dict(geomag or {})
        self.seed = seed
        self.traceability = dict(traceability or {})
        params = {
            "tid": dict(tid_params or {}),
            "bubble": dict(bubble_params or {}),
            "spreadf": dict(spreadf_params or {}),
            "es": dict(es_params or {}),
        }
        self.components = []
        if params["tid"].get("enable", False):
            self.components.append(TID3DComponent(params["tid"], self.geomag, seed))
        if params["bubble"].get("enable", False):
            self.components.append(Bubble3DComponent(
                params["bubble"], self.geomag, seed))
        if params["spreadf"].get("enable", False):
            self.components.append(SpreadF3DComponent(
                params["spreadf"], self.geomag, seed))
        if params["es"].get("enable", False):
            es_params_merged = params["es"]
            try:
                import config as cfg
                es_params_merged = {**getattr(cfg, "ES_BG", {}),
                                    **es_params_merged}
            except Exception:
                pass
            self.components.append(Es3DComponent(
                es_params_merged, self.geomag, seed))
        self._tx_up = geocentric_to_ecef(
            self.tx_lat_deg, self.tx_lon_deg, 0.0, self.earth_radius_km)
        self._tx_up /= np.linalg.norm(self._tx_up)
        self._tx_enu = enu_basis(
            self.tx_lat_deg, self.tx_lon_deg, self.earth_radius_km)
        self._b_hat = geomag_unit_vector(self.geomag)

    def derive_link_coordinates(self, lat_deg: float, lon_deg: float,
                                height_km: float) -> np.ndarray:
        return _derive_link_coordinates(
            lat_deg, lon_deg, height_km, self.tx_lat_deg, self.tx_lon_deg,
            self.link_bearing_deg, self.earth_radius_km)

    def get_Ne_geocentric(self, lat_deg: float, lon_deg: float,
                          height_km: float, t: float = 0.0) -> float:
        x, y, z = self.derive_link_coordinates(lat_deg, lon_deg, height_km)
        background = max(float(self.background.get_Ne_geocentric(
            lat_deg, lon_deg, height_km, t=t)), 0.0)
        density = background
        for component in self.components:
            density += component.delta_Ne(float(x), float(y), float(z), t, background)
        return max(float(density), 0.0)

    def domain_margin_geocentric(self, lat_deg: float, lon_deg: float,
                                  height_km: float) -> float:
        return self.background.domain_margin_geocentric(
            lat_deg, lon_deg, height_km)

    def get_metadata(self) -> dict:
        frozen = {
            "tx_lat_deg": self.tx_lat_deg,
            "tx_lon_deg": self.tx_lon_deg,
            "link_bearing_deg": self.link_bearing_deg,
            "earth_radius_km": self.earth_radius_km,
            "coordinate_mapping": "sphere_to_link_azimuthal_equidistant",
            "component_coordinate_semantics": (
                "derived_along_link_cross_link_physical_height_km"),
            "seed": self.seed,
            "components": [component.metadata() for component in self.components],
            "background_provider": self.background.get_metadata(),
            "geomag_b_hat_link_frame": [float(value) for value in self._b_hat],
        }
        return {**super().get_metadata(), **self.traceability, **frozen,
                "provider_hash": _hash_value(frozen)}


__all__ = [
    "AXIS_ORDER",
    "COORDINATE_SYSTEM",
    "GeocentricNeProvider3D",
    "GeocentricGridNeProvider3D",
    "GeocentricLinkGridNeProvider3D",
    "ExternalGeocentricNeProvider3D",
    "VacuumGeocentricNeProvider3D",
    "ChapmanGeocentricNeProvider3D",
    "GeocentricParameterizedIonosphere3DProvider",
]
