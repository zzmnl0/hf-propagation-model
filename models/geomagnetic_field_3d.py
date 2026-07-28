"""Traceable three-dimensional geomagnetic fields for ECEF AH ray tracing."""

from __future__ import annotations

from datetime import datetime
import hashlib
import json

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from utils import ecef_to_geodetic, enu_basis_wgs84

GYRO_MHZ_PER_T = 27992.489872333


def _field_hash(metadata: dict, *arrays) -> str:
    digest = hashlib.sha256(json.dumps(
        metadata, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("ascii"))
    for value in arrays:
        array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.tobytes())
    return digest.hexdigest()


class ConstantGeomagneticField3D:
    """Constant ECEF magnetic field retained for analytic regression only."""

    def __init__(self, b_ecef_T):
        field = np.asarray(b_ecef_T, dtype=float)
        if field.shape != (3,) or not np.all(np.isfinite(field)):
            raise ValueError("b_ecef_T must be a finite length-3 vector")
        self.b_ecef_T = field
        self.field_hash = _field_hash(
            {"model": "constant_ecef", "units": "T"}, field)

    def get_B_ecef_T(self, ecef_km, t: float = 0.0) -> np.ndarray:
        del t
        point = np.asarray(ecef_km, dtype=float)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("ECEF position must be a finite length-3 vector")
        return self.b_ecef_T.copy()

    def domain_margin_ecef(self, ecef_km) -> float:
        self.get_B_ecef_T(ecef_km)
        return float("inf")

    def get_metadata(self) -> dict:
        return {
            "name": "constant_geomagnetic_field_3d",
            "model": "constant_ecef",
            "coordinate_system": "ecef",
            "units": "T",
            "field_hash": self.field_hash,
            "validation_only": True,
        }


class IGRFGridGeomagneticField3D:
    """Trilinear WGS84 IGRF grid queried in physical ECEF coordinates."""

    def __init__(self, latitude_deg, longitude_deg, height_km, b_ecef_T,
                 *, epoch: datetime, igrf_version: str = "IGRF-14"):
        self.latitude_deg = np.asarray(latitude_deg, dtype=float)
        self.longitude_deg = np.asarray(longitude_deg, dtype=float)
        self.height_km = np.asarray(height_km, dtype=float)
        self.b_ecef_T = np.asarray(b_ecef_T, dtype=float)
        self.epoch = epoch
        self.igrf_version = str(igrf_version)
        for name, axis in (
            ("latitude_deg", self.latitude_deg),
            ("longitude_deg", self.longitude_deg),
            ("height_km", self.height_km),
        ):
            if (axis.ndim != 1 or axis.size < 2
                    or not np.all(np.isfinite(axis))
                    or not np.all(np.diff(axis) > 0.0)):
                raise ValueError(f"{name} must be finite and strictly increasing")
        expected = (self.latitude_deg.size, self.longitude_deg.size,
                    self.height_km.size, 3)
        if self.b_ecef_T.shape != expected:
            raise ValueError(
                f"b_ecef_T shape {self.b_ecef_T.shape} != {expected}")
        if not np.all(np.isfinite(self.b_ecef_T)):
            raise ValueError("b_ecef_T must be finite")
        metadata = {
            "model": self.igrf_version,
            "epoch": self.epoch.isoformat(),
            "coordinate_system": "wgs84_geodetic_grid_ecef_vectors",
            "units": "T",
            "interpolation": "trilinear",
        }
        self.field_hash = _field_hash(
            metadata, self.latitude_deg, self.longitude_deg, self.height_km,
            self.b_ecef_T)
        self._interpolator = RegularGridInterpolator(
            (self.latitude_deg, self.longitude_deg, self.height_km),
            self.b_ecef_T, method="linear",
            bounds_error=True,
        )

    @classmethod
    def from_ppigrf(cls, latitude_deg, longitude_deg, height_km,
                    *, epoch: datetime):
        import ppigrf

        lat = np.asarray(latitude_deg, dtype=float)
        lon = np.asarray(longitude_deg, dtype=float)
        height = np.asarray(height_km, dtype=float)
        lat_grid, lon_grid, height_grid = np.meshgrid(
            lat, lon, height, indexing="ij")
        east_nT, north_nT, up_nT = ppigrf.igrf(
            lon_grid.ravel(), lat_grid.ravel(), height_grid.ravel(), epoch)
        east_nT = np.asarray(east_nT).reshape(-1)
        north_nT = np.asarray(north_nT).reshape(-1)
        up_nT = np.asarray(up_nT).reshape(-1)
        vectors = np.empty((lat_grid.size, 3), dtype=float)
        for index, (la, lo, east, north, up) in enumerate(zip(
                lat_grid.ravel(), lon_grid.ravel(), east_nT, north_nT, up_nT)):
            vectors[index] = (
                enu_basis_wgs84(float(la), float(lo))
                @ np.array([east, north, up], dtype=float)
                * 1.0e-9
            )
        return cls(lat, lon, height,
                   vectors.reshape(lat.size, lon.size, height.size, 3),
                   epoch=epoch)

    def _geodetic(self, ecef_km) -> np.ndarray:
        point = np.asarray(ecef_km, dtype=float)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("ECEF position must be a finite length-3 vector")
        return ecef_to_geodetic(point)

    def get_B_ecef_T(self, ecef_km, t: float = 0.0) -> np.ndarray:
        del t
        geodetic = self._geodetic(ecef_km)
        point = [[float(geodetic[0]), float(geodetic[1]),
                  float(geodetic[2])]]
        try:
            field = np.asarray(self._interpolator(point)[0], dtype=float)
        except ValueError as exc:
            raise ValueError("IGRF query is outside the frozen field grid") from exc
        if not np.all(np.isfinite(field)) or np.linalg.norm(field) <= 0.0:
            raise ValueError("IGRF interpolation returned an invalid field")
        return field

    def domain_margin_ecef(self, ecef_km) -> float:
        lat, lon, height = self._geodetic(ecef_km)
        return float(min(
            lat - self.latitude_deg[0], self.latitude_deg[-1] - lat,
            lon - self.longitude_deg[0], self.longitude_deg[-1] - lon,
            height - self.height_km[0], self.height_km[-1] - height,
        ))

    def gyrofrequency_MHz(self, ecef_km, t: float = 0.0) -> float:
        return float(GYRO_MHZ_PER_T
                     * np.linalg.norm(self.get_B_ecef_T(ecef_km, t=t)))

    def get_metadata(self) -> dict:
        return {
            "name": "igrf_grid_geomagnetic_field_3d",
            "model": self.igrf_version,
            "epoch": self.epoch.isoformat(),
            "coordinate_system": "wgs84_geodetic_grid_ecef_vectors",
            "units": "T",
            "grid_shape": list(self.b_ecef_T.shape),
            "latitude_range_deg": [float(self.latitude_deg[0]),
                                   float(self.latitude_deg[-1])],
            "longitude_range_deg": [float(self.longitude_deg[0]),
                                    float(self.longitude_deg[-1])],
            "height_range_km": [float(self.height_km[0]),
                                float(self.height_km[-1])],
            "interpolation": "trilinear",
            "field_hash": self.field_hash,
            "validation_only": False,
        }


__all__ = [
    "GYRO_MHZ_PER_T",
    "ConstantGeomagneticField3D",
    "IGRFGridGeomagneticField3D",
]
