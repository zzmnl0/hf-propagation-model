"""Electron-density provider interfaces for stage-1 3-D hybrid propagation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from typing import Iterable

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class BaseNeProvider(ABC):
    """Abstract electron-density provider returning Ne in m^-3."""

    name = "base"

    @abstractmethod
    def get_Ne(self, pos_km, t: float = 0.0) -> float:
        """Return Ne at (x,z) or (x,y,z) in km."""

    def get_Ne_grid(
        self,
        x_km: Iterable[float],
        z_km: Iterable[float],
        t: float = 0.0,
        y_km: float = 0.0,
    ) -> np.ndarray:
        """Return an (Nx, Nz) Ne grid in m^-3 at fixed y."""
        x_arr = np.asarray(x_km, dtype=float)
        z_arr = np.asarray(z_km, dtype=float)
        out = np.empty((x_arr.size, z_arr.size), dtype=float)
        for i, x_val in enumerate(x_arr):
            for j, z_val in enumerate(z_arr):
                out[i, j] = self.get_Ne((x_val, y_km, z_val), t=t)
        return out

    def get_grad_Ne(self, pos_km, t: float = 0.0, step_km: float = 1.0) -> np.ndarray:
        """Return central-difference grad(Ne) in m^-3 km^-1 at (x,y,z)."""
        arr = np.asarray(pos_km, dtype=float).ravel()
        if arr.size == 2:
            arr = np.array([arr[0], 0.0, arr[1]], dtype=float)
        elif arr.size >= 3:
            arr = arr[:3].astype(float)
        else:
            raise ValueError("pos_km must contain (x,z) or (x,y,z)")
        h = float(step_km)
        grad = np.empty(3, dtype=float)
        for i in range(3):
            pp = arr.copy()
            pm = arr.copy()
            pp[i] += h
            pm[i] -= h
            grad[i] = (self.get_Ne(pp, t=t) - self.get_Ne(pm, t=t)) / (2.0 * h)
        return grad

    def get_metadata(self) -> dict:
        """Return provider metadata for traceability."""
        return {"name": self.name, "dimensions": "2d-compatible", "units": "m^-3"}

    def sample_grid3d(self, x_km, y_km, z_km, t: float = 0.0) -> np.ndarray:
        """Return an (Nx, Ny, Nz) Ne grid in m^-3."""
        x_arr = np.asarray(x_km, dtype=float)
        y_arr = np.asarray(y_km, dtype=float)
        z_arr = np.asarray(z_km, dtype=float)
        out = np.empty((x_arr.size, y_arr.size, z_arr.size), dtype=float)
        for i, x_val in enumerate(x_arr):
            for j, y_val in enumerate(y_arr):
                for k, z_val in enumerate(z_arr):
                    out[i, j, k] = self.get_Ne((x_val, y_val, z_val), t=t)
        return out

    def sample_ray_path(self, points_km, t: float = 0.0) -> np.ndarray:
        """Sample Ne along ray points with columns x,y,z."""
        pts = np.asarray(points_km, dtype=float)
        return np.asarray([self.get_Ne(p[:3], t=t) for p in pts], dtype=float)


class BaseNeProvider3D(BaseNeProvider):
    """Abstract 3-D electron-density provider returning Ne in m^-3."""

    name = "base3d"


class GridNeProvider(BaseNeProvider):
    """Provider backed by an existing Ne(x,z) grid."""

    name = "grid"

    def __init__(
        self,
        Ne_2d: np.ndarray,
        x_km: Iterable[float],
        z_km: Iterable[float],
        fill_value: float = 0.0,
    ):
        self.Ne_2d = np.asarray(Ne_2d, dtype=float)
        self.x_km = np.asarray(x_km, dtype=float)
        self.z_km = np.asarray(z_km, dtype=float)
        self._interp = RegularGridInterpolator(
            (self.x_km, self.z_km),
            self.Ne_2d,
            method="linear",
            bounds_error=False,
            fill_value=fill_value,
        )

    def get_Ne(self, pos_km, t: float = 0.0) -> float:
        arr = np.asarray(pos_km, dtype=float).ravel()
        if arr.size == 2:
            x_val, z_val = arr
        elif arr.size >= 3:
            x_val, z_val = arr[0], arr[2]
        else:
            raise ValueError("pos_km must contain (x,z) or (x,y,z)")
        return max(float(self._interp([[x_val, z_val]])[0]), 0.0)

    def get_Ne_grid(
        self,
        x_km: Iterable[float],
        z_km: Iterable[float],
        t: float = 0.0,
        y_km: float = 0.0,
    ) -> np.ndarray:
        x_arr = np.asarray(x_km, dtype=float)
        z_arr = np.asarray(z_km, dtype=float)
        xx, zz = np.meshgrid(x_arr, z_arr, indexing="ij")
        pts = np.column_stack([xx.ravel(), zz.ravel()])
        return np.maximum(self._interp(pts).reshape(x_arr.size, z_arr.size), 0.0)

    def get_metadata(self) -> dict:
        return {
            "name": self.name,
            "dimensions": "2d",
            "units": "m^-3",
            "x_range_km": [float(self.x_km[0]), float(self.x_km[-1])],
            "z_range_km": [float(self.z_km[0]), float(self.z_km[-1])],
        }


ExistingGridProvider = GridNeProvider


class GridNeProvider3D(BaseNeProvider3D):
    """Provider backed by a regular Ne(x,y,z) grid."""

    name = "grid3d"

    def __init__(self, Ne_3d, x_km, y_km, z_km, fill_value: float = 0.0):
        self.Ne_3d = np.asarray(Ne_3d, dtype=float)
        self.x_km = np.asarray(x_km, dtype=float)
        self.y_km = np.asarray(y_km, dtype=float)
        self.z_km = np.asarray(z_km, dtype=float)
        self._interp = RegularGridInterpolator(
            (self.x_km, self.y_km, self.z_km),
            self.Ne_3d,
            method="linear",
            bounds_error=False,
            fill_value=fill_value,
        )

    def get_Ne(self, pos_km, t: float = 0.0) -> float:
        arr = np.asarray(pos_km, dtype=float).ravel()
        if arr.size == 2:
            arr = np.array([arr[0], 0.0, arr[1]], dtype=float)
        elif arr.size >= 3:
            arr = arr[:3].astype(float)
        else:
            raise ValueError("pos_km must contain (x,z) or (x,y,z)")
        return max(float(self._interp([arr])[0]), 0.0)

    def get_metadata(self) -> dict:
        return {
            "name": self.name,
            "dimensions": "3d",
            "units": "m^-3",
            "x_range_km": [float(self.x_km[0]), float(self.x_km[-1])],
            "y_range_km": [float(self.y_km[0]), float(self.y_km[-1])],
            "z_range_km": [float(self.z_km[0]), float(self.z_km[-1])],
        }


class ParameterizedIonosphere3DProvider(BaseNeProvider3D):
    """3-D parameterized complex ionosphere provider."""

    name = "parameterized_ionosphere_3d"

    def __init__(
        self,
        background_provider: BaseNeProvider,
        tid_params: dict | None = None,
        bubble_params: dict | None = None,
        spreadf_params: dict | None = None,
        es_params: dict | None = None,
        geomag: dict | None = None,
        seed: int | None = None,
    ):
        self.background = background_provider
        self.tid = dict(tid_params or {})
        self.bubble = dict(bubble_params or {})
        self.spreadf = dict(spreadf_params or {})
        self.es = dict(es_params or {})
        self.geomag = dict(geomag or {})
        self.seed = seed
        self._b_hat = self._geomag_unit_vector()

    def _geomag_unit_vector(self) -> np.ndarray:
        dip = np.radians(float(self.geomag.get("dip_deg", 48.7)))
        decl = np.radians(float(self.geomag.get("decl_deg", -5.5)))
        link = np.radians(float(self.geomag.get("link_bearing_deg", 0.0)))
        b = np.array([
            np.cos(dip) * np.cos(decl - link),
            np.cos(dip) * np.sin(decl - link),
            -np.sin(dip),
        ], dtype=float)
        return b / max(float(np.linalg.norm(b)), 1e-30)

    def _tid_delta(self, x: float, y: float, z: float, t: float, Ne_bg: float) -> float:
        if not self.tid.get("enable", False):
            return 0.0
        n_comp = int(self.tid.get("n_components", 1))
        amp_list = self.tid.get("amplitude_list", [float(self.tid.get("amplitude", 0.1))])
        az_list = self.tid.get("az_deg_list", [0.0])
        per_list = self.tid.get("period_s_list", [float(self.tid.get("T_s", 2400.0))])
        lam_list = self.tid.get("lambda_h_km_list", [float(self.tid.get("lambda_h_km", 300.0))])
        lam_z = float(self.tid.get("lambda_z_km", 180.0))
        z0 = float(self.tid.get("z0_km", 250.0))
        sigma_z = float(self.tid.get("sigma_z_km", 90.0))
        link = float(self.tid.get("link_bearing_deg", 0.0))
        envelope = np.exp(-((z - z0) / max(sigma_z, 1e-6)) ** 2)
        delta = 0.0
        for i in range(n_comp):
            amp = float(amp_list[i % len(amp_list)])
            az = float(az_list[i % len(az_list)]) - link
            period = float(per_list[i % len(per_list)])
            lam_h = float(lam_list[i % len(lam_list)])
            kh = 2.0 * np.pi / max(lam_h, 1e-6)
            kz = 2.0 * np.pi / max(lam_z, 1e-6)
            phase = kh * (x * np.cos(np.radians(az)) + y * np.sin(np.radians(az))) + kz * z - 2.0 * np.pi * t / period
            delta += amp * Ne_bg * envelope * np.cos(phase)
        return float(delta)

    def _bubble_factor(self, x: float, y: float, z: float) -> float:
        if not self.bubble.get("enable", False):
            return 0.0
        r = np.array([
            x - float(self.bubble.get("x0_km", 0.0)),
            y - float(self.bubble.get("y0_km", 0.0)),
            z - float(self.bubble.get("z0_km", 300.0)),
        ], dtype=float)
        r_para = float(np.dot(r, self._b_hat))
        r_perp = r - r_para * self._b_hat
        sig_perp = max(float(self.bubble.get("sigma_perp_km", self.bubble.get("Lx_km", 50.0))), 1e-6)
        sig_para = max(float(self.bubble.get("sigma_para_km", 500.0)), 1e-6)
        frac = float(self.bubble.get("delta_max", 0.0)) * np.exp(
            -float(np.dot(r_perp, r_perp)) / sig_perp ** 2 - r_para ** 2 / sig_para ** 2
        )
        return float(np.clip(frac, 0.0, 1.0))

    def get_Ne(self, pos_km, t: float = 0.0) -> float:
        arr = np.asarray(pos_km, dtype=float).ravel()
        if arr.size == 2:
            x, y, z = float(arr[0]), 0.0, float(arr[1])
        elif arr.size >= 3:
            x, y, z = map(float, arr[:3])
        else:
            raise ValueError("pos_km must contain (x,z) or (x,y,z)")
        Ne_bg = max(float(self.background.get_Ne((x, y, z), t=t)), 0.0)
        Ne = max(Ne_bg + self._tid_delta(x, y, z, t, Ne_bg), 0.0)
        Ne *= 1.0 - self._bubble_factor(x, y, z)
        return max(float(Ne), 0.0)

    def get_metadata(self) -> dict:
        return {
            "name": self.name,
            "dimensions": "3d",
            "units": "m^-3",
            "background": self.background.get_metadata(),
            "tid_enabled": bool(self.tid.get("enable", False)),
            "bubble_enabled": bool(self.bubble.get("enable", False)),
            "spreadf_diagnostic_only": bool(self.spreadf.get("enable", False)),
            "es_diagnostic_only": bool(self.es.get("enable", False)),
            "seed": self.seed,
            "geomag_b_hat_xyz": [float(v) for v in self._b_hat],
        }

    def perturbation_summary(self, x_km, y_km, z_km, t: float = 0.0) -> dict:
        ratios = []
        for x in np.asarray(x_km, dtype=float):
            for y in np.asarray(y_km, dtype=float):
                for z in np.asarray(z_km, dtype=float):
                    bg = max(float(self.background.get_Ne((x, y, z), t=t)), 0.0)
                    if bg <= 0.0:
                        continue
                    ratios.append((self.get_Ne((x, y, z), t=t) - bg) / bg)
        arr = np.asarray(ratios, dtype=float) if ratios else np.zeros(1)
        return {
            "max_relative_perturbation": float(np.max(np.abs(arr))),
            "min_relative_perturbation": float(np.min(arr)),
            "max_positive_relative_perturbation": float(np.max(arr)),
            "metadata": self.get_metadata(),
        }


class IriNeProvider(GridNeProvider):
    """IRI-backed provider sampled once onto an (x,z) grid."""

    name = "iri2016_grid"

    def __init__(
        self,
        x_km: Iterable[float],
        z_km: Iterable[float],
        lat: float,
        lon: float,
        dt: datetime,
    ):
        import iri2016

        x_arr = np.asarray(x_km, dtype=float)
        z_arr = np.asarray(z_km, dtype=float)
        dz = float(z_arr[1] - z_arr[0]) if z_arr.size > 1 else 1.0
        iri = iri2016.IRI(dt, (float(z_arr[0]), float(z_arr[-1]), dz), lat, lon)
        profile = np.asarray(iri["ne"].values, dtype=float)
        if profile.size != z_arr.size:
            src_z = np.linspace(float(z_arr[0]), float(z_arr[-1]), profile.size)
            profile = np.interp(z_arr, src_z, profile)
        Ne_2d = np.repeat(profile[None, :], x_arr.size, axis=0)
        super().__init__(Ne_2d, x_arr, z_arr)
