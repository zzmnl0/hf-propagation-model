"""Electron-density provider interfaces for update2.1."""

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
        """Return Ne at ``(x,z)`` or ``(x,y,z)`` in km."""

    def get_Ne_grid(
        self,
        x_km: Iterable[float],
        z_km: Iterable[float],
        t: float = 0.0,
        y_km: float = 0.0,
    ) -> np.ndarray:
        """Return an ``(Nx, Nz)`` Ne grid in m^-3."""
        x_arr = np.asarray(x_km, dtype=float)
        z_arr = np.asarray(z_km, dtype=float)
        out = np.empty((x_arr.size, z_arr.size), dtype=float)
        for i, x_val in enumerate(x_arr):
            for j, z_val in enumerate(z_arr):
                out[i, j] = self.get_Ne((x_val, y_km, z_val), t=t)
        return out


class GridNeProvider(BaseNeProvider):
    """Provider backed by an existing ``Ne(x,z)`` grid."""

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


ExistingGridProvider = GridNeProvider


class IriNeProvider(GridNeProvider):
    """IRI-backed provider sampled once onto an ``(x,z)`` grid."""

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
