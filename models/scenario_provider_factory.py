"""Deterministic construction of parameterized 3-D ionosphere scenarios."""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import asdict, dataclass
from datetime import date, datetime
from typing import Iterable

import numpy as np

from .ionosphere_model import IonosphereModel
from .ne_provider_geocentric import (
    ExternalGeocentricNeProvider3D,
    GeocentricLinkGridNeProvider3D,
    GeocentricParameterizedIonosphere3DProvider,
)


SCENARIO_COMPONENTS = {
    "baseline": (),
    "tid": ("tid",),
    "spreadf": ("spreadf",),
    "bubble": ("bubble",),
    "es": ("es",),
}

SCENARIO_ALIASES = {
    "iri": "baseline",
    "spread-f": "spreadf",
    "sporadic-e": "es",
}


def normalize_scenario_label(label: str) -> str:
    key = str(label).strip().lower().replace("_", "-")
    key = SCENARIO_ALIASES.get(key, key)
    if key not in SCENARIO_COMPONENTS:
        raise ValueError("unknown scenario {!r}; choose from {}".format(
            label, ", ".join(SCENARIO_COMPONENTS)))
    return key


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, dict):
        return {str(k): _jsonable(value[k]) for k in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _hash_json(value) -> str:
    payload = json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _hash_provider(config_hash: str, x_km, z_km, ne_background) -> str:
    digest = hashlib.sha256(config_hash.encode("ascii"))
    for value in (x_km, z_km, ne_background):
        arr = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest.update(str(arr.shape).encode("ascii"))
        digest.update(arr.tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class ScenarioSpec:
    scenario: str
    components: tuple[str, ...]
    seed: int
    t_s: float
    coordinate_system: str = "geocentric_spherical"
    ne_units: str = "m^-3"
    time_units: str = "s"
    interpolation: str = "linear_horizontal_pchip_vertical"
    bounds_error: bool = False
    extrapolation: str = "fill_zero"


class ScenarioProviderFactory:
    """Create one reproducible 3-D provider from a frozen config snapshot."""

    def __init__(self, config_module):
        self.cfg = config_module
        self._background_cache = None

    def _background_field(self, x_km, z_km):
        signature = _hash_json({
            "x_km": x_km, "z_km": z_km,
            "iri_dt": getattr(self.cfg, "IRI_DT", None),
            "iri_lat": getattr(self.cfg, "IRI_LAT", None),
            "iri_lon": getattr(self.cfg, "IRI_LON", None),
        })
        if self._background_cache is None or self._background_cache[0] != signature:
            values, _ = IonosphereModel().build_Ne_field(x_km, z_km)
            frozen = np.asarray(values, dtype=float).copy()
            frozen.setflags(write=False)
            self._background_cache = (signature, frozen)
        return self._background_cache[1]

    def create(self, scenario: str = "baseline", *, components: Iterable[str] | None = None,
               seed: int | None = 0, t_s: float = 0.0,
               validation: bool = True) -> GeocentricParameterizedIonosphere3DProvider:
        label = normalize_scenario_label(scenario)
        requested = SCENARIO_COMPONENTS[label] if components is None else tuple(
            str(v).strip().lower().replace("-", "") for v in components)
        allowed = {"tid", "spreadf", "bubble", "es"}
        unknown = set(requested) - allowed
        if unknown:
            raise ValueError("unknown ionosphere components: {}".format(sorted(unknown)))
        selected = tuple(v for v in ("tid", "spreadf", "bubble", "es")
                         if v in requested)
        if validation and seed is None:
            raise ValueError("validation providers require an explicit integer seed")
        seed_value = 0 if seed is None else int(seed)

        external_cfg = self.cfg.UPDATE2.get("EXTERNAL_NE_PROVIDER", {})
        if external_cfg.get("enable", False):
            if selected:
                raise ValueError(
                    "external Ne already contains unknown irregularities; scenario components must be empty")
            path = external_cfg.get("path")
            if not path:
                raise ValueError("external geocentric Ne requires an explicit NPZ path")
            provider = ExternalGeocentricNeProvider3D(path)
            metadata = provider.get_metadata()
            if float(metadata["earth_radius_km"]) != float(self.cfg.RE_KM):
                raise ValueError("external Ne Earth radius differs from config EARTH_MODEL")
            return provider

        x_km = np.asarray(self.cfg.BG_X, dtype=float).copy()
        z_km = np.asarray(self.cfg.BG_Z, dtype=float).copy()
        # The current IRI background is range-height sampled and y-invariant. All
        # perturbations and all propagation queries are evaluated by the 3-D provider.
        ne_background = self._background_field(x_km, z_km)
        background = GeocentricLinkGridNeProvider3D(
            ne_background, x_km, z_km,
            tx_lat_deg=self.cfg.TX_GEOCENTRIC[0],
            tx_lon_deg=self.cfg.TX_GEOCENTRIC[1],
            link_bearing_deg=self.cfg.LINK_BEARING_DEG,
            earth_radius_km=self.cfg.EARTH_MODEL["radius_km"],
            fill_value_m3=0.0,
            vacuum_floor_height_km=0.0,
            vertical_interpolation="pchip",
            source_metadata={
                "model": "IRI-2016",
                "datetime": self.cfg.IRI_DT.isoformat(),
                "formal_input_coordinates": "geocentric_spherical",
                "derived_grid_axes": "great_circle_ground_range,physical_height",
            })

        spec = ScenarioSpec(label, tuple(selected), seed_value, float(t_s))
        params = {
            "tid": copy.deepcopy(self.cfg.TID_3D),
            "spreadf": copy.deepcopy(self.cfg.SPREAD_F),
            "bubble": copy.deepcopy(self.cfg.BUBBLE_3D),
            "es": copy.deepcopy(self.cfg.ES),
            "geomag": copy.deepcopy(self.cfg.GEOMAG),
        }
        active = {
            name: {**params[name], "enable": name in selected}
            for name in ("tid", "spreadf", "bubble", "es")
        }
        config_snapshot = {"spec": asdict(spec), "active_parameters": active,
                           "geomag": params["geomag"]}
        config_hash = _hash_json(config_snapshot)
        physics_hash = _hash_json({
            "components": list(selected), "seed": seed_value, "t_s": float(t_s),
            "active_parameters": active, "geomag": params["geomag"],
        })
        provider_hash = _hash_provider(physics_hash, x_km, z_km, ne_background)
        spec_metadata = asdict(spec)
        spec_metadata["scenario_components"] = spec_metadata.pop("components")
        traceability = {
            **spec_metadata,
            "provider_name": GeocentricParameterizedIonosphere3DProvider.name,
            "config_hash": config_hash,
            "provider_hash": provider_hash,
            "background_semantics": "range_height_sampled_y_invariant",
            "background_grid_shape": list(ne_background.shape),
        }
        return GeocentricParameterizedIonosphere3DProvider(
            background,
            tx_lat_deg=self.cfg.TX_GEOCENTRIC[0],
            tx_lon_deg=self.cfg.TX_GEOCENTRIC[1],
            link_bearing_deg=self.cfg.LINK_BEARING_DEG,
            earth_radius_km=self.cfg.EARTH_MODEL["radius_km"],
            tid_params=active["tid"],
            spreadf_params=active["spreadf"],
            bubble_params=active["bubble"],
            es_params=active["es"],
            geomag=params["geomag"],
            seed=seed_value,
            traceability=traceability,
        )
