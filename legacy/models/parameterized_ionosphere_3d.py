"""Stage-1 component factory for parameterized 3-D ionospheres."""

from __future__ import annotations

from models.ionosphere_components_3d import (
    BaseIonosphere3DComponent,
    Bubble3DComponent,
    Es3DComponent,
    SpreadF3DComponent,
    TID3DComponent,
    build_components,
)
from .ne_provider import (
    BaseNeProvider3D,
    ExternalNeProvider3D,
    GridNeProvider3D,
    ParameterizedIonosphere3DProvider,
)
from models.ne_provider_geocentric import (
    ChapmanGeocentricNeProvider3D,
    ExternalGeocentricNeProvider3D,
    GeocentricGridNeProvider3D,
    GeocentricNeProvider3D,
    GeocentricParameterizedIonosphere3DProvider,
    VacuumGeocentricNeProvider3D,
)


class ParameterizedIonosphere3DFactory:
    """Build a 3-D Ne provider from a background provider and components."""

    def __init__(self, background_provider, geomag: dict | None = None,
                 seed: int | None = None):
        self.background_provider = background_provider
        self.geomag = dict(geomag or {})
        self.seed = seed

    def build(self, components: list[dict]):
        tid = None
        bubble = None
        spreadf = None
        es = None
        for spec in components:
            kind = str(spec.get("type", spec.get("component_type", ""))).lower()
            params = dict(spec.get("params", spec))
            params.setdefault("enable", True)
            if kind in ("tid", "tid3d"):
                tid = params
            elif kind in ("bubble", "bubble3d", "field_aligned_bubble"):
                bubble = params
            elif kind in ("spreadf", "spread-f", "spreadf3d"):
                spreadf = params
            elif kind in ("es", "es3d", "sporadic_e"):
                es = params
            else:
                raise ValueError("unknown ionosphere component type: {}".format(kind))
        provider = ParameterizedIonosphere3DProvider(
            self.background_provider,
            tid_params=tid,
            bubble_params=bubble,
            spreadf_params=spreadf,
            es_params=es,
            geomag=self.geomag,
            seed=self.seed,
        )
        return provider, provider.get_metadata()


__all__ = [
    "BaseNeProvider3D",
    "GridNeProvider3D",
    "ExternalNeProvider3D",
    "ParameterizedIonosphere3DProvider",
    "ParameterizedIonosphere3DFactory",
    "BaseIonosphere3DComponent",
    "TID3DComponent",
    "Bubble3DComponent",
    "SpreadF3DComponent",
    "Es3DComponent",
    "build_components",
    "GeocentricNeProvider3D",
    "GeocentricGridNeProvider3D",
    "ExternalGeocentricNeProvider3D",
    "VacuumGeocentricNeProvider3D",
    "ChapmanGeocentricNeProvider3D",
    "GeocentricParameterizedIonosphere3DProvider",
]
