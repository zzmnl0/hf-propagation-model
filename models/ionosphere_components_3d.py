"""Parameterized 3-D ionospheric electron-density components.

The components in this module modify Ne(x,y,z,t) only.  Propagation,
reflection, absorption, and scattering remain in the unified propagation
pipeline and downstream solvers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


LITERATURE_REFERENCES = {
    "tid": (
        "Hines, C. O. (1960). Internal atmospheric gravity waves at ionospheric "
        "heights. Canadian Journal of Physics, 38(11), 1441-1481."
    ),
    "spreadf": (
        "Rino, C. L. (1979). A power law phase screen model for ionospheric "
        "scintillation. Radio Science, 14(6), 1135-1145."
    ),
    "scintillation": (
        "Yeh, K. C., & Liu, C. H. (1982). Radio wave scintillations in the "
        "ionosphere. Proceedings of the IEEE, 70(4), 324-360."
    ),
    "es": (
        "Haldoupis, C. (2011). A tutorial review on sporadic E layers. In "
        "Aeronomy of the Earth's Atmosphere and Ionosphere, 381-394. Springer."
    ),
    "booker_gordon": (
        "Booker, H. G., & Gordon, W. E. (1950). A theory of radio scattering "
        "in the troposphere. Proceedings of the IRE, 38(4), 401-412."
    ),
    "bubble": (
        "Kelley, M. C. (2009). The Earth's Ionosphere: Plasma Physics and "
        "Electrodynamics (2nd ed.). Academic Press."
    ),
}


def geomag_unit_vector(geomag: dict | None = None) -> np.ndarray:
    g = dict(geomag or {})
    dip = np.radians(float(g.get("dip_deg", 48.7)))
    decl = np.radians(float(g.get("decl_deg", -5.5)))
    link = np.radians(float(g.get("link_bearing_deg", 0.0)))
    b = np.array([
        np.cos(dip) * np.cos(decl - link),
        np.cos(dip) * np.sin(decl - link),
        -np.sin(dip),
    ], dtype=float)
    return b / max(float(np.linalg.norm(b)), 1e-30)


def ne_from_plasma_frequency_mhz(freq_MHz: float) -> float:
    """Convert plasma frequency in MHz to Ne in m^-3."""
    return float((float(freq_MHz) * 1.0e6 / 8.98) ** 2)


@dataclass
class BaseIonosphere3DComponent:
    params: dict
    geomag: dict | None = None
    seed: int | None = None

    component_type = "base"
    physics_basis = "component"
    validity = "parameterized_reference"

    def __post_init__(self) -> None:
        pass

    def enabled(self) -> bool:
        return bool(self.params.get("enable", False))

    def delta_Ne(self, x: float, y: float, z: float, t: float, Ne_bg: float) -> float:
        return 0.0

    def metadata(self) -> dict[str, Any]:
        return {
            "component_type": self.component_type,
            "enabled": self.enabled(),
            "units": "m^-3",
            "seed": self.seed,
            "coordinate_system": "local_cartesian_km",
            "validity": self.validity,
            "physics_basis": self.physics_basis,
        }


class TID3DComponent(BaseIonosphere3DComponent):
    component_type = "TID3D"
    physics_basis = LITERATURE_REFERENCES["tid"]

    def delta_Ne(self, x: float, y: float, z: float, t: float, Ne_bg: float) -> float:
        if not self.enabled():
            return 0.0
        n_comp = int(self.params.get("n_components", 1))
        amps = self.params.get("amplitude_list", [self.params.get("amplitude", 0.1)])
        azs = self.params.get("az_deg_list", [0.0])
        periods = self.params.get("period_s_list", [self.params.get("T_s", 2400.0)])
        lambdas = self.params.get("lambda_h_km_list", [self.params.get("lambda_h_km", 300.0)])
        phases = self.params.get("phase_rad_list", [self.params.get("phase_rad", 0.0)])
        # H_km is the legacy TID scale-height parameter, not the F-region
        # envelope centre.  Keep the two physical quantities distinct.
        z0 = float(self.params.get("z0_km", 280.0))
        sigma_z = float(self.params.get("sigma_z_km", 70.0))
        lambda_z = float(self.params.get("lambda_z_km", 1.0e30))
        link = float(self.params.get("link_bearing_deg", 0.0))
        envelope = np.exp(-0.5 * ((z - z0) / max(sigma_z, 1e-6)) ** 2)
        total = 0.0
        for i in range(n_comp):
            amp = float(amps[i % len(amps)])
            az = np.radians(float(azs[i % len(azs)]) - link)
            period = max(float(periods[i % len(periods)]), 1e-9)
            lam_h = max(float(lambdas[i % len(lambdas)]), 1e-9)
            kh = 2.0 * np.pi / lam_h
            kz = 0.0 if lambda_z > 1.0e20 else 2.0 * np.pi / max(lambda_z, 1e-9)
            phase0 = float(phases[i % len(phases)])
            phase = (
                kh * (x * np.cos(az) + y * np.sin(az))
                + kz * z
                - 2.0 * np.pi * float(t) / period
                + phase0
            )
            total += amp * Ne_bg * envelope * np.cos(phase)
        return float(total)

    def metadata(self) -> dict[str, Any]:
        meta = super().metadata()
        lambdas = self.params.get("lambda_h_km_list", [self.params.get("lambda_h_km", 300.0)])
        periods = self.params.get("period_s_list", [self.params.get("T_s", 2400.0)])
        meta.update({
            "n_components": int(self.params.get("n_components", 1)),
            "max_relative_perturbation": float(max(abs(v) for v in self.params.get(
                "amplitude_list", [self.params.get("amplitude", 0.0)]))),
            "lambda_h_km_list": [float(v) for v in lambdas],
            "period_s_list": [float(v) for v in periods],
            "phase_speed_ms_list": [
                float(lambdas[i % len(lambdas)]) * 1000.0 / max(float(periods[i % len(periods)]), 1e-9)
                for i in range(int(self.params.get("n_components", 1)))
            ],
            "effective_height_km": float(self.params.get("z0_km", 280.0)),
            "sigma_z_km": float(self.params.get("sigma_z_km", 70.0)),
        })
        return meta


class Bubble3DComponent(BaseIonosphere3DComponent):
    component_type = "field_aligned_bubble"
    physics_basis = LITERATURE_REFERENCES["bubble"]

    def __post_init__(self) -> None:
        self._b_hat = geomag_unit_vector(self.geomag)

    def depletion_fraction(self, x: float, y: float, z: float, t: float) -> float:
        if not self.enabled():
            return 0.0
        vx = float(self.params.get("drift_x_kms", 0.0))
        vy = float(self.params.get("drift_y_kms", 0.0))
        vz = float(self.params.get("drift_z_kms", 0.0))
        center = np.array([
            float(self.params.get("x0_km", 0.0)) + vx * float(t),
            float(self.params.get("y0_km", 0.0)) + vy * float(t),
            float(self.params.get("z0_km", 300.0)) + vz * float(t),
        ], dtype=float)
        r = np.array([x, y, z], dtype=float) - center
        r_para = float(np.dot(r, self._b_hat))
        r_perp = r - r_para * self._b_hat
        sig_perp = max(float(self.params.get("sigma_perp_km", self.params.get("Lx_km", 50.0))), 1e-6)
        sig_para = max(float(self.params.get("sigma_para_km", self.params.get("Lz_km", 500.0))), 1e-6)
        frac = float(self.params.get("delta_max", 0.0)) * np.exp(
            -0.5 * float(np.dot(r_perp, r_perp)) / sig_perp ** 2
            -0.5 * r_para ** 2 / sig_para ** 2
        )
        return float(np.clip(frac, 0.0, 1.0))

    def delta_Ne(self, x: float, y: float, z: float, t: float, Ne_bg: float) -> float:
        return -Ne_bg * self.depletion_fraction(x, y, z, t)

    def metadata(self) -> dict[str, Any]:
        meta = super().metadata()
        meta.update({
            "delta_max": float(self.params.get("delta_max", 0.0)),
            "sigma_perp_km": float(self.params.get("sigma_perp_km", self.params.get("Lx_km", 50.0))),
            "sigma_para_km": float(self.params.get("sigma_para_km", self.params.get("Lz_km", 500.0))),
            "anisotropy_ratio": float(self.params.get("sigma_para_km", self.params.get("Lz_km", 500.0)))
            / max(float(self.params.get("sigma_perp_km", self.params.get("Lx_km", 50.0))), 1e-12),
            "b_hat": [float(v) for v in self._b_hat],
        })
        return meta


class _SpectralModeField:
    def __init__(self, params: dict, geomag: dict | None, seed: int | None):
        self.params = params
        self.geomag = geomag
        self.seed = seed
        self.b_hat = geomag_unit_vector(geomag)
        self.n_modes = int(params.get("n_modes", 96))
        rng = np.random.default_rng(seed)
        L0 = max(float(params.get("L0_km", params.get("outer_scale_km", 50.0))), 1e-6)
        l0_km = max(float(params.get("l0_m", params.get("inner_scale_m", 100.0))) / 1000.0, 1e-6)
        k_min = 2.0 * np.pi / L0
        k_max = 2.0 * np.pi / l0_km
        k_mag = np.exp(rng.uniform(np.log(k_min), np.log(max(k_max, k_min * 1.01)), self.n_modes))
        dirs = rng.normal(size=(self.n_modes, 3))
        dirs /= np.maximum(np.linalg.norm(dirs, axis=1)[:, None], 1e-30)
        anis = max(float(params.get("anisotropy_ratio", params.get("field_aligned_anisotropy", 20.0))), 1.0)
        k_para = dirs @ self.b_hat
        dirs = dirs - (1.0 - 1.0 / anis) * k_para[:, None] * self.b_hat[None, :]
        dirs /= np.maximum(np.linalg.norm(dirs, axis=1)[:, None], 1e-30)
        self.k_vec = dirs * k_mag[:, None]
        p = float(params.get("p", params.get("spectral_index", 3.0)))
        k0 = 2.0 * np.pi / L0
        weights = (k_mag ** 2 + k0 ** 2) ** (-0.25 * (p + 3.0))
        self.weights = weights / max(float(np.sqrt(np.mean(weights ** 2))), 1e-30)
        self.phase = rng.uniform(0.0, 2.0 * np.pi, self.n_modes)

    def value(self, x: float, y: float, z: float, t: float) -> float:
        drift = np.array([
            float(self.params.get("drift_x_kms", 0.0)),
            float(self.params.get("drift_y_kms", 0.0)),
            float(self.params.get("drift_z_kms", 0.0)),
        ], dtype=float)
        r = np.array([x, y, z], dtype=float) - drift * float(t)
        arg = self.k_vec @ r + self.phase
        val = np.sum(self.weights * np.cos(arg))
        return float(val / max(np.sqrt(0.5 * np.sum(self.weights ** 2)), 1e-30))


class SpreadF3DComponent(BaseIonosphere3DComponent):
    component_type = "SpreadF3D"
    physics_basis = LITERATURE_REFERENCES["spreadf"] + " " + LITERATURE_REFERENCES["scintillation"]

    def __post_init__(self) -> None:
        self._field = _SpectralModeField(self.params, self.geomag, self.seed)

    def delta_Ne(self, x: float, y: float, z: float, t: float, Ne_bg: float) -> float:
        if not self.enabled():
            return 0.0
        h = float(self.params.get("h_screen_km", self.params.get("z0_km", 300.0)))
        width = max(float(self.params.get("sigma_z_km", self.params.get("thickness_km", 40.0))), 1e-6)
        envelope = np.exp(-0.5 * ((z - h) / width) ** 2)
        # Cs is a spectral-strength parameter and is not a dimensionless
        # density fluctuation.  Density RMS must be supplied independently.
        rel = float(self.params.get("delta_Ne_rel_rms", 0.05))
        val = rel * self._field.value(x, y, z, t) * envelope
        return float(Ne_bg * np.clip(val, -0.95, 2.0))

    def metadata(self) -> dict[str, Any]:
        meta = super().metadata()
        meta.update({
            "Cs": float(self.params.get("Cs", self.params.get("Cs_rel", 0.01))),
            "delta_Ne_rel_rms": float(self.params.get("delta_Ne_rel_rms", 0.05)),
            "spectral_index": float(self.params.get("p", self.params.get("spectral_index", 3.0))),
            "outer_scale_km": float(self.params.get("L0_km", self.params.get("outer_scale_km", 50.0))),
            "inner_scale_m": float(self.params.get("l0_m", self.params.get("inner_scale_m", 100.0))),
            "sigma_phi": float(self.params.get("sigma_phi", 0.0)),
            "anisotropy_ratio": float(self.params.get("anisotropy_ratio", self.params.get("field_aligned_anisotropy", 20.0))),
            "n_modes": int(self.params.get("n_modes", 96)),
            "field_model": "3d_anisotropic_von_karman_random_modes",
        })
        return meta


class Es3DComponent(BaseIonosphere3DComponent):
    component_type = "Es3D"
    physics_basis = LITERATURE_REFERENCES["es"] + " " + LITERATURE_REFERENCES["booker_gordon"]

    def __post_init__(self) -> None:
        self._b_hat = geomag_unit_vector(self.geomag)
        self._field = _SpectralModeField({
            "L0_km": float(self.params.get("L_perp_m", self.params.get("L1_m", 300.0))) / 1000.0,
            "l0_m": max(float(self.params.get("L3_m", 30.0)), 1.0),
            "p": float(self.params.get("p_es", 1.5)),
            "anisotropy_ratio": max(
                float(self.params.get("L_para_m", 10000.0))
                / max(float(self.params.get("L_perp_m", self.params.get("L1_m", 300.0))), 1.0),
                1.0,
            ),
            "n_modes": int(self.params.get("n_modes", 48)),
        }, self.geomag, self.seed)

    def delta_Ne(self, x: float, y: float, z: float, t: float, Ne_bg: float) -> float:
        if not self.enabled():
            return 0.0
        h = float(self.params.get("h_Es_km", 110.0))
        sigma_z = max(float(self.params.get("delta_h_m", 115.0)) / 1000.0, 1e-4)
        peak = ne_from_plasma_frequency_mhz(float(self.params.get("foEs_MHz", 5.0)))
        layer = peak * np.exp(-0.5 * ((z - h) / sigma_z) ** 2)
        x0 = float(self.params.get("x0_km", self.params.get("patch_x0_km", x)))
        y0 = float(self.params.get("y0_km", self.params.get("patch_y0_km", y)))
        sx = max(float(self.params.get("patch_sigma_x_km", 500.0)), 1e-6)
        sy = max(float(self.params.get("patch_sigma_y_km", 500.0)), 1e-6)
        patch = np.exp(-0.5 * (((x - x0) / sx) ** 2 + ((y - y0) / sy) ** 2))
        fluct = 1.0 + float(self.params.get("delta_N_N", 0.0)) * self._field.value(x, y, z, t)
        return float(max(layer * patch * max(fluct, 0.0) - Ne_bg * 0.0, 0.0))

    def metadata(self) -> dict[str, Any]:
        l_perp_m = float(self.params.get("L_perp_m", self.params.get("L1_m", 300.0)))
        l_para_m = float(self.params.get("L_para_m", 10000.0))
        meta = super().metadata()
        meta.update({
            "foEs_MHz": float(self.params.get("foEs_MHz", 5.0)),
            "Ne_peak_m3": ne_from_plasma_frequency_mhz(float(self.params.get("foEs_MHz", 5.0))),
            "h_Es_km": float(self.params.get("h_Es_km", 110.0)),
            "delta_h_km": float(self.params.get("delta_h_m", 115.0)) / 1000.0,
            "L_perp_m": l_perp_m,
            "L_para_m": l_para_m,
            "aspect_sensitivity": "Booker-Gordon field-aligned Gaussian factor",
            "anisotropy_ratio": l_para_m / max(l_perp_m, 1e-12),
            "b_hat": [float(v) for v in self._b_hat],
        })
        return meta


def build_components(component_specs: list[dict], geomag: dict | None = None,
                     seed: int | None = None) -> list[BaseIonosphere3DComponent]:
    out: list[BaseIonosphere3DComponent] = []
    for i, spec in enumerate(component_specs):
        kind = str(spec.get("type", spec.get("component_type", ""))).lower()
        params = dict(spec.get("params", spec))
        comp_seed = params.get("seed", seed if seed is None else int(seed) + i)
        if kind in ("tid", "tid3d"):
            out.append(TID3DComponent(params, geomag, comp_seed))
        elif kind in ("bubble", "bubble3d", "field_aligned_bubble"):
            out.append(Bubble3DComponent(params, geomag, comp_seed))
        elif kind in ("spreadf", "spread-f", "spreadf3d"):
            out.append(SpreadF3DComponent(params, geomag, comp_seed))
        elif kind in ("es", "es3d", "sporadic_e"):
            out.append(Es3DComponent(params, geomag, comp_seed))
        else:
            raise ValueError("unknown ionosphere component type: {}".format(kind))
    return out
