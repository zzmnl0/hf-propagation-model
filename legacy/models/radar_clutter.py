"""Single-station OTH radar clutter physics helpers.

The module keeps surface clutter separate from ray tracing.  A propagation mode
provides the one-way group path, delay, grazing geometry, and illuminated tube
area; this module supplies sigma0, Bragg/Doppler diagnostics, confidence, and a
monostatic radar-equation power estimate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from config import C_MS, FREQ_MHZ, GT, GR, PT_W
from utils import radar_equation_W, to_dBW


_SURFACE_SIGMA0_DB = {
    "sea": -20.0,
    "land": -25.0,
    "ice": -30.0,
    "desert": -28.0,
    "custom": -20.0,
}


@dataclass(frozen=True)
class RadarClutterConfig:
    """Physical inputs for single-station monostatic clutter."""

    surface_type: str = "sea"
    sigma0_dB: float | None = None
    radial_velocity_ms: float = 0.0
    grazing_angle_deg: float | None = None
    clutter_class: str = "surface_bragg"


class RadarClutterModel:
    """Compute monostatic OTH surface clutter diagnostics and power."""

    def __init__(self, freq_MHz: float = FREQ_MHZ, params: dict | None = None):
        self.freq_MHz = float(freq_MHz)
        self.params = dict(params or {})
        self.wavelength_m = C_MS / (self.freq_MHz * 1e6)
        self.k0_rad_m = 2.0 * math.pi / self.wavelength_m

    def sigma0_linear(
        self,
        surface_type: str = "sea",
        sigma0_dB: float | None = None,
        grazing_angle_deg: float | None = None,
    ) -> float:
        """Return normalized RCS sigma0 [linear], including grazing weighting."""
        s_type = str(surface_type or "sea").lower()
        base_dB = _SURFACE_SIGMA0_DB.get(s_type, _SURFACE_SIGMA0_DB["custom"])
        if sigma0_dB is not None:
            base_dB = float(sigma0_dB)
        sigma0 = 10.0 ** (base_dB / 10.0)

        if grazing_angle_deg is None:
            return float(max(sigma0, 0.0))
        grazing = max(float(grazing_angle_deg), 0.0)
        # Surface illuminated projection vanishes at zero grazing and tends to
        # the configured sigma0 at normal incidence.
        weight = math.sin(math.radians(min(grazing, 90.0)))
        return float(max(sigma0 * weight, 0.0))

    def bragg_wavenumber_rad_m(self) -> float:
        """Backscatter Bragg surface wavenumber |q| = 2 k0 [rad/m]."""
        return float(2.0 * self.k0_rad_m)

    def bragg_wavelength_m(self) -> float:
        """Surface Bragg wavelength lambda_B = pi/k0 = lambda_radar/2 [m]."""
        return float(self.wavelength_m / 2.0)

    def doppler_hz(self, radial_velocity_ms: float = 0.0) -> float:
        """Monostatic Doppler shift fD = 2 v_r / lambda [Hz]."""
        return float(2.0 * float(radial_velocity_ms) / self.wavelength_m)

    def sea_bragg_doppler_hz(self, radial_velocity_ms: float = 0.0) -> tuple[float, float]:
        """Deep-water first-order Bragg Doppler pair for monostatic sea clutter.

        Barrick first-order HF sea clutter peaks occur at +/-f_B, shifted by
        surface radial current.  Deep-water dispersion gives
        f_B = sqrt(g k_B) / (2*pi), k_B = 4*pi/lambda.
        """
        g_ms2 = 9.80665
        f_bragg = math.sqrt(g_ms2 * self.bragg_wavenumber_rad_m()) / (2.0 * math.pi)
        f_current = self.doppler_hz(radial_velocity_ms)
        return float(f_current - f_bragg), float(f_current + f_bragg)

    def estimate_grazing_angle_deg(self, mode: dict) -> float:
        """Estimate ground grazing angle from final ray segment or launch beta."""
        explicit = mode.get("grazing_angle_deg")
        if explicit is not None:
            return float(explicit)

        pts = mode.get("points")
        if pts is not None:
            arr = np.asarray(pts, dtype=float)
            if arr.ndim == 2 and arr.shape[0] >= 2:
                p0 = arr[-2]
                p1 = arr[-1]
                dx = float(p1[0] - p0[0])
                dz = float(p1[-1] - p0[-1])
                if abs(dx) + abs(dz) > 0.0:
                    return abs(float(np.degrees(np.arctan2(dz, dx))))

        return abs(float(mode.get("beta_deg", 0.0)))

    def mode_confidence(self, mode: dict, sigma0: float, grazing_angle_deg: float) -> float:
        """Return a bounded confidence diagnostic for clutter attribution."""
        if sigma0 <= 0.0:
            return 0.0
        if float(mode.get("group_path_km", 0.0)) <= 0.0:
            return 0.0
        area = float(mode.get("A_tube_km2", mode.get("scatter_area_km2", 0.0)))
        if area <= 0.0:
            return 0.0
        g = max(min(float(grazing_angle_deg), 90.0), 0.0)
        return float(max(0.05, min(1.0, math.sin(math.radians(g)) ** 0.5)))

    def evaluate_mode(
        self,
        mode: dict,
        Pt_W: float = PT_W,
        Gt: float = GT,
        Gr: float = GR,
        params: dict | None = None,
    ) -> dict:
        """Return clutter diagnostics and radar-equation received power."""
        cfg = dict(self.params)
        if params:
            cfg.update(params)

        surface_type = str(cfg.get("surface_type", "sea"))
        radial_velocity_ms = float(cfg.get("radial_velocity_ms", 0.0))
        grazing = cfg.get("grazing_angle_deg")
        grazing_angle_deg = (
            self.estimate_grazing_angle_deg(mode) if grazing is None else float(grazing)
        )
        sigma0 = self.sigma0_linear(
            surface_type=surface_type,
            sigma0_dB=cfg.get("sigma0_dB"),
            grazing_angle_deg=grazing_angle_deg,
        )

        area_km2 = float(mode.get("A_tube_km2", mode.get("scatter_area_km2", 0.0)))
        sigma_eff_m2 = sigma0 * area_km2 * 1e6
        group_path_km = float(mode.get("group_path_km", 0.0))
        Pr_W = 0.0
        if sigma_eff_m2 > 0.0 and group_path_km > 0.0:
            Pr_W = radar_equation_W(Pt_W, Gt, Gr, self.freq_MHz, group_path_km, sigma_eff_m2)

        tau_oneway = float(mode.get("tau_ms", 0.0))
        return {
            "clutter_enabled": True,
            "clutter_class": str(cfg.get("clutter_class", "surface_bragg")),
            "surface_type": surface_type,
            "sigma0_linear": float(sigma0),
            "sigma0_dB": to_dBW(sigma0),
            "sigma_eff_m2": float(sigma_eff_m2),
            "grazing_angle_deg": float(grazing_angle_deg),
            "bragg_wavenumber_rad_m": self.bragg_wavenumber_rad_m(),
            "bragg_wavelength_m": self.bragg_wavelength_m(),
            "radial_velocity_ms": radial_velocity_ms,
            "doppler_hz": self.doppler_hz(radial_velocity_ms),
            "range_km": group_path_km,
            "tau_oneway_ms": tau_oneway,
            "tau_2way_ms": 2.0 * tau_oneway,
            "clutter_Pr_W": float(Pr_W),
            "clutter_Pr_dBW": to_dBW(Pr_W),
            "mode_confidence": self.mode_confidence(mode, sigma0, grazing_angle_deg),
        }


def apply_clutter_to_mode(
    mode: dict,
    freq_MHz: float = FREQ_MHZ,
    params: dict | None = None,
    Pt_W: float = PT_W,
    Gt: float = GT,
    Gr: float = GR,
    apply_to_power: bool = False,
) -> dict:
    """Attach clutter diagnostics; optionally replace mode power."""
    model = RadarClutterModel(freq_MHz=freq_MHz, params=params)
    diag = model.evaluate_mode(mode, Pt_W=Pt_W, Gt=Gt, Gr=Gr)
    out = dict(mode)
    out.update(diag)
    if apply_to_power:
        out["Pr_W"] = diag["clutter_Pr_W"]
        out["Pr_dBW"] = diag["clutter_Pr_dBW"]
    return out


def build_surface_clutter_pd_spectrum(
    modes: list,
    freq_MHz: float = FREQ_MHZ,
    params: dict | None = None,
    tau_res_ms: float = 0.02,
    gate_span_km: float = 220.0,
    n_range_bins: int = 49,
    doppler_bins_hz: int = 81,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Build physically structured OTH surface-clutter P-D projection.

    The spectrum is a power-delay projection of a monostatic range-Doppler
    clutter surface.  For each flux-tube mode, neighbouring ground range gates
    contribute with R^-4 propagation, sigma0*area backscatter, grazing/focusing
    weighting, and first-order Bragg Doppler doublets.  Projecting this
    range-Doppler surface onto delay produces the peak-bearing curved envelope
    expected for OTH surface clutter, instead of a few isolated mode spikes.
    """
    if not modes:
        return np.array([0.0]), np.array([0.0]), {}

    cfg = dict(params or {})
    model = RadarClutterModel(freq_MHz=freq_MHz, params=cfg)
    surface_type = str(cfg.get("surface_type", "sea"))
    sigma0_dB = cfg.get("sigma0_dB")
    radial_velocity = float(cfg.get("radial_velocity_ms", 0.0))
    gate_span_km = max(float(cfg.get("gate_span_km", gate_span_km)), 1.0)
    n_range_bins = max(int(cfg.get("n_range_bins", n_range_bins)), 5)
    doppler_bins_hz = max(int(cfg.get("doppler_bins_hz", doppler_bins_hz)), 11)

    records = []
    offsets = np.linspace(-0.5 * gate_span_km, 0.5 * gate_span_km, n_range_bins)
    for mode_index, m in enumerate(modes):
        tau_c = float(m.get("tau_2way_ms", 2.0 * float(m.get("tau_ms", 0.0))))
        R_c = max(float(m.get("range_km", m.get("group_path_km", 0.0))), 1e-6)
        base_area = max(float(m.get("A_tube_km2", m.get("scatter_area_km2", 0.0))), 0.0)
        if tau_c <= 0.0 or base_area <= 0.0:
            continue

        tau_per_km = tau_c / R_c
        grazing0 = float(m.get("grazing_angle_deg", model.estimate_grazing_angle_deg(m)))
        focus = max(float(m.get("F_focus_3d", m.get("F_focus", 1.0))), 0.0)
        base_Pr = max(float(m.get("clutter_Pr_W", m.get("Pr_W", 0.0))), 0.0)
        gate_sigma_km = max(0.22 * gate_span_km, 1.0)
        bragg_pair = model.sea_bragg_doppler_hz(radial_velocity)
        doppler_sigma = max(float(cfg.get("bragg_width_hz", 0.08)), 0.01)
        d_min = min(bragg_pair) - 4.0 * doppler_sigma
        d_max = max(bragg_pair) + 4.0 * doppler_sigma
        dopplers = np.linspace(d_min, d_max, doppler_bins_hz)
        doppler_kernel = (
            np.exp(-0.5 * ((dopplers - bragg_pair[0]) / doppler_sigma) ** 2)
            + np.exp(-0.5 * ((dopplers - bragg_pair[1]) / doppler_sigma) ** 2)
        )
        doppler_kernel /= max(float(doppler_kernel.sum()), 1e-30)

        for dx in offsets:
            R = max(R_c + float(dx), 1e-6)
            tau = tau_c + tau_per_km * float(dx)
            grazing = max(0.0, min(90.0, grazing0 + 0.02 * float(dx)))
            sigma0 = model.sigma0_linear(surface_type, sigma0_dB, grazing)
            if sigma0 <= 0.0:
                continue
            range_weight = math.exp(-0.5 * (float(dx) / gate_sigma_km) ** 2)
            area_gate = base_area / n_range_bins
            range_loss_weight = (R_c / R) ** 4
            Pr_gate = base_Pr * range_weight * range_loss_weight
            if base_Pr <= 0.0:
                sigma_eff = sigma0 * area_gate * 1e6 * max(focus, 0.05)
                Pr_gate = radar_equation_W(PT_W, GT, GR, freq_MHz, R, sigma_eff)
            for doppler, dw in zip(dopplers, doppler_kernel):
                records.append((tau, Pr_gate * float(dw) / n_range_bins, float(doppler), mode_index))

    if not records:
        taus = [float(m.get("tau_2way_ms", 2.0 * float(m.get("tau_ms", 0.0)))) for m in modes]
        tau_axis = np.arange(max(0.0, min(taus) - 2.0), max(taus) + 2.0, tau_res_ms)
        return tau_axis, np.zeros(len(tau_axis)), {}

    tau_values = np.array([r[0] for r in records], dtype=float)
    tau_axis = np.arange(max(0.0, float(tau_values.min()) - 1.0),
                         float(tau_values.max()) + 1.0, tau_res_ms)
    pd = np.zeros(len(tau_axis))
    tau_sigma = max(float(cfg.get("delay_width_ms", 0.12)), tau_res_ms)
    for tau, power, _, _ in records:
        pd += float(power) * np.exp(-0.5 * ((tau_axis - float(tau)) / tau_sigma) ** 2)

    scale = sum(max(float(m.get("clutter_Pr_W", m.get("Pr_W", 0.0))), 0.0) for m in modes)
    pd_sum = float(pd.sum())
    if scale > 0.0 and pd_sum > 0.0:
        pd *= scale / pd_sum

    meta = {
        "clutter_pd_model": "range_doppler_bragg_projection",
        "surface_type": surface_type,
        "bragg_doppler_minus_hz": model.sea_bragg_doppler_hz(radial_velocity)[0],
        "bragg_doppler_plus_hz": model.sea_bragg_doppler_hz(radial_velocity)[1],
        "range_gate_span_km": gate_span_km,
        "n_range_bins": n_range_bins,
        "n_doppler_bins": doppler_bins_hz,
    }
    return tau_axis, pd, meta
