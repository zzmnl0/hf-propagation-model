"""Unified stage-1 3-D reflection/scattering propagation pipeline.

P0 scope: centralize the existing 3-D propagation chain and output semantics.
It does not claim new physics beyond the existing ray, absorption, PE, and
scatter modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import config as cfg
from utils import (enu_basis, free_space_loss_dB, haversine_km,
                   radar_equation_W, to_dBW)
from .pd_spectrum import build_pd_spectrum, identify_main_mode
from .ne_provider_geocentric import GeocentricNeProvider3D
from .point_to_point_ah_3d_ecef import find_seeded_rays_ah3d_ecef
from .ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF,
    trace_ray_ah3d_ecef,
)


@dataclass
class UnifiedPropagationRequest:
    scenario: str
    provider: GeocentricNeProvider3D
    tx_geocentric: tuple = cfg.TX_GEOCENTRIC
    rx_geocentric: tuple = cfg.RX_GEOCENTRIC
    bearing_deg: float = cfg.LINK_BEARING_DEG
    earth_radius_km: float = cfg.RE_KM
    freq_MHz: float = cfg.FREQ_MHZ
    x_km: Any = field(default_factory=lambda: cfg.BG_X)
    z_km: Any = field(default_factory=lambda: cfg.BG_Z)
    rt_config: dict = field(default_factory=lambda: dict(cfg.RT_3D))
    p2p_config: dict = field(default_factory=lambda: dict(cfg.P2P_3D))
    solver: str = "ah3d_ecef"
    observable_type: str = "communication_oneway"
    pd_method: str = "gaussian"
    pd_delay_key: str = "tau_ms"
    apply_scatter_to_pd: bool = False
    apply_ray_pe_reference: bool = False
    pe_irregularity_params: dict | None = None
    environment: str | None = None
    metadata: dict = field(default_factory=dict)
    make_ray_fan: bool = True
    make_p2p: bool = True
    make_pd: bool = True


@dataclass
class UnifiedPropagationResult:
    scenario: str
    provider_metadata: dict
    Ne_2d: np.ndarray
    n3d: Any
    rays3d: list
    modes: list
    tau_axis: np.ndarray
    pd_W: np.ndarray
    main_mode: dict | None
    mode_summary: dict
    traceability: dict
    warnings: list[str]


def _annotate_mode_causes(modes: list, p2p_params: dict,
                          environment: str | None) -> None:
    """Attach coordinate-independent branch semantics to formal ECEF modes."""
    low_max = float(p2p_params.get("launch_low_max_deg", 25.0))
    high_min = float(p2p_params.get("launch_high_min_deg", 55.0))
    seed = str(p2p_params.get("seed_source", "internal_scan"))
    env = str(environment or "").lower()
    for mode in modes:
        beta = float(mode.get("elevation_tx_deg", mode.get("beta_deg_tx", 0.0)))
        launch = ("low_elevation" if beta <= low_max else
                  "high_elevation" if beta >= high_min else "medium_elevation")
        height = float(mode.get("apex_height_km", mode.get("h_reflect_km", 0.0)))
        label = str(mode.get("label", ""))
        if label.startswith("Es") or height < float(cfg.MODE["h_Es_km"]):
            layer = "Es"
        elif height < float(cfg.MODE["h_E_km"]):
            layer = "E"
        else:
            layer = "F1" if height < 300.0 else "F2"
        wave_mode = str(mode.get("wave_mode", "iso"))
        if wave_mode in ("O", "X"):
            mechanism = "magnetoionic_split"
        elif bool(mode.get("pe_coupled", False)):
            mechanism = "ray_pe_coupled"
        elif layer == "Es":
            mechanism = "layered_reflection"
        elif env in ("bubble", "spreadf", "full") and mode.get("scatter_method"):
            mechanism = "irregularity_scatter"
        else:
            mechanism = "elevation_branch"
        summary = f"{layer} {launch.replace('_', '-')} branch"
        if mechanism in ("magnetoionic_split", "irregularity_scatter", "ray_pe_coupled"):
            summary += f" ({mechanism})"
        mode.update({
            "cause_layer": layer, "cause_launch": launch,
            "cause_mechanism": mechanism,
            "cause_summary": summary,
            "p2p_seed_source": seed,
        })


def _fixed_geomagnetic_direction_ecef(tx_geocentric) -> np.ndarray:
    """Map the configured local declination/dip direction into TX ECEF."""
    decl = np.radians(float(cfg.GEOMAG.get("decl_deg", 0.0)))
    dip = np.radians(float(cfg.GEOMAG.get("dip_deg", 0.0)))
    east, north, up = enu_basis(
        float(tx_geocentric[0]), float(tx_geocentric[1]), cfg.RE_KM).T
    direction = (np.cos(dip) * np.sin(decl) * east
                 + np.cos(dip) * np.cos(decl) * north
                 - np.sin(dip) * up)
    return direction / np.linalg.norm(direction)


def _ah_model(request: UnifiedPropagationRequest, wave_mode: str,
              profile: str) -> RefractiveIndexAH3DECEF:
    return RefractiveIndexAH3DECEF(
        request.provider, earth_radius_km=request.earth_radius_km,
        freq_MHz=request.freq_MHz, wave_mode=wave_mode,
        gyro_MHz=float(cfg.GEOMAG.get("fH_MHz", 0.0)),
        b_ecef_hat=_fixed_geomagnetic_direction_ecef(request.tx_geocentric),
        integration_profile=profile)


def _trace_fan(request: UnifiedPropagationRequest, wave_mode: str,
               profile: str = "smoke") -> list:
    model = _ah_model(request, wave_mode, profile)
    elevations = np.linspace(
        float(request.rt_config["beta_min"]),
        float(request.rt_config["beta_max"]),
        int(request.rt_config["n_fan"]))
    max_path_km = float(request.rt_config.get(
        "max_path_km",
        request.rt_config.get("max_steps", 4000)
        * request.rt_config.get("ds_max_km", 2.0)))
    rays = []
    for elevation in elevations:
        try:
            rays.append(trace_ray_ah3d_ecef(
                model, tx_lat_deg=float(request.tx_geocentric[0]),
                tx_lon_deg=float(request.tx_geocentric[1]),
                tx_height_km=float(request.tx_geocentric[2]),
                elevation_deg=float(elevation),
                bearing_deg=float(request.bearing_deg),
                h_stop_km=float(request.rt_config["z_stop_km"]),
                max_path_km=max_path_km, profile=profile))
        except (ValueError, RuntimeError, FloatingPointError):
            continue
    return rays


def _nearest_ground_seeds(rays: list, rx_geocentric, count: int) -> list:
    ranked = []
    for ray in rays:
        landing = ray.get("landing_geocentric_lat_lon_height")
        if ray.get("termination_status") != "ground" or landing is None:
            continue
        error = haversine_km(
            float(landing[0]), float(landing[1]),
            float(rx_geocentric[0]), float(rx_geocentric[1]))
        ranked.append((error, ray))
    ranked.sort(key=lambda item: item[0])
    return [{
        "seed_id": "ecef_fan_{:02d}".format(index),
        "elevation_deg": float(ray["elevation_tx_deg"]),
        "bearing_deg": float(ray["bearing_tx_deg"]),
        "semantics": "new_coordinate_coarse_fan_seed",
    } for index, (_, ray) in enumerate(ranked[:max(int(count), 1)])]


def _provider_components(meta: dict) -> str:
    components = meta.get("components", [])
    if components:
        names = [
            str(c.get("component_type", "component"))
            for c in components
            if c.get("enabled", False)
        ]
        if names:
            return "+".join(names)
    parts = []
    for key in ("tid_enabled", "bubble_enabled", "spreadf_enabled", "es_enabled",
                "spreadf_diagnostic_only", "es_diagnostic_only"):
        if meta.get(key):
            parts.append(key)
    if not parts:
        parts.append("background")
    return "+".join(parts)


def _feature_gates() -> str:
    u2 = cfg.UPDATE2
    gates = {
        "AH3D_ENABLE": bool(u2.get("AH3D_ENABLE", False)),
        "AH3D_SOLVER": u2.get("AH3D_SOLVER", "v2"),
        "AH3D_PROFILE": u2.get("AH3D_INTEGRATION_PROFILE", "research"),
        "OX_MODE": u2.get("OX_MODE", "postproc"),
        "FAN_X": bool(u2.get("RAY_FAN_INCLUDE_X", False)),
        "ABS_COMPLEX": bool(u2.get("ABSORPTION_COMPLEX", {}).get("enable", False)),
        "PE": bool(u2.get("PE_COUPLING", {}).get("enable", False)),
        "SCATTER_PD": bool(u2.get("SCATTER_PD", {}).get("apply_to_pd", False)),
    }
    return ";".join("{}={}".format(k, v) for k, v in gates.items())


def _tag_common_semantics(modes: list, request: UnifiedPropagationRequest,
                          provider_meta: dict) -> None:
    components = _provider_components(provider_meta)
    solver_chain = "{}->power->scatter_diag->pd".format(request.solver)
    one_or_two = "two_way" if request.pd_delay_key == "tau_2way_ms" else "one_way"
    for m in modes:
        m["provider_name"] = provider_meta.get("name", "")
        m["provider_components"] = components
        m["provider_hash"] = provider_meta.get("provider_hash", "")
        m["provider_config_hash"] = provider_meta.get("config_hash", "")
        m["solver_chain"] = solver_chain
        m["feature_gates"] = _feature_gates()
        m["pd_kernel"] = "scatter_coupled" if request.apply_scatter_to_pd else request.pd_method
        m["observable_type"] = request.observable_type
        m["delay_key"] = request.pd_delay_key
        m["oneway_or_twoway"] = one_or_two
        m["ray_pe_dimensionality"] = m.get("pe_solver_dimensionality", "off")
        m["scatter_status"] = m.get("scatter_method", "diagnostic/off")
        m["calibration_status"] = "reference_uncalibrated"
        m["solver_is_legacy"] = False
        m.setdefault("ox_solver", "ah3d_ecef_hamiltonian")
        m.setdefault("group_index_model", "ah3d_dispersion_Dk_Domega")
        m.setdefault("group_path_definition", "c*t_group")
        m.setdefault("earth_model", "sphere")
        m.setdefault("coordinate_model", "geocentric_spherical_ECEF_state")
        m.setdefault("coordinate_semantics", (
            "input=geocentric_lat_lon_height;state=ECEF;plot=derived_ground_range_height"))
        m.setdefault("plot_coordinate_semantics",
                     "derived_great_circle_ground_range_physical_height")
        m.setdefault("theory_reference", "Appleton-Hartree; Haselgrove Hamiltonian; Jones-Stephenson 3D ray tracing")


def _tag_communication_modes(modes: list) -> None:
    for m in modes:
        tau = float(m.get("tau_ms", 0.0))
        m["tau_oneway_ms"] = tau
        m["tau_2way_ms"] = None
        m["range_km"] = float(m.get("group_path_km", 0.0))
        m["radar_observable"] = "none"
        m.setdefault("physics_basis", "3d_haselgrove_friis")


def _tag_point_target_modes(modes: list) -> None:
    for m in modes:
        tau = float(m.get("tau_ms", 0.0))
        m["tau_oneway_ms"] = tau
        m["tau_2way_ms"] = 2.0 * tau
        m["range_km"] = float(m.get("group_path_km", 0.0))
        m["radar_observable"] = "two_way_point_target"
        m["radar_equation"] = "monostatic_point_target_R4"
        m["physics_basis"] = "Skolnik_2008_point_target_radar_equation"


def _add_power(modes: list, n3d, request: UnifiedPropagationRequest) -> None:
    radar_mode = request.observable_type == "radar_point_target_twoway"
    for m in modes:
        gp = float(m.get("group_path_km", 0.0))
        if radar_mode:
            sigma = float(request.metadata.get("sigma_rcs_m2",
                                               cfg.RADAR.get("sigma_rcs_m2", 1.0)))
            m["Pr_W"] = radar_equation_W(cfg.PT_W, cfg.GT, cfg.GR,
                                          request.freq_MHz, gp, sigma)
        else:
            loss_dB = free_space_loss_dB(gp, request.freq_MHz)
            m["Pr_W"] = cfg.PT_W * cfg.GT * cfg.GR * 10.0 ** (-loss_dB / 10.0)

        abs_oneway_dB = 0.0
        abs_total_dB = 0.0
        abs_model = "none"
        abs_diag = {"abs_status": "disabled"}
        complex_cfg = cfg.UPDATE2.get("ABSORPTION_COMPLEX", {})
        if complex_cfg.get("enable", False):
            from .absorption_collision import complex_absorption_ecef_ray_dB
            abs_diag = complex_absorption_ecef_ray_dB(
                m, n3d, complex_cfg, two_way=bool(radar_mode))
            abs_oneway_dB = float(abs_diag.get("abs_complex_dB", 0.0))
            abs_total_dB = float(abs_diag.get("abs_total_dB", abs_oneway_dB))
            abs_model = "appleton_hartree_collision_path"
            m["Pr_W"] *= 10.0 ** (-abs_total_dB / 10.0)

        m["abs_dB"] = float(abs_oneway_dB)
        m["abs_total_dB"] = float(abs_total_dB)
        m["abs_oneway_dB"] = float(abs_oneway_dB)
        m["abs_model"] = abs_model
        m["abs_complex_dB"] = float(abs_diag.get("abs_complex_dB", 0.0))
        m["abs_status"] = abs_diag.get("abs_status", "")
        m["collision_model"] = abs_diag.get("collision_model", "")
        m["abs_index_model"] = abs_diag.get("abs_index_model", "")
        m["absorption_applied_once"] = bool(abs_model != "none")
        m["legacy_d_layer_absorption_applied"] = False
        m["Pr_dBW"] = to_dBW(m["Pr_W"])
        m["delta_tau_ms"] = float(m.get("delta_tau_ms", 0.0))

    if radar_mode:
        _tag_point_target_modes(modes)
    else:
        _tag_communication_modes(modes)


def _add_scatter_diagnostics(modes: list, n3d, request: UnifiedPropagationRequest) -> None:
    if not cfg.SCATTER.get("enable", True):
        return
    from .scatter_coupler import apply_scatter_regime_3d
    apply_scatter_regime_3d(
        modes, n3d, cfg.SCATTER, request.freq_MHz, cfg.PT_W, cfg.GT, cfg.GR,
        enable_mps=cfg.SCATTER.get("enable_mps", False))


def _apply_p3_energy_closure(modes: list, apply_scatter_to_pd: bool) -> None:
    for m in modes:
        before = float(m.get("Pr_before_scatter", m.get("Pr_W", 0.0)))
        after = float(m.get("Pr_W", before))
        m["Pr_before_scatter"] = before
        m["Pr_after_scatter"] = after
        if before > 0.0 and after > 0.0:
            m["scatter_loss_dB"] = float(-10.0 * np.log10(after / before))
            rel_err = abs(after - before) / before
        else:
            m["scatter_loss_dB"] = 0.0
            rel_err = 0.0
        m["energy_conservation_error"] = float(m.get("energy_conservation_error", rel_err))
        m.setdefault("kernel_normalization", 1.0 if apply_scatter_to_pd else 0.0)
        m.setdefault("scatter_power_applied", False)
        m.setdefault("scatter_solver_dimensionality", m.get("pe_solver_dimensionality", "diagnostic"))
        m.setdefault("scatter_validity", "diagnostic_only")


def _build_pd_for_delay(modes: list, request: UnifiedPropagationRequest):
    if request.apply_scatter_to_pd:
        old = cfg.UPDATE2["SCATTER_PD"].get("apply_to_pd", False)
        old_kernel = cfg.UPDATE2["SCATTER_PD"].get("kernel")
        try:
            cfg.UPDATE2["SCATTER_PD"]["apply_to_pd"] = True
            cfg.UPDATE2["SCATTER_PD"]["kernel"] = "physical"
            return _build_pd_for_delay(
                modes,
                UnifiedPropagationRequest(
                    **{**request.__dict__, "apply_scatter_to_pd": False}))
        finally:
            cfg.UPDATE2["SCATTER_PD"]["apply_to_pd"] = old
            cfg.UPDATE2["SCATTER_PD"]["kernel"] = old_kernel
    if request.pd_delay_key == "tau_ms":
        return build_pd_spectrum(modes, method=request.pd_method)
    shifted = []
    for m in modes:
        mm = dict(m)
        mm["tau_ms"] = float(m.get(request.pd_delay_key, m.get("tau_ms", 0.0)))
        if request.pd_delay_key == "tau_2way_ms":
            if "tau_spread_ms" in mm:
                mm["tau_spread_ms"] = 2.0 * float(mm.get("tau_spread_ms", 0.0))
            if "delta_tau_ms" in mm:
                mm["delta_tau_ms"] = 2.0 * float(mm.get("delta_tau_ms", 0.0))
        shifted.append(mm)
    return build_pd_spectrum(shifted, method=request.pd_method)


def run_unified_pipeline(request: UnifiedPropagationRequest) -> UnifiedPropagationResult:
    if request.solver != "ah3d_ecef":
        raise ValueError("C0R production pipeline accepts only solver='ah3d_ecef'")
    if not isinstance(request.provider, GeocentricNeProvider3D):
        raise TypeError("C0R production pipeline requires GeocentricNeProvider3D")
    if float(request.earth_radius_km) != float(cfg.EARTH_MODEL["radius_km"]):
        raise ValueError("request Earth radius differs from config EARTH_MODEL")
    if request.apply_ray_pe_reference:
        raise RuntimeError(
            "formal ECEF pipeline rejects archived local_2p5d PE; migrate it before use")
    provider_meta = request.provider.get_metadata()
    Ne_2d = request.provider.sample_link_plane(
        request.tx_geocentric, request.bearing_deg,
        request.x_km, request.z_km)
    n3d = request.provider
    fan_by_mode = {}
    rays3d = []
    if request.make_ray_fan:
        wave_modes = (("O", "X") if cfg.UPDATE2.get("RAY_FAN_INCLUDE_X", False)
                      else ("O",))
        for wave_mode in wave_modes:
            fan_by_mode[wave_mode] = _trace_fan(request, wave_mode)
            rays3d.extend(fan_by_mode[wave_mode])
        for ray in rays3d:
            ray["coordinate_semantics"] = (
                "input=geocentric_lat_lon_height;state=ECEF;plot=derived_ground_range_height")
            ray["plot_coordinate_semantics"] = (
                "derived_great_circle_ground_range_physical_height")
            ray["oneway_or_twoway"] = (
                "two_way" if request.pd_delay_key == "tau_2way_ms" else "one_way")
            ray["ray_pe_dimensionality"] = (
                "archived_local_2p5d_not_in_production_chain")
            ray["scatter_status"] = (
                "power_coupled" if request.apply_scatter_to_pd else "diagnostic_only")

    if not request.make_p2p:
        modes = []
    else:
        modes = []
        for wave_mode in ("O", "X"):
            candidate_fan = fan_by_mode.get(wave_mode)
            if candidate_fan is None:
                candidate_fan = _trace_fan(request, wave_mode)
            seeds = _nearest_ground_seeds(
                candidate_fan, request.rx_geocentric,
                request.p2p_config.get("candidate_nearest", 3))
            if not seeds:
                continue
            model = _ah_model(request, wave_mode, "reference")
            roots = find_seeded_rays_ah3d_ecef(
                model, tx_geocentric=request.tx_geocentric,
                rx_geocentric=request.rx_geocentric, seeds=seeds,
                h_stop_km=float(request.rt_config["z_stop_km"]),
                max_path_km=float(request.rt_config.get("max_path_km", 8000.0)),
                profile="reference", tol_km=float(request.p2p_config["tol_km"]),
                max_nfev=int(request.p2p_config["max_iter"]),
                d_elevation_deg=float(request.p2p_config["d_beta"]),
                d_bearing_deg=float(request.p2p_config["d_az"]),
                jacobian_scales=request.p2p_config.get("jacobian_scales", (0.5, 1.0, 2.0)),
                jacobian_scale_rtol=float(request.p2p_config.get("jacobian_scale_rtol", 1.0e-2)),
                jacobian_condition_max=float(request.p2p_config.get("jacobian_condition_max", 1.0e10)),
                dedup_launch_deg=float(request.p2p_config.get("clust_launch_deg", 0.02)),
                dedup_trajectory_rms_km=float(request.p2p_config.get("clust_trajectory_rms_km", 1.0)),
                dedup_apex_km=float(request.p2p_config.get("clust_h_km", 10.0)),
                dedup_tau_ms=float(request.p2p_config.get("clust_tau_ms", 0.05)))
            for index, m in enumerate(roots, start=1):
                m["label"] = "{}-ECEF-{:02d}".format(
                    wave_mode, len(modes) + index)
                m["h_reflect_km"] = float(m.get("apex_height_km", 0.0))
                m["beta_deg_tx"] = float(m["elevation_tx_deg"])
                m["az_deg_tx"] = float(m["bearing_tx_deg"])
                m["branch_status"] = "ok"
                m["ox_solver"] = "ah3d_ecef_hamiltonian"
                m["physics_basis"] = "3d_appleton_hartree_hamiltonian_ECEF"
                m["freq_MHz"] = float(request.freq_MHz)
            modes.extend(roots)

    if modes:
        _add_power(modes, n3d, request)
        _add_scatter_diagnostics(modes, n3d, request)
        _apply_p3_energy_closure(modes, request.apply_scatter_to_pd)
        _annotate_mode_causes(modes, request.p2p_config, request.environment)
        _tag_common_semantics(modes, request, provider_meta)

    if request.make_pd and modes:
        tau_axis, pd_W = _build_pd_for_delay(modes, request)
        main_mode, ranked, summary = identify_main_mode(tau_axis, pd_W, modes)
    else:
        tau_axis = np.asarray([], dtype=float)
        pd_W = np.asarray([], dtype=float)
        main_mode, ranked, summary = None, [], {}
    traceability = {
        "scenario": request.scenario,
        "earth_model": "sphere",
        "earth_radius_km": float(request.earth_radius_km),
        "coordinate_model": "geocentric_spherical_ECEF_state",
        "coordinate_semantics": (
            "TX_RX=lat_lon_height;launch=elevation_bearing;state=ECEF;"
            "plot=derived_great_circle_range_physical_height"),
        "oneway_or_twoway": (
            "two_way" if request.pd_delay_key == "tau_2way_ms" else "one_way"),
        "reflection_model": "AH3D_Appleton_Hartree_Hamiltonian_ECEF",
        "ray_pe_coupling": "off_archived_local_2p5d_not_in_production_chain",
        "scatter_power_coupling": (
            "on" if request.apply_scatter_to_pd else "diagnostic_only"),
        "provider_name": provider_meta.get("name", ""),
        "provider_components": _provider_components(provider_meta),
        "provider_hash": provider_meta.get("provider_hash", ""),
        "provider_config_hash": provider_meta.get("config_hash", ""),
        "provider_seed": provider_meta.get("seed"),
        "provider_t_s": provider_meta.get("t_s", 0.0),
        "solver_chain": modes[0].get("solver_chain", "") if modes else request.solver,
        "feature_gates": _feature_gates(),
        "ranked_modes": len(ranked),
        "ray_fan_wave_modes": sorted({str(r.get("wave_mode", "iso")) for r in rays3d}),
        "p2p_wave_modes": sorted({str(m.get("wave_mode", "iso")) for m in modes}),
    }
    warnings = []
    bad_status = sorted({
        str(m.get("branch_status", ""))
        for m in modes
        if str(m.get("branch_status", "")) not in ("", "ok", "near_cutoff")
    })
    if bad_status:
        warnings.append("AH3D rejected/non-mainline branch statuses present: {}".format(",".join(bad_status)))
    if "spreadf_diagnostic_only" in traceability.get("provider_components", ""):
        warnings.append("Spread-F provider component remains diagnostic in P0.")
    if "es_diagnostic_only" in traceability.get("provider_components", ""):
        warnings.append("Es provider component remains diagnostic in P0.")

    return UnifiedPropagationResult(
        scenario=request.scenario,
        provider_metadata=provider_meta,
        Ne_2d=Ne_2d,
        n3d=n3d,
        rays3d=rays3d,
        modes=modes,
        tau_axis=tau_axis,
        pd_W=pd_W,
        main_mode=main_mode,
        mode_summary=summary,
        traceability=traceability,
        warnings=warnings,
    )
