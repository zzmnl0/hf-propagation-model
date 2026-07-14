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
from utils import free_space_loss_dB, radar_equation_W, to_dBW
from .hybrid_model import build_pd_spectrum, identify_main_mode
from .ne_provider import BaseNeProvider
from .point_to_point_3d import annotate_mode_causes_3d, find_all_rays_3d
from .ray_tracer_3d import RefractiveIndex3D, shoot_rays_fan_3d


@dataclass
class UnifiedPropagationRequest:
    scenario: str
    provider: BaseNeProvider
    tx_pos: tuple = cfg.TX_POS
    rx_pos: tuple = cfg.RX_POS
    freq_MHz: float = cfg.FREQ_MHZ
    x_km: Any = field(default_factory=lambda: cfg.BG_X)
    z_km: Any = field(default_factory=lambda: cfg.BG_Z)
    rt_config: dict = field(default_factory=lambda: dict(cfg.RT_3D))
    p2p_config: dict = field(default_factory=lambda: dict(cfg.P2P_3D))
    solver: str = "iso3d_legacy"
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
        m["solver_is_legacy"] = bool(request.solver == "iso3d_legacy")
        if request.solver == "ah3d":
            m.setdefault("ox_solver", "ah3d_hamiltonian")
            m.setdefault("group_index_model", "ah3d_fd_dn_df")
            m.setdefault("group_path_definition", "integral_n_g_ds")
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
            from .absorption_collision import complex_absorption_path_dB
            path = m.get("trajectory_3d", m.get("points", []))
            sampler = lambda x, y, z: n3d.Ne_at(x, y, z)
            abs_diag = complex_absorption_path_dB(
                path, request.freq_MHz, sampler, complex_cfg,
                two_way=bool(radar_mode), wave_mode=m.get("wave_mode", "O"),
                geomag=cfg.GEOMAG)
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


def _apply_ray_pe_reference(modes: list, Ne_2d: np.ndarray,
                            irregularity_params: dict | None) -> None:
    if not modes:
        return
    from .hybrid_ray_pe import apply_ray_pe_coupling_3d
    from .pe_propagator import PEPropagator
    pe_cfg = dict(cfg.UPDATE2.get("PE_COUPLING", {}))
    pe_cfg.update({
        "enable": True,
        "apply_to_power": bool(cfg.UPDATE2.get("PE_COUPLING", {}).get("apply_to_power", False)),
        "mode": "local_2p5d",
        "q_boundary": 1.0,
    })
    pe = PEPropagator(freq_MHz=cfg.FREQ_MHZ)
    params = dict(irregularity_params or {})
    for m in modes:
        Pr_in = float(m.get("Pr_W", 0.0))
        Pr_out, diag = apply_ray_pe_coupling_3d(
            pe, m, Ne_2d, cfg.BG_X, cfg.BG_Z, Pr_in,
            pe.k0, params, pe_cfg)
        m.update(diag)
        m["pe_solver_dimensionality"] = "local_2p5d_reference"
        m["pe_validity"] = "reference_only"
        if bool(diag.get("pe_apply_to_power", False)):
            m["Pr_W"] = float(Pr_out)
            m["Pr_dBW"] = to_dBW(m["Pr_W"])
        m["delta_tau_ms"] = (
            float(m.get("delta_tau_ms", 0.0)) +
            float(diag.get("pe_delta_tau_ms", 0.0)))


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
    provider_meta = request.provider.get_metadata()
    Ne_2d = request.provider.get_Ne_grid(request.x_km, request.z_km)
    n3d = RefractiveIndex3D(
        Ne_2d, request.x_km, request.z_km, request.freq_MHz,
        earth_flat=True, ne_provider=request.provider)
    rays3d = []
    if request.make_ray_fan:
        if request.solver == "ah3d":
            from .ray_tracer_ah_3d import RefractiveIndexAH3D, shoot_rays_fan_ah_3d
            rays3d = []
            wave_modes = (("O", "X") if cfg.UPDATE2.get("RAY_FAN_INCLUDE_X", False)
                          else ("O",))
            for wave_mode in wave_modes:
                n_ah = RefractiveIndexAH3D(
                    n3d, request.freq_MHz, wave_mode=wave_mode, geomag=cfg.GEOMAG,
                    integration_profile="smoke")
                rays3d.extend(shoot_rays_fan_ah_3d(
                    request.tx_pos, n_ah, request.freq_MHz, request.rt_config))
        else:
            rays3d = shoot_rays_fan_3d(
                request.tx_pos, n3d, request.freq_MHz, request.rt_config)

    if not request.make_p2p:
        modes = []
    elif request.solver == "ah3d":
        from .point_to_point_ah_3d import find_all_rays_ah_3d
        modes = find_all_rays_ah_3d(
            request.tx_pos, request.rx_pos, n3d, request.freq_MHz,
            request.p2p_config, wave_mode="both", geomag=cfg.GEOMAG,
            rays_fan=rays3d)
        for m in modes:
            m["ox_solver"] = "ah3d_hamiltonian"
            m["physics_basis"] = "3d_appleton_hartree_hamiltonian"
            m["freq_MHz"] = float(request.freq_MHz)
    else:
        modes = find_all_rays_3d(
            request.tx_pos, request.rx_pos, n3d, request.freq_MHz,
            request.p2p_config)
        for m in modes:
            m["ox_solver"] = "iso3d_legacy"
            m["physics_basis"] = m.get("physics_basis", "3d_haselgrove_isotropic_legacy")
            m["group_index_model"] = m.get("group_index_model", "cold_plasma_ng_approx_legacy")
            m["freq_MHz"] = float(request.freq_MHz)

    if modes:
        _add_power(modes, n3d, request)
        if request.apply_ray_pe_reference:
            _apply_ray_pe_reference(modes, Ne_2d, request.pe_irregularity_params)
        _add_scatter_diagnostics(modes, n3d, request)
        _apply_p3_energy_closure(modes, request.apply_scatter_to_pd)
        annotate_mode_causes_3d(modes, request.p2p_config,
                                environment=request.environment)
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
    if request.apply_ray_pe_reference:
        warnings.append("Ray-PE is local_2p5d_reference, not full 3-D PE/full-wave.")
    if request.solver == "iso3d_legacy":
        warnings.append("iso3d_legacy is a baseline/performance fallback, not the P2 AH3D research mainline.")
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
