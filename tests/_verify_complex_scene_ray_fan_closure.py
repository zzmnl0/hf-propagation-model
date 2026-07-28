"""Verify and package the formal parameterized-scene O-wave ray fans."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config as cfg
from models.ionosphere_components_3d import geomag_unit_vector
from models.ne_provider_geocentric import GeocentricLinkGridNeProvider3D
from models.scenario_provider_factory import ScenarioProviderFactory


SCENES = ("baseline", "tid", "spreadf", "es", "bubble")
OUT = ROOT / "output2"
ASSETS = OUT / "report_assets" / "reflection_route_2026-07-26"
ASSETS.mkdir(parents=True, exist_ok=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sample_fields(factory):
    providers = {
        scene: factory.create(scene, seed=0, t_s=0.0, validation=True)
        for scene in SCENES
    }
    z_es = np.linspace(
        cfg.ES["h_Es_km"] - 0.6, cfg.ES["h_Es_km"] + 0.6, 121)
    z_plot = np.unique(np.r_[cfg.BG_Z, z_es])
    fields = {
        scene: provider.sample_link_plane(
            cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, z_plot)
        for scene, provider in providers.items()
    }
    return providers, fields, z_plot


def _scene_figures(fields, z_plot):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True)
    positive = np.concatenate([
        value[value > 0.0] for value in fields.values()])
    vmin = float(np.log10(max(np.percentile(positive, 0.2), 1.0)))
    vmax = float(np.log10(max(value.max() for value in fields.values())))
    image = None
    for ax, scene in zip(axes.flat, SCENES):
        image = ax.pcolormesh(
            cfg.BG_X, z_plot, np.log10(np.maximum(fields[scene].T, 1.0)),
            cmap="jet", shading="auto", alpha=0.5, vmin=vmin, vmax=vmax)
        ax.set_title(scene)
        ax.set_xlabel("Ground range (km)")
        ax.set_ylabel("Physical height (km)")
        if scene == "es":
            ax.axhline(cfg.ES["h_Es_km"], color="darkorange", ls="--", lw=1.2)
    axes.flat[-1].axis("off")
    fig.colorbar(image, ax=axes.ravel().tolist(), label="log10 Ne (m^-3)",
                 fraction=0.025, pad=0.02)
    fig.suptitle("Formal geocentric Provider electron-density scenes")
    fig.savefig(ASSETS / "02_ne_scene_comparison.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)

    baseline = fields["baseline"]
    peak = max(float(baseline.max()), 1.0)
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True, sharey=True)
    image = None
    for ax, scene in zip(axes.flat, ("tid", "spreadf", "es", "bubble")):
        anomaly = (fields[scene] - baseline) / peak
        limit = max(float(np.max(np.abs(anomaly))), 1.0e-6)
        image = ax.pcolormesh(
            cfg.BG_X, z_plot, anomaly.T, cmap="seismic", shading="auto",
            alpha=0.5, vmin=-limit, vmax=limit)
        ax.set_title("{} - baseline".format(scene))
        ax.set_xlabel("Ground range (km)")
        ax.set_ylabel("Physical height (km)")
        fig.colorbar(image, ax=ax, label="dNe / baseline peak")
        if scene == "es":
            h_es = float(cfg.ES["h_Es_km"])
            zoom = (z_plot >= h_es - 0.6) & (z_plot <= h_es + 0.6)
            ax.axhline(h_es, color="darkorange", ls="--", lw=1.2)
            inset = ax.inset_axes([0.42, 0.08, 0.46, 0.30])
            zoom_data = anomaly[:, zoom].T
            zoom_limit = max(float(np.max(np.abs(zoom_data))), 1.0e-8)
            zoom_image = inset.pcolormesh(
                cfg.BG_X, z_plot[zoom], zoom_data, cmap="seismic",
                shading="auto", alpha=0.5, vmin=-zoom_limit,
                vmax=zoom_limit)
            inset.axhline(h_es, color="black", ls="--", lw=0.7)
            inset.set_title(
                "Es display zoom: dNe / baseline peak; hEs +/- 0.6 km",
                fontsize=6.5)
            inset.set_xlabel("x (km)", fontsize=6)
            inset.set_ylabel("z (km)", fontsize=6)
            inset.tick_params(labelsize=6)
            colorbar_ax = inset.inset_axes([1.04, 0.0, 0.045, 1.0])
            colorbar = fig.colorbar(zoom_image, cax=colorbar_ax)
            colorbar.ax.tick_params(labelsize=5)
    fig.suptitle(
        "Physical electron-density perturbations "
        "(normalized by the common baseline peak)")
    fig.savefig(ASSETS / "03_ne_relative_perturbations.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)


def _bubble_figure(provider):
    x = np.linspace(350.0, 900.0, 300)
    z = np.linspace(60.0, 600.0, 300)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    params = cfg.BUBBLE_3D
    b_hat = geomag_unit_vector(cfg.GEOMAG)
    center = np.array([
        params["x0_km"], params["y0_km"], params["z0_km"]], dtype=float)
    points = np.stack((xx - center[0], np.zeros_like(xx), zz - center[2]),
                      axis=-1)
    parallel = points @ b_hat
    perpendicular2 = np.sum(points * points, axis=-1) - parallel * parallel
    old = params["delta_max"] * np.exp(
        -0.5 * perpendicular2 / params["sigma_perp_km"] ** 2
        -0.5 * parallel ** 2 / params["sigma_para_km"] ** 2)
    gate = 0.5 * (1.0 + np.tanh(
        (zz - params["bottomside_height_km"])
        / params["bottomside_transition_km"]))
    new = old * gate

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, data, title in (
            (axes[0], old, "Before: ungated depletion"),
            (axes[1], new, "After: F-region bottomside gate")):
        image = ax.pcolormesh(
            x, z, data.T, cmap="jet", shading="auto", vmin=0.0,
            vmax=params["delta_max"], alpha=0.5)
        ax.axhline(params["bottomside_height_km"], color="black",
                   ls="--", lw=1.0)
        ax.set_title(title)
        ax.set_xlabel("Ground range (km)")
        ax.set_ylabel("Physical height (km)")
        fig.colorbar(image, ax=ax, label="Relative depletion")
    ix = np.argmin(np.abs(x - params["x0_km"]))
    axes[2].plot(old[ix], z, color="0.5", label="before")
    axes[2].plot(new[ix], z, color="crimson", label="after")
    axes[2].axhline(params["bottomside_height_km"], color="black", ls="--")
    axes[2].set_xlabel("Relative depletion")
    axes[2].set_ylabel("Physical height (km)")
    axes[2].set_title("Centre-column morphology")
    axes[2].legend()
    fig.suptitle(
        "Parameterized Bubble: field-aligned Gaussian depletion tube; "
        "not a dynamical EPB plume")
    fig.savefig(ASSETS / "10_bubble_morphology_closure.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)


def _es_figure(baseline, es):
    x = cfg.BG_X
    h_es = float(cfg.ES["h_Es_km"])
    z_full = np.unique(np.r_[cfg.BG_Z, np.linspace(h_es - 0.6, h_es + 0.6, 121)])
    ne_base = baseline.sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, x, z_full)
    ne_es = es.sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, x, z_full)
    delta = ne_es - ne_base
    zoom = (z_full >= h_es - 0.6) & (z_full <= h_es + 0.6)
    midpoint = int(np.argmin(np.abs(x - cfg.RX_RANGE / 2.0)))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    image = axes[0].pcolormesh(
        x, z_full, np.log10(np.maximum(ne_es.T, 1.0)), cmap="jet",
        shading="auto", alpha=0.5)
    axes[0].axhline(h_es, color="darkorange", ls="--")
    axes[0].set_title("Full-height physical Ne")
    fig.colorbar(image, ax=axes[0], label="log10 Ne (m^-3)")
    image = axes[1].pcolormesh(
        x, z_full[zoom], delta[:, zoom].T / 1.0e10, cmap="jet",
        shading="auto", alpha=0.5)
    axes[1].axhline(h_es, color="black", ls="--")
    axes[1].set_title("Es physical dNe zoom (unscaled)")
    fig.colorbar(image, ax=axes[1], label="dNe (x10^10 m^-3)")
    axes[2].plot(delta[midpoint, zoom] / 1.0e10, z_full[zoom])
    axes[2].axhline(h_es, color="black", ls="--")
    axes[2].set_title("Mid-link vertical profile")
    axes[2].set_xlabel("dNe (x10^10 m^-3)")
    for ax in axes[:2]:
        ax.set_xlabel("Ground range (km)")
    for ax in axes:
        ax.set_ylabel("Physical height (km)")
    fig.suptitle(
        "Sporadic-E visualization: display zoom only; propagation Ne unscaled")
    fig.savefig(ASSETS / "11_es_physical_zoom.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    return ne_es, delta


def _ray_fan_outputs():
    summaries = {}
    for scene in SCENES:
        path = OUT / "ray_fan_unified_{}.json".format(scene)
        with path.open("r", encoding="utf-8") as stream:
            summaries[scene] = json.load(stream)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for ax, scene in zip(axes.flat, SCENES):
        image = plt.imread(OUT / "tx_fan_unified_{}.png".format(scene))
        ax.imshow(image)
        ax.set_title(scene)
        ax.axis("off")
    axes.flat[-1].axis("off")
    fig.suptitle(
        "Formal ECEF AH3D O-wave ray fans | 10 MHz | one-way | "
        "P2P/P-D/Ray-PE/scatter-power off")
    fig.savefig(ASSETS / "12_formal_complex_scene_ray_fans.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)

    statuses = sorted({
        status for summary in summaries.values()
        for status in summary["termination_counts"]
    })
    fig, ax = plt.subplots(figsize=(10, 5))
    base = np.zeros(len(SCENES))
    for status in statuses:
        values = np.array([
            summary["termination_counts"].get(status, 0)
            for summary in summaries.values()])
        ax.bar(SCENES, values, bottom=base, label=status)
        base += values
    ax.set_ylabel("Returned rays")
    ax.set_title("Formal ray-fan termination-state distribution")
    ax.legend(fontsize=8)
    fig.savefig(ASSETS / "13_formal_ray_fan_statuses.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    return summaries


def run() -> int:
    factory = ScenarioProviderFactory(cfg)
    providers, fields, z_plot = _sample_fields(factory)
    _scene_figures(fields, z_plot)
    _bubble_figure(providers["bubble"])
    h_es = float(cfg.ES["h_Es_km"])
    z_es_check = np.unique(np.r_[
        cfg.BG_Z, np.linspace(h_es - 0.6, h_es + 0.6, 121)])
    es_before_hash = providers["es"].get_metadata()["provider_hash"]
    es_before_values = providers["es"].sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, z_es_check)
    es_values, _ = _es_figure(providers["baseline"], providers["es"])
    es_after_hash = providers["es"].get_metadata()["provider_hash"]
    summaries = _ray_fan_outputs()

    z_check = np.array([0.0, 30.0, 59.9, 60.0, 110.0, 350.0])
    ne_base = providers["baseline"].sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, z_check)
    ne_bubble = providers["bubble"].sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X, z_check)
    depletion = np.divide(
        ne_base - ne_bubble, ne_base, out=np.zeros_like(ne_base),
        where=ne_base > 0.0)
    bubble_component = providers["bubble"].components[0]
    bubble_center = np.array([
        cfg.BUBBLE_3D["x0_km"], cfg.BUBBLE_3D["y0_km"],
        cfg.BUBBLE_3D["z0_km"]], dtype=float)
    b_hat = geomag_unit_vector(cfg.GEOMAG)
    perpendicular = np.cross(b_hat, np.array([0.0, 0.0, 1.0]))
    perpendicular /= np.linalg.norm(perpendicular)
    offset_km = 100.0
    along_fraction = max(
        bubble_component.depletion_fraction(
            *(bubble_center + sign * offset_km * b_hat), 0.0)
        for sign in (-1.0, 1.0))
    across_fraction = max(
        bubble_component.depletion_fraction(
            *(bubble_center + sign * offset_km * perpendicular), 0.0)
        for sign in (-1.0, 1.0))
    background = providers["baseline"].background
    strict = GeocentricLinkGridNeProvider3D(
        background.Ne_m3, background.ground_range_km, background.height_km,
        tx_lat_deg=background.tx_lat_deg, tx_lon_deg=background.tx_lon_deg,
        link_bearing_deg=background.link_bearing_deg,
        earth_radius_km=background.earth_radius_km,
        fill_value_m3=background.fill_value_m3,
        vertical_interpolation="pchip")
    above = np.array([60.0, 110.0, 250.0, 350.0, 600.0])
    extended_values = providers["baseline"].sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X[::20], above)
    strict_values = strict.sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, cfg.BG_X[::20], above)
    node_ix = np.flatnonzero(
        (background.ground_range_km > background.ground_range_km[0])
        & (background.ground_range_km < background.ground_range_km[-1]))[::20]
    node_x = background.ground_range_km[node_ix]
    node_z = background.height_km[::20]
    sampled_nodes = providers["baseline"].sample_link_plane(
        cfg.TX_GEOCENTRIC, cfg.LINK_BEARING_DEG, node_x, node_z)
    stored_nodes = background.Ne_m3[np.ix_(
        node_ix,
        np.arange(background.height_km.size)[::20])]

    checks = {
        "bubble_110km_depletion_le_1pct":
            float(depletion[:, 4].max()) <= 0.01,
        "bubble_350km_depletion_near_0p6":
            abs(float(depletion[:, 5].max()) - 0.6) <= 1.0e-3,
        "bubble_density_nonnegative":
            all(np.all(value >= 0.0) for value in fields.values()),
        "bubble_field_aligned_extent_exceeds_cross_field":
            along_fraction > across_fraction,
        "ground_tx_inside_provider":
            providers["baseline"].domain_margin_geocentric(
                *cfg.TX_GEOCENTRIC) >= 0.0,
        "ground_to_60km_strict_vacuum":
            float(ne_base[:, z_check < 60.0].max()) == 0.0,
        "above_60km_values_unchanged":
            np.array_equal(extended_values, strict_values),
        "pchip_preserves_interior_stored_grid_nodes":
            np.allclose(sampled_nodes, stored_nodes, rtol=1.0e-11, atol=1.0),
        "horizontal_outside_rejected":
            providers["baseline"].domain_margin_geocentric(
                cfg.TX_LAT, cfg.TX_LON, 100.0) >= 0.0
            and providers["baseline"].domain_margin_geocentric(
                cfg.TX_LAT - 40.0, cfg.TX_LON, 100.0) < 0.0,
        "es_visualization_provider_hash_unchanged":
            es_before_hash == es_after_hash,
        "es_visualization_provider_values_unchanged":
            np.array_equal(es_before_values, es_values),
        "es_visualization_values_finite":
            np.all(np.isfinite(es_values)),
        "current_provider_hashes_match_formal_outputs":
            all(providers[scene].get_metadata()["provider_hash"]
                == summaries[scene]["provider_hash"] for scene in SCENES),
        "five_formal_o_wave_fans_present":
            all((OUT / "tx_fan_unified_{}.png".format(scene)).stat().st_size
                > 10_000 for scene in SCENES),
        "all_fans_attempt_33":
            all(summary["ray_fan_diagnostics"]["O"]["attempted"] == 33
                for summary in summaries.values()),
        "all_fans_return_nonempty":
            all(summary["ray_fan_diagnostics"]["O"]["returned"] > 0
                for summary in summaries.values()),
        "all_fans_have_nonrejected_trajectory":
            all(any(status != "physics_domain_rejected" and count > 0
                    for status, count in summary["termination_counts"].items())
                for summary in summaries.values()),
        "all_attempts_accounted_for":
            all(
                summary["ray_fan_diagnostics"]["O"]["returned"]
                + sum(summary["ray_fan_diagnostics"]["O"][
                    "exception_counts"].values()) == 33
                for summary in summaries.values()),
        "no_tx_outside_provider_failure":
            all(not any("TX is outside" in key for key in
                        summary["ray_fan_diagnostics"]["O"][
                            "exception_counts"])
                for summary in summaries.values()),
        "only_o_wave":
            all(summary["wave_modes"] == ["O"]
                for summary in summaries.values()),
    }
    checks = {name: bool(passed) for name, passed in checks.items()}
    metrics = {
        "bubble_110km_max_relative_depletion":
            float(depletion[:, 4].max()),
        "bubble_350km_max_relative_depletion":
            float(depletion[:, 5].max()),
        "bubble_100km_offset_depletion": {
            "along_field": float(along_fraction),
            "cross_field": float(across_fraction),
        },
        "vacuum_max_Ne_m3": float(ne_base[:, :4].max()),
        "provider_hashes": {
            scene: summary["provider_hash"]
            for scene, summary in summaries.items()},
        "solver_configuration": {
            "integrator": "DOP853",
            "integration_profile": "smoke",
            "spatial_fd_base_km": {
                scene: min(
                    0.5,
                    min([
                        float(component["minimum_structure_scale_km"])
                        for component in summary["provider_metadata"].get(
                            "components", [])
                        if component.get("enabled", False)
                        and float(component.get(
                            "minimum_structure_scale_km", 0.0)) > 0.0
                    ], default=30.0) / 60.0)
                for scene, summary in summaries.items()
            },
            "spatial_fd_scales": [0.5, 1.0, 2.0],
            "stability_thresholds_unchanged": True,
        },
        "termination_counts": {
            scene: summary["termination_counts"]
            for scene, summary in summaries.items()},
    }
    report = {
        "schema_version": "complex_scene_ray_fan_closure_v1",
        "checks": checks,
        "metrics": metrics,
        "status": "passed" if all(checks.values()) else "failed",
    }
    report_path = ASSETS / "formal_complex_scene_ray_fan_acceptance.json"
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    artifact_paths = sorted(
        path for path in ASSETS.iterdir()
        if path.is_file() and path.name != "formal_complex_scene_manifest.json")
    manifest = {
        "schema_version": "complex_scene_ray_fan_manifest_v1",
        "artifacts": {
            path.name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in artifact_paths
        },
        "inputs": {
            str(path.relative_to(ROOT)).replace("\\", "/"): {
                "sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in [
                ROOT / "config.py",
                ROOT / "main2.py",
                ROOT / "models" / "ionosphere_components_3d.py",
                ROOT / "models" / "ne_provider_geocentric.py",
                ROOT / "models" / "scenario_provider_factory.py",
                ROOT / "models" / "unified_hybrid_pipeline.py",
                ROOT / "models" / "ray_tracer_ah_3d_ecef.py",
                *[OUT / "ray_fan_unified_{}.json".format(scene)
                  for scene in SCENES],
                *[OUT / "tx_fan_unified_{}.png".format(scene)
                  for scene in SCENES],
            ]
        },
    }
    (ASSETS / "formal_complex_scene_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    failures = [name for name, passed in checks.items() if not passed]
    print("complex_scene_ray_fan_closure: {} fail(s)".format(len(failures)))
    for name in failures:
        print(" FAIL {}".format(name))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(run())
