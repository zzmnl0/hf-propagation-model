"""Check that the explicit ray-wave route remains separate from AH entry points."""

from __future__ import annotations

import copy
import json
import sys
import types
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.constants import c

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config
import main
import models.unified_hybrid_pipeline
from models.geomagnetic_field_3d import ConstantGeomagneticField3D
from models.ne_provider_geocentric import GeocentricNeProvider3D
from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF, integrate_fixed_path_ah3d_ecef


OUT = Path(__file__).resolve().parents[1] / "output2" / "ray_wave_oneway" / "s4" / "isolation"
RADIUS = 6371.0
COEFFICIENT = 80.61638587963631
FREQUENCY = 10e6
TIME = 37.0


class UniformProvider(GeocentricNeProvider3D):
    earth_radius_km = RADIUS

    def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
        return .2 * FREQUENCY**2 / COEFFICIENT

    def get_metadata(self):
        return {**super().get_metadata(), "earth_radius_km": self.earth_radius_km}


def config_values():
    return copy.deepcopy({name: value for name, value in vars(config).items()
                          if not name.startswith("_")
                          and not isinstance(value, (types.ModuleType, type))
                          and not callable(value)})


def model():
    return RefractiveIndexAH3DECEF(
        UniformProvider(), earth_radius_km=RADIUS, freq_MHz=FREQUENCY / 1e6,
        wave_mode="O", geomagnetic_provider=ConstantGeomagneticField3D(
            45e-6 * np.ones(3) / np.sqrt(3)), plasma_coefficient=COEFFICIENT,
    )


def original_ah_trace():
    return integrate_fixed_path_ah3d_ecef(
        model(), tx_lat_deg=30., tx_lon_deg=120., tx_height_km=250.,
        elevation_deg=60., bearing_deg=40., path_length_km=.1,
        medium_time_s=TIME, profile="L2", stop_on_physical_events=False,
    )


def main_test():
    metrics = {}
    try:
        loaded_before = sorted(name for name in sys.modules
                               if name == "ray_wave_oneway" or name.startswith("ray_wave_oneway."))
        metrics["original_entry_imported_new_route"] = loaded_before
        if loaded_before:
            raise AssertionError("original AH entry imported the new route")
        before_config = config_values()
        before = original_ah_trace()
        from ray_wave_oneway.local_wave import MediumSnapshot, solve_slab
        from ray_wave_oneway.periodic import solve_periodic_slab
        from ray_wave_oneway.propagation import prepare_windows
        from models.local_wave import MediumSnapshot as CompatibilitySnapshot
        from models.local_wave_periodic import solve_periodic_slab as compatibility_periodic
        from models.oneway_hybrid import prepare_windows as compatibility_prepare

        current = model()
        position = current.position_to_ecef(30., 120., 250.)
        medium = MediumSnapshot.from_ah_model(current, position, TIME)
        layer = replace(medium, Ne_m3=1.1 * medium.Ne_m3)
        wave = solve_slab(medium, layer, medium, medium, np.eye(3), 0., 0., 2 * c / FREQUENCY)
        if not np.all(np.isfinite(wave.S)):
            raise AssertionError("new route did not produce a finite scattering matrix")
        if (CompatibilitySnapshot is not MediumSnapshot
                or compatibility_periodic is not solve_periodic_slab
                or compatibility_prepare is not prepare_windows):
            raise AssertionError("compatibility exports do not resolve to the isolated implementation")
        after = original_ah_trace()
        after_config = config_values()
        np.testing.assert_equal(before_config, after_config)
        metrics["config_equal"] = True
        for key in ("trajectory_ecef_km", "wavevector_ecef_km^-1", "trajectory_group_time_s"):
            np.testing.assert_array_equal(before[key], after[key])
            metrics[key + "_max_abs_difference"] = float(np.max(np.abs(before[key] - after[key])))
        if before["termination_status"] != after["termination_status"]:
            raise AssertionError("original AH trace termination changed")
        metrics["termination_status"] = after["termination_status"]
        metrics["compatibility_exports_same_objects"] = True
        metrics["new_route_finite_S"] = True

        # Mutate the original AH settings only in this test process.  A new
        # route model and an actual AH plane crossing must remain identical.
        from models.ray_tracer_ah_3d_ecef import _integration_profile
        from ray_wave_oneway.cases import _context
        from ray_wave_oneway.config import ControlledCase, ONEWAY_AH_FD_SETTINGS
        from ray_wave_oneway.propagation import _trace_to_plane

        def isolated_snapshot():
            case = ControlledCase(grid_size=1)
            models, medium, window = _context(case, case.frequency_hz, 0.0)
            port = next(p for p in window["slabs"][0].left_ports.forward if p.label == "O")
            target = window["entrance_ecef_km"]
            ray = _trace_to_plane(models["O"], port, target - np.array([0.1, 0, 0]),
                                  target, window["basis"], medium.time_s, window["left_medium"])
            return models["O"], ray

        route_before_model, route_before_ray = isolated_snapshot()
        update2_before = copy.deepcopy(config.UPDATE2)
        try:
            config.UPDATE2["AH3D_SPATIAL_FD_BASE_KM"] = 0.2
            config.UPDATE2["AH3D_INTEGRATION_PROFILE"] = "smoke"
            config.UPDATE2["AH3D_INTEGRATION"]["L2"]["rtol"] = 3.0e-5
            legacy = model()
            assert legacy.spatial_fd_base_km == 0.2 and legacy.integration_profile == "smoke"
            assert _integration_profile("L2")["rtol"] == 3.0e-5
            route_after_model, route_after_ray = isolated_snapshot()
            for key, value in ONEWAY_AH_FD_SETTINGS.items():
                assert getattr(route_before_model, key) == value == getattr(route_after_model, key)
            assert route_before_model.integration_profile == route_after_model.integration_profile == "L2"
            for key in ("phase_rad", "ray_phase_rad", "group_time_s", "endpoint_ecef_km",
                        "endpoint_k_ecef_rad_km", "path_length_km", "terminal_plane_error_km"):
                np.testing.assert_array_equal(route_before_ray[key], route_after_ray[key])
            metrics["legacy_update2_mutation_changes_legacy_settings"] = True
            metrics["new_context_and_actual_ray_unchanged"] = True
        finally:
            config.UPDATE2.clear()
            config.UPDATE2.update(update2_before)
        np.testing.assert_equal(before_config, config_values())
        status, detail = "PASS", None
    except Exception as exc:
        status, detail = "FAIL", str(exc)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "isolation.json").write_text(json.dumps({"status": status, "detail": detail,
                                                     "metrics": metrics}, indent=2), encoding="utf-8")
    print(f"Isolation {status}" + ("" if detail is None else f": {detail}"))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main_test())
