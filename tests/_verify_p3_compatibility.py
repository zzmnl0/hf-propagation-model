"""P3 production-contract checks; no PH agreement claim is made here."""
from pathlib import Path
import sys
import platform

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _verify_production_boundary import install_guard, check_loaded
platform.uname()
install_guard()
import numpy as np
import config
from models.dispersion_ah import appleton_hartree_roots, dispersion_D, group_velocity, C_KM_S
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    COORDINATE_MODEL, RefractiveIndexAH3DECEF, integrate_fixed_path_ah3d_ecef,
    launch_state_ah3d_ecef, hamiltonian_residual_ah3d_ecef,
)
from models.unified_hybrid_pipeline import UnifiedPropagationRequest, _ah_model


def rejected(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError("unsupported input was accepted")


def provider(profile=(1e11, 1e11, 1e11), **kwargs):
    return GeocentricGridNeProvider3D(
        np.broadcast_to(profile, (2, 2, 3)).copy(), [20, 40], [110, 130],
        [80, 130, 180], 6371, interpolation="radial_hermite", **kwargs)


def main():
    p = provider()
    mode = dict(earth_radius_km=6371, freq_MHz=10,
                compatibility_mode="pharlap_radial_b0")
    m = RefractiveIndexAH3DECEF(p, **mode)
    assert m.plasma_coefficient == 80.61638494312905
    assert m.spatial_fd_base_km == 0.125
    launch = dict(tx_lat_deg=30, tx_lon_deg=120, tx_height_km=100,
                  elevation_deg=22, bearing_deg=0)
    state = launch_state_ah3d_ecef(m, **launch)
    assert abs(hamiltonian_residual_ah3d_ecef(state, m)) <= 1e-10
    root = lambda k, w: dispersion_D(k, w, appleton_hartree_roots(
        w, 1e11, k, [0, 0, 1], 0, plasma_coefficient=m.plasma_coefficient)["O"])
    n = np.sqrt(1 - m.plasma_coefficient * 1e11 / 1e14)
    velocity = group_velocity(root, state[3:6], m.omega)
    assert np.isclose(np.linalg.norm(velocity), C_KM_S * n, rtol=2e-12)
    ray = integrate_fixed_path_ah3d_ecef(m, **launch, path_length_km=1)
    assert np.isclose(ray["group_path_km"], 1 / n, rtol=2e-12)
    assert ray["trajectory_group_path_km"][-1] == ray["group_path_km"]
    np.testing.assert_array_equal(ray["trajectory_group_path_km"],
                                  C_KM_S * ray["trajectory_group_time_s"])
    assert ray["trajectory_group_tangent_ecef_semantics"] == "normalized_position_chord_estimate"
    rejected(lambda: RefractiveIndexAH3DECEF(p, **dict(mode, gyro_MHz=0.1)))
    rejected(lambda: RefractiveIndexAH3DECEF(p, **dict(mode, collision_s=1)))
    rejected(lambda: RefractiveIndexAH3DECEF(p, **dict(mode, earth_radius_km=6370)))
    rejected(lambda: RefractiveIndexAH3DECEF(p, **dict(mode, plasma_coefficient=np.nan)))
    rejected(lambda: p.get_Ne_geocentric(30, 120, 79))
    rejected(lambda: p.get_Ne_geocentric(30, 120, 181))
    rejected(lambda: p.get_Ne_geocentric(19, 120, 100))
    rejected(lambda: p.get_Ne_geocentric(30, 109, 100))
    rejected(lambda: provider((0, 0, 1e8)))
    from scipy.io import loadmat
    from scipy.interpolate import CubicHermiteSpline
    original_grid = loadmat(Path(__file__).resolve().parents[1] /
        "output2/a4_pharlap/c1/c1m/a2_20260928_v1/a2/grids/smooth_s5.mat",
        simplify_cells=True)
    height, density = original_grid["height_km"], original_grid["ne_profile_m3"]
    spline = CubicHermiteSpline(height, density, np.gradient(density, height, edge_order=2))
    assert np.isclose(spline(19 + 1/3), -228.3014205749986, rtol=1e-12)
    rejected(lambda: GeocentricGridNeProvider3D(
        original_grid["Ne_m3"], original_grid["latitude_deg"], original_grid["longitude_deg"],
        height, 6371, interpolation="radial_hermite"))
    rejected(lambda: provider(fill_value_m3=1))
    bad = np.full((2, 2, 3), 1e11)
    bad[1, 1, 1] += 1
    rejected(lambda: GeocentricGridNeProvider3D(
        bad, [20, 40], [110, 130], [80, 130, 180], 6371,
        interpolation="radial_hermite"))
    zero = provider((0, 0, 0))
    assert zero.get_Ne_geocentric(30, 120, 79) == 0
    assert zero.get_Ne_geocentric(30, 120, 181) == 0
    ordinary = GeocentricGridNeProvider3D(
        np.full((2, 2, 3), 1e11), [20, 40], [110, 130], [80, 130, 180], 6371)
    rejected(lambda: RefractiveIndexAH3DECEF(ordinary, **mode))
    original = RefractiveIndexAH3DECEF(ordinary, earth_radius_km=6371, freq_MHz=10)
    explicit = RefractiveIndexAH3DECEF(ordinary, earth_radius_km=6371, freq_MHz=10,
                                     plasma_coefficient=config.K_FP**2,
                                     spatial_fd_base_km=0.5)
    a = integrate_fixed_path_ah3d_ecef(original, **launch, path_length_km=1)
    b = integrate_fixed_path_ah3d_ecef(explicit, **launch, path_length_km=1)
    for key in ("trajectory_ecef_km", "wavevector_ecef_km^-1", "trajectory_group_path_km"):
        np.testing.assert_array_equal(a[key], b[key])
    compatibility = config.UPDATE2["PHARLAP_RADIAL_B0"]
    enabled, gyro = compatibility["enable"], config.GEOMAG["fH_MHz"]
    try:
        compatibility["enable"] = True
        config.GEOMAG["fH_MHz"] = 0
        request = UnifiedPropagationRequest(scenario="p3_contract", provider=p,
                                           coordinate_model=COORDINATE_MODEL,
                                           earth_radius_km=6371, freq_MHz=10)
        pipeline_model = _ah_model(request, "O", "reference")
        assert pipeline_model.plasma_coefficient == m.plasma_coefficient
        assert pipeline_model.spatial_fd_base_km == m.spatial_fd_base_km
        np.testing.assert_array_equal(launch_state_ah3d_ecef(pipeline_model, **launch), state)
        pipeline_ray = integrate_fixed_path_ah3d_ecef(pipeline_model, **launch, path_length_km=1)
        for key in ("trajectory_ecef_km", "wavevector_ecef_km^-1", "trajectory_group_path_km"):
            np.testing.assert_array_equal(pipeline_ray[key], ray[key])
        for key in ("plasma_coefficient", "compatibility_mode", "compatibility_scope",
                    "coordinate_model", "spatial_fd_base_km"):
            assert pipeline_model.metadata()[key] == m.metadata()[key]
        request.provider = ordinary
        rejected(lambda: _ah_model(request, "O", "reference"))
    finally:
        compatibility["enable"], config.GEOMAG["fH_MHz"] = enabled, gyro
    assert check_loaded()
    print("P3 production compatibility: all contract checks passed; 0 fail(s)")


if __name__ == "__main__":
    main()
