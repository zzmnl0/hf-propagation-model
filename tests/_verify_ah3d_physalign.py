"""Bounded S4 coefficient routing and independent B=0 physics checks."""
import argparse
import ast
import hashlib
import json
import math
import os
from pathlib import Path
import sys

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import config
from models.absorption_collision import (
    complex_absorption_ecef_ray_dB, complex_refractive_index_ah,
)
from models.dispersion_ah import appleton_hartree_roots, C_KM_S
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
    raise AssertionError("invalid coefficient accepted")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    output = Path(parser.parse_args().output)
    output.mkdir(parents=True, exist_ok=True)
    result_path = output / "engineering_result.json"
    if result_path.exists():
        raise FileExistsError(result_path)
    # CODATA 2022 literals: scipy/constants/_codata.py txt2022, frozen locally.
    derived = 1.602176634e-19**2 / (4 * math.pi**2 * 8.8541878188e-12 * 9.1093837139e-31)
    si = config.AH3D_PLASMA_COEFFICIENT_SI
    assert math.isclose(si, derived, rel_tol=1e-15)
    assert config.K_FP == 8.980
    legacy, ph = config.K_FP**2, 80.61638494312905
    provider = GeocentricGridNeProvider3D(
        np.full((2, 2, 3), 1e11), [20., 40.], [110., 130.],
        [80., 130., 180.], 6371., interpolation="radial_hermite")
    base = dict(earth_radius_km=6371., freq_MHz=10., gyro_MHz=0.)
    compatibility = config.UPDATE2["PHARLAP_RADIAL_B0"]
    saved = dict(compatibility)
    old_gyro = config.GEOMAG["fH_MHz"]
    rows = []
    try:
        config.GEOMAG["fH_MHz"] = 0.
        for compat in (False, True):
            compatibility["enable"] = compat
            for value in (None, legacy, si, ph):
                expected = (ph if compat else si) if value is None else value
                kwargs = dict(base, plasma_coefficient=value,
                              compatibility_mode="pharlap_radial_b0" if compat else None)
                direct = RefractiveIndexAH3DECEF(provider, **kwargs)
                request = UnifiedPropagationRequest(
                    scenario="S4_coefficient_contract", provider=provider,
                    coordinate_model=COORDINATE_MODEL, earth_radius_km=6371., freq_MHz=10.,
                    plasma_coefficient=value)
                factory = _ah_model(request, "O", "L2")
                assert direct.plasma_coefficient == factory.plasma_coefficient == expected
                assert factory.compatibility_mode == direct.compatibility_mode
                assert direct.spatial_fd_base_km == (
                    saved["spatial_fd_base_km"] if compat else config.AH3D_L2_SPATIAL_FD_BASE_KM)
                rows.append(dict(compatibility=compat, requested=value, actual=expected))
        compatibility["enable"] = False
        for bad in (0., -1., float("nan"), float("inf")):
            rejected(lambda: RefractiveIndexAH3DECEF(provider, **base, plasma_coefficient=bad))
            rejected(lambda: complex_refractive_index_ah(
                10., 1e11, 1e5, [1., 0., 0.], plasma_coefficient=bad))

        launch = dict(tx_lat_deg=30., tx_lon_deg=120., tx_height_km=100.,
                      elevation_deg=35., bearing_deg=90.)
        state_errors = []
        for coefficient in (legacy, si, ph):
            model = RefractiveIndexAH3DECEF(provider, **base, plasma_coefficient=coefficient)
            state = launch_state_ah3d_ecef(model, **launch)
            expected_n2 = 1 - coefficient * 1e11 / 1e14
            assert abs(model.n2_raw(state[:3], state[3:6]) - expected_n2) <= 2e-16
            assert abs(np.dot(state[3:6], state[3:6]) / model.k0**2 - expected_n2) <= 3e-16
            residual = hamiltonian_residual_ah3d_ecef(state, model)
            assert residual <= 1e-12
            state_errors.append(float(residual))
            for nu in (0., 1e5):
                expected = np.sqrt(complex(1 - coefficient * 1e11 / 1e14 /
                                          complex(1., -nu / (2 * math.pi * 1e7))))
                expected = complex(expected.real, abs(expected.imag))
                actual = complex_refractive_index_ah(
                    10., 1e11, nu, [1., 0., 0.], geomag={"fH_MHz": 0.},
                    plasma_coefficient=coefficient)
                assert abs(actual - expected) <= 2e-16
            vacuum = complex_refractive_index_ah(
                10., 0., 1e5, [1., 0., 0.], geomag={"fH_MHz": 0.},
                plasma_coefficient=coefficient)
            assert vacuum == 1.

        omega = 2 * math.pi * 1e7
        old_root = appleton_hartree_roots(omega, 1e11, [1., 0., 0.], [0., 0., 1.], 0.)["O"]
        assert old_root == 1 - legacy * 1e11 / 1e14
        assert complex_refractive_index_ah(10., 1e11, 1e5, [1., 0., 0.], geomag={"fH_MHz": 0.}) == complex_refractive_index_ah(
            10., 1e11, 1e5, [1., 0., 0.], geomag={"fH_MHz": 0.}, plasma_coefficient=legacy)

        model = RefractiveIndexAH3DECEF(provider, **base)
        explicit = RefractiveIndexAH3DECEF(provider, **base, plasma_coefficient=si)
        ray = integrate_fixed_path_ah3d_ecef(model, **launch, path_length_km=.05)
        other = integrate_fixed_path_ah3d_ecef(explicit, **launch, path_length_km=.05)
        for field in ("trajectory_ecef_km", "wavevector_ecef_km^-1", "trajectory_group_path_km"):
            np.testing.assert_array_equal(ray[field], other[field])
        segment = ray["trajectory_ecef_km"][-1] - ray["trajectory_ecef_km"][0]
        path_ray = dict(ray, trajectory_ecef_km=np.array([
            ray["trajectory_ecef_km"][0], ray["trajectory_ecef_km"][-1]]))
        absorption_rows = []
        for coefficient in (legacy, si, ph):
            path_ray["plasma_coefficient"] = coefficient
            params = {"collision_model": "constant", "nu0_s": 1e5}
            absorption = complex_absorption_ecef_ray_dB(path_ray, provider, params)
            doubled = complex_absorption_ecef_ray_dB(path_ray, provider, params, two_way=True)
            index = complex_refractive_index_ah(
                10., 1e11, 1e5, segment, geomag={"fH_MHz": 0.}, plasma_coefficient=coefficient)
            expected = 8.686 * 2 * (omega / C_KM_S / 1000) * index.imag * np.linalg.norm(segment) * 1000
            assert math.isclose(absorption["abs_total_dB"], expected, rel_tol=2e-14)
            assert absorption["abs_total_dB"] > 0.
            assert doubled["abs_total_dB"] == 2 * absorption["abs_total_dB"]
            assert absorption["plasma_coefficient"] == coefficient
            assert absorption["plasma_coefficient_source"] == "ray_metadata"
            zero = complex_absorption_ecef_ray_dB(path_ray, provider, dict(params, nu0_s=0.))
            assert zero["abs_total_dB"] == 0.
            absorption_rows.append(absorption)
        for bad in (0., -1., float("nan"), float("inf")):
            rejected(lambda: complex_absorption_ecef_ray_dB(dict(path_ray, plasma_coefficient=bad), provider))
        del path_ray["plasma_coefficient"]
        fallback = complex_absorption_ecef_ray_dB(path_ray, provider, params)
        assert fallback["plasma_coefficient"] == legacy
        assert fallback["plasma_coefficient_source"] == "legacy_K_FP_squared"
        assert fallback["abs_total_dB"] == absorption_rows[0]["abs_total_dB"]
        assert provider.get_Ne_geocentric(30., 120., 100.) == 1e11
        assert config.K_FP == 8.980 and saved["plasma_coefficient"] == ph
        p31 = ROOT / "tests/_verify_radial_pchip_physics.py"
        tree = ast.parse(p31.read_text(encoding="utf-8-sig"))
        operational = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in
                       {"physics_trace", "physics_case", "physics_register", "continuous_control", "physics_preconditions"}]
        constructors = [node for function in operational for node in ast.walk(function) if isinstance(node, ast.Call)
                        and isinstance(node.func, ast.Name) and node.func.id == "RefractiveIndexAH3DECEF"]
        assert constructors and all(any(k.arg == "plasma_coefficient" for k in node.keywords) for node in constructors)
        codata_path = Path(np.__file__).parents[1] / "scipy/constants/_codata.py"
        codata = codata_path.read_text(encoding="utf-8")
        assert "txt2022" in codata and "9.109 383 7139" in codata and "8.854 187 8188" in codata
        result = dict(task="AH3D-PH474-S4-PHYSALIGN-20261029", decision="PASS",
                      coefficient_rows=rows, independent_CODATA2022=derived,
                      state_residuals=state_errors, absorption=absorption_rows,
                      legacy_fallback=fallback, short_path_km=.05,
                      short_path_role="routing regression only", p31_static_explicit_coefficient="PASS",
                      codata_source=str(codata_path), codata_sha256=hashlib.sha256(codata_path.read_bytes()).hexdigest(),
                      codata_formula_minus_frozen=derived-si,
                      natural_trajectories=0, full_scientific_acceptance=False)
        result_path.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    finally:
        compatibility.clear()
        compatibility.update(saved)
        config.GEOMAG["fH_MHz"] = old_gyro
    print("PASS: S4 coefficient routing, B=0 physics, ECEF absorption, legacy fallback")


if __name__ == "__main__":
    main()
