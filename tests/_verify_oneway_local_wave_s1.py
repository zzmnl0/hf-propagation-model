"""S1 bounded acceptance for the local one-way Maxwell interface."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import psutil
from scipy.constants import c

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from models.dispersion_ah import appleton_hartree_roots, select_root
from ray_wave_oneway.local_wave import (
    MediumSnapshot,
    modal_ports,
    position_ecef_km,
    position_local_m,
    solve_slab,
    wavevector_ecef_rad_km,
    wavevector_local_rad_m,
)
from models.ne_provider_geocentric import GeocentricNeProvider3D, WGS84NeProviderAdapter
from models.ray_tracer_ah_3d_ecef import (
    COORDINATE_MODEL, WGS84_COORDINATE_MODEL, RefractiveIndexAH3DECEF,
)
from models.vector_pe_3d import cold_plasma_dielectric_tensor
from config import K_FP
from utils import ecef_to_geocentric, ecef_to_geodetic


OUT = Path(__file__).resolve().parents[1] / "output2" / "ray_wave_oneway" / "s4" / "isolation" / "s1"
FREQUENCY = 10e6
COEFFICIENT = 80.61638587963631
K0 = 2 * np.pi * FREQUENCY / c
THICKNESS = 2 * c / FREQUENCY
TIME = 37.0
POSITION = np.array([6371., 10., -20.])
BASIS = np.eye(3)
EARTH_RADIUS_KM = 6371.0
MAX_RSS = 256 * 1024 ** 2
MIN_AVAILABLE = 1024 ** 3
MAX_SECONDS = 300.0
START = time.monotonic()
PROCESS = psutil.Process()
RESOURCES = []


def resource_check():
    available = psutil.virtual_memory().available
    memory = PROCESS.memory_info()
    rss = memory.rss
    peak_rss = getattr(memory, "peak_wset", rss)
    elapsed = time.monotonic() - START
    RESOURCES.append({"available_bytes": available, "rss_bytes": rss,
                      "peak_rss_bytes": peak_rss, "elapsed_s": elapsed})
    if available < MIN_AVAILABLE or peak_rss > MAX_RSS or elapsed > MAX_SECONDS:
        raise RuntimeError(
            f"resource limit: available={available}, peak_rss={peak_rss}, elapsed={elapsed:.2f}"
        )
    return RESOURCES[-1]


def medium(x=0.0, field=None, collision=0.0, position=POSITION):
    return MediumSnapshot(
        FREQUENCY, TIME, position, x * FREQUENCY**2 / COEFFICIENT,
        np.zeros(3) if field is None else field, collision, COEFFICIENT,
        EARTH_RADIUS_KM, COORDINATE_MODEL, 0.02,
    )


def compare(metrics, name, actual, limit):
    value = float(actual)
    metrics[name] = {"value": value, "limit": limit, "pass": bool(np.isfinite(value) and value <= limit)}


def scalar_fresnel(result, eps_left, eps_slab, eps_right, ku, kv):
    """Independent TE/TM Fresnel-Fabry-Perot amplitudes at the port planes."""
    k_t_squared = ku * ku + kv * kv
    q = [np.sqrt(complex(K0 * K0 * ep - k_t_squared)) for ep in (eps_left, eps_slab, eps_right)]
    q = [value if value.imag >= 0 else -value for value in q]
    expected = np.zeros((4, 4), dtype=complex)
    for index, label in enumerate(("TE", "TM")):
        if label == "TE":
            admittance = q
            component = 1 if k_t_squared == 0 else 2
        else:
            admittance = [ep / wave for ep, wave in zip((eps_left, eps_slab, eps_right), q)]
            component = 2 if k_t_squared == 0 else 1
        y1, y2, y3 = admittance
        r12, r23 = (y1 - y2) / (y1 + y2), (y2 - y3) / (y2 + y3)
        r32, r21 = -r23, -r12
        t12, t23 = 2 * y1 / (y1 + y2), 2 * y2 / (y2 + y3)
        t32, t21 = 2 * y3 / (y3 + y2), 2 * y2 / (y2 + y1)
        phase = np.exp(1j * q[1] * THICKNESS)
        denominator = 1 + r12 * r23 * phase**2
        rleft = (r12 + r23 * phase**2) / denominator
        rright = (r32 + r21 * phase**2) / denominator
        tleft = t12 * t23 * phase / denominator
        tright = t32 * t21 * phase / denominator
        lf = result.left_ports.forward[index].electric_suv[component]
        lb = result.left_ports.backward[index].electric_suv[component]
        rf = result.right_ports.forward[index].electric_suv[component]
        rb = result.right_ports.backward[index].electric_suv[component]
        expected[index, index] = rleft * lf / lb
        expected[index + 2, index] = tleft * lf / rf
        expected[index, index + 2] = tright * rb / lb
        expected[index + 2, index + 2] = rright * rb / rf
    return expected


def matrix_properties(result, metrics, name):
    gram = result.S.conj().T @ result.S
    compare(metrics, name + "_unitary", np.linalg.norm(gram - np.eye(4), ord=2), 1e-5)
    compare(metrics, name + "_ports_1w", max(
        abs(abs(port.power_W) - 1) for group in
        (result.left_ports.forward, result.left_ports.backward,
         result.right_ports.forward, result.right_ports.backward) for port in group
    ), 1e-10)


def write_result(metrics, arrays, status, detail=None):
    OUT.mkdir(parents=True, exist_ok=True)
    data = {
        "status": status,
        "detail": detail,
        "frequency_hz": FREQUENCY,
        "medium_time_s": TIME,
        "plasma_coefficient": COEFFICIENT,
        "reference_area_m2": 1.0,
        "thickness_m": THICKNESS,
        "metrics": metrics,
        "resource_at_end": {
            "available_bytes": psutil.virtual_memory().available,
            "rss_bytes": PROCESS.memory_info().rss,
            "elapsed_s": time.monotonic() - START,
        },
        "resource_samples": RESOURCES,
    }
    (OUT / "metrics.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    if arrays:
        np.savez_compressed(OUT / "representative_fields.npz", **arrays)


def main():
    metrics = {}
    arrays = {}
    try:
        resource_check()
        vacuum = medium()
        base = medium(0.2)
        layer = medium(0.3)
        field = 45e-6 * np.ones(3) / np.sqrt(3)

        # C01 and C02 check transparent phase and both incident directions.
        for name, material, ku, kv in (
            ("C01", vacuum, 0., 0.),
            ("C02", base, .2 * K0, .1 * K0),
        ):
            resource_check()
            result = solve_slab(material, material, material, material, BASIS, ku, kv, THICKNESS)
            q = np.sqrt(K0**2 * (1 - COEFFICIENT * material.Ne_m3 / FREQUENCY**2) - ku**2 - kv**2)
            phase = np.exp(1j * q * THICKNESS)
            expected = np.block([[np.zeros((2, 2)), phase * np.eye(2)],
                                 [phase * np.eye(2), np.zeros((2, 2))]])
            compare(metrics, name + "_uniform_phase", np.max(np.abs(result.S - expected)), 1e-6)
            compare(metrics, name + "_background", np.max(np.abs(result.S - result.S0)), 1e-8)
            matrix_properties(result, metrics, name)

        # C03-C05 compare all four two-polarization blocks with Fresnel theory.
        for name, right, slab, background, ku in (
            ("C03", vacuum, layer, vacuum, 0.),
            ("C04", base, layer, base, .3 * K0),
            ("C05", base, medium(.3, collision=1e5), base, .3 * K0),
        ):
            resource_check()
            result = solve_slab(vacuum, slab, right, background, BASIS, ku, 0, THICKNESS)
            eps_left = vacuum.epsilon()[0, 0]
            eps_layer = slab.epsilon()[0, 0]
            eps_right = right.epsilon()[0, 0]
            analytic = scalar_fresnel(result, eps_left, eps_layer, eps_right, ku, 0)
            compare(metrics, name + "_fresnel_complex", np.max(np.abs(result.S - analytic)), 1e-6)
            arrays[name + "_S"] = result.S
            arrays[name + "_S0"] = result.S0
            if name != "C05":
                matrix_properties(result, metrics, name)
            else:
                gram = result.S.conj().T @ result.S
                compare(metrics, name + "_passivity", max(0., np.linalg.eigvalsh(gram - np.eye(4)).max()), 1e-5)
                for index, incident in enumerate((np.array([1, 0, 0, 0], complex),
                                                   np.array([1, 1j, .5, -.3j], complex))):
                    incident = incident / np.linalg.norm(incident)
                    outgoing = result.S @ incident
                    absorption = result.absorbed_power_W(incident, quadrature_order=64)
                    coarse = result.absorbed_power_W(incident, quadrature_order=32)
                    compare(metrics, f"{name}_closure_{index}", abs(1 - np.vdot(outgoing, outgoing).real - absorption), 1e-5)
                    compare(metrics, f"{name}_quadrature_{index}", abs(absorption - coarse), 1e-5)
                    arrays[f"{name}_incident_{index}"] = incident
                    arrays[f"{name}_electric_mid_{index}"] = result.field_suv(THICKNESS / 2, incident)[0]
            if name == "C04":
                compare(metrics, name + "_background_interface_nonzero", max(0., 1e-4 - np.max(np.abs(result.S0[:2, :2]))), 0.)

        resource_check()
        # C06 compares circular-polarization eigenvalues and handedness analytically.
        for sign in (1, -1):
            magnetized = medium(.2, np.array([sign * 45e-6, 0., 0.]))
            result = solve_slab(magnetized, magnetized, magnetized, magnetized, BASIS, 0, 0, THICKNESS)
            gyro_ratio = 27_992_489_872.333 * 45e-6 / FREQUENCY
            expected = {1: 1 - .2 / (1 - gyro_ratio), -1: 1 - .2 / (1 + gyro_ratio)}
            for port in result.left_ports.forward:
                plus = abs(port.electric_suv[1] + 1j * sign * port.electric_suv[2])
                minus = abs(port.electric_suv[1] - 1j * sign * port.electric_suv[2])
                handedness = 1 if minus > plus else -1
                n2 = np.dot(port.k_suv_rad_m, port.k_suv_rad_m) / K0**2
                compare(metrics, f"C06_{sign}_{port.label}_circular_root", abs(n2 - expected[handedness]), 1e-8)
                compare(metrics, f"C06_{sign}_{port.label}_circular_polarization", min(plus, minus) / max(plus, minus), 1e-10)
            matrix_properties(result, metrics, f"C06_{sign}")

        resource_check()
        # C07 uses a time-sensitive provider in the actual snapshot path.
        class Provider(GeocentricNeProvider3D):
            def __init__(self):
                self.calls = []
                self.earth_radius_km = EARTH_RADIUS_KM

            def get_Ne_geocentric(self, lat_deg, lon_deg, height_km, t=0.0):
                self.calls.append(("Ne", float(t), np.array([lat_deg, lon_deg, height_km])))
                return (.2 if t == TIME else .4) * FREQUENCY**2 / COEFFICIENT

            def get_B_ecef_T(self, position, t):
                self.calls.append(("B", float(t), np.asarray(position).copy()))
                return field if t == TIME else -2 * field

            def get_metadata(self):
                return {**super().get_metadata(), "earth_radius_km": self.earth_radius_km}

        provider = Provider()
        ah_model = RefractiveIndexAH3DECEF(
            provider, earth_radius_km=EARTH_RADIUS_KM, freq_MHz=FREQUENCY / 1e6,
            geomagnetic_provider=provider, plasma_coefficient=COEFFICIENT,
        )
        magnetized = MediumSnapshot.from_ah_model(ah_model, POSITION, TIME)
        spherical_calls = list(provider.calls)
        compare(metrics, "C07_provider_time", max(abs(call[1] - TIME) for call in spherical_calls), 0.)
        ne_call = next(call for call in spherical_calls if call[0] == "Ne")
        b_call = next(call for call in spherical_calls if call[0] == "B")
        compare(metrics, "C07_provider_spherical_coordinates", np.linalg.norm(
            ne_call[2] - ecef_to_geocentric(POSITION, EARTH_RADIUS_KM)), 1e-10)
        compare(metrics, "C07_provider_position_km", np.linalg.norm(b_call[2] - POSITION), 0.)
        compare(metrics, "C07_provider_earth_radius", abs(magnetized.earth_radius_km - ah_model.earth_radius_km), 0.)
        compare(metrics, "C07_provider_model_frequency", abs(magnetized.frequency_hz - ah_model.freq_MHz * 1e6), 0.)
        compare(metrics, "C07_provider_model_coefficient", abs(magnetized.plasma_coefficient - ah_model.plasma_coefficient), 0.)
        compare(metrics, "C07_provider_model_cutoff", abs(magnetized.near_cutoff_n2 - ah_model.near_cutoff_n2), 0.)
        provider.calls.clear()
        wgs_source = WGS84NeProviderAdapter(provider)
        wgs_model = RefractiveIndexAH3DECEF(
            wgs_source, earth_radius_km=EARTH_RADIUS_KM, freq_MHz=FREQUENCY / 1e6,
            geomagnetic_provider=provider, plasma_coefficient=COEFFICIENT,
            coordinate_model=WGS84_COORDINATE_MODEL,
        )
        wgs_snapshot = MediumSnapshot.from_ah_model(wgs_model, POSITION, TIME)
        wgs_calls = list(provider.calls)
        wgs_ne_call = next(call for call in wgs_calls if call[0] == "Ne")
        compare(metrics, "C07_WGS84_time", max(abs(call[1] - TIME) for call in wgs_calls), 0.)
        compare(metrics, "C07_WGS84_geodetic_coordinates", np.linalg.norm(
            wgs_ne_call[2] - ecef_to_geodetic(POSITION)), 1e-10)
        compare(metrics, "C07_WGS84_density", abs(wgs_snapshot.Ne_m3 - magnetized.Ne_m3), 0.)
        compare(metrics, "C07_WGS84_tensor", np.max(np.abs(wgs_snapshot.epsilon() - magnetized.epsilon())), 1e-12)
        compare(metrics, "C07_WGS84_model", float(wgs_snapshot.coordinate_model != wgs_model.coordinate_model), 0.)
        try:
            MediumSnapshot.from_providers(
                FREQUENCY, TIME, POSITION, wgs_source, provider,
                plasma_coefficient=COEFFICIENT, earth_radius_km=EARTH_RADIUS_KM,
                coordinate_model=COORDINATE_MODEL, near_cutoff_n2=ah_model.near_cutoff_n2,
            )
        except ValueError:
            rejected_mismatch = True
        else:
            rejected_mismatch = False
        compare(metrics, "C07_provider_coordinate_mismatch_rejected", float(not rejected_mismatch), 0.)
        at_zero = MediumSnapshot.from_ah_model(ah_model, POSITION, 0.)
        compare(metrics, "C07_time_changes_density", abs(at_zero.Ne_m3 - 2 * magnetized.Ne_m3), 0.)
        compare(metrics, "C07_time_changes_B", np.linalg.norm(at_zero.B_ecef_T + 2 * magnetized.B_ecef_T), 0.)
        near_cutoff = medium(.99)
        near_roots = appleton_hartree_roots(
            near_cutoff.omega, near_cutoff.Ne_m3, [1, 0, 0], [1, 0, 0], 0,
            plasma_coefficient=near_cutoff.plasma_coefficient,
        )
        compare(metrics, "G1_AH_near_cutoff_status", float(
            select_root(near_roots, "O", near_cutoff_n2=near_cutoff.near_cutoff_n2).status != "near_cutoff"
        ), 0.)
        try:
            modal_ports(near_cutoff, BASIS, 0, 0)
        except ValueError as exc:
            rejected_near_cutoff = "near cutoff" in str(exc)
        else:
            rejected_near_cutoff = False
        compare(metrics, "G1_near_cutoff_port_rejected", float(not rejected_near_cutoff), 0.)
        internal = solve_slab(vacuum, near_cutoff, vacuum, vacuum, BASIS, 0, 0, THICKNESS)
        compare(metrics, "G1_near_cutoff_internal_layer", float(not np.all(np.isfinite(internal.S))), 0.)
        slab = medium(.3, field)
        ku, kv = .2 * K0, .1 * K0
        result = solve_slab(magnetized, slab, magnetized, magnetized, BASIS, ku, kv, THICKNESS)
        compare(metrics, "C07_provider_tensor", np.max(np.abs(magnetized.epsilon() - medium(.2, field).epsilon())), 1e-12)
        uniform = solve_slab(magnetized, magnetized, magnetized, magnetized, BASIS, ku, kv, THICKNESS)
        compare(metrics, "C07_S0_background_recomputed", np.max(np.abs(result.S0 - uniform.S)), 1e-8)
        matrix_properties(result, metrics, "C07")
        for group in (result.left_ports.forward, result.left_ports.backward):
            for port in group:
                roots = appleton_hartree_roots(
                    magnetized.omega, magnetized.Ne_m3, port.k_suv_rad_m,
                    field / np.linalg.norm(field),
                    27_992_489_872.333 * np.linalg.norm(field) / 1e6,
                    plasma_coefficient=COEFFICIENT,
                )
                n2 = np.dot(port.k_suv_rad_m, port.k_suv_rad_m) / K0**2
                compare(metrics, f"C07_AH_root_{port.direction}_{port.label}", abs(roots[port.label] - n2) / abs(n2), 1e-8)
                ah_roots = ah_model.roots(POSITION, port.k_ecef_rad_km, t=TIME)
                wgs_roots = wgs_model.roots(POSITION, port.k_ecef_rad_km, t=TIME)
                compare(metrics, f"C07_AH_model_root_{port.direction}_{port.label}", abs(ah_roots[port.label] - n2) / abs(n2), 1e-8)
                compare(metrics, f"C07_WGS84_AH_model_root_{port.direction}_{port.label}", abs(wgs_roots[port.label] - n2) / abs(n2), 1e-8)
                compare(metrics, f"C07_k_roundtrip_{port.direction}_{port.label}", np.linalg.norm(
                    wavevector_local_rad_m(port.k_ecef_rad_km, BASIS) - port.k_suv_rad_m
                ) / np.linalg.norm(port.k_suv_rad_m), 1e-10)
                compare(metrics, f"C07_E_roundtrip_{port.direction}_{port.label}", np.linalg.norm(
                    BASIS.T @ port.electric_ecef - port.electric_suv
                ) / np.linalg.norm(port.electric_suv), 1e-10)
                compare(metrics, f"C07_H_roundtrip_{port.direction}_{port.label}", np.linalg.norm(
                    BASIS.T @ port.magnetic_ecef - port.magnetic_suv
                ) / np.linalg.norm(port.magnetic_suv), 1e-10)
        legacy = cold_plasma_dielectric_tensor(vacuum.omega, base.Ne_m3, np.zeros(3))
        legacy_explicit = cold_plasma_dielectric_tensor(
            vacuum.omega, base.Ne_m3, field, plasma_coefficient=K_FP**2,
        )
        legacy_implicit = cold_plasma_dielectric_tensor(vacuum.omega, base.Ne_m3, field)
        expected_legacy = 1 - K_FP**2 * base.Ne_m3 / FREQUENCY**2
        compare(metrics, "G1_legacy_coefficient", abs(legacy[0, 0] - expected_legacy), 1e-12)
        compare(metrics, "G1_legacy_explicit_equivalence", np.max(np.abs(legacy_implicit - legacy_explicit)), 0.)
        compare(metrics, "G1_AH_coefficient", abs(base.epsilon()[0, 0] - .8), 1e-12)
        arrays["C07_S"] = result.S
        arrays["C07_S0"] = result.S0
        arrays["C07_E_ecef"] = np.array([p.electric_ecef for p in result.right_ports.forward])
        arrays["C07_H_ecef"] = np.array([p.magnetic_ecef for p in result.right_ports.forward])
        arrays["C07_k_ecef_rad_km"] = np.array([p.k_ecef_rad_km for p in result.right_ports.forward])

        # C08 rotates the complete geometry and tests port rephasing on the fields.
        theta = .71
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                             [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        rotated_base = medium(.2, rotation @ field, position=rotation @ POSITION)
        rotated_slab = medium(.3, rotation @ field, position=rotation @ POSITION)
        rotated = solve_slab(rotated_base, rotated_slab, rotated_base, rotated_base,
                             rotation @ BASIS, ku, kv, THICKNESS)
        compare(metrics, "C08_rotation_S", np.max(np.abs(rotated.S - result.S)), 1e-8)
        compare(metrics, "C08_rotation_S0", np.max(np.abs(rotated.S0 - result.S0)), 1e-8)
        for index, (original, transformed) in enumerate(zip(result.right_ports.forward, rotated.right_ports.forward)):
            compare(metrics, f"C08_rotation_E_{index}", np.linalg.norm(transformed.electric_ecef - rotation @ original.electric_ecef) / np.linalg.norm(original.electric_ecef), 1e-10)
            compare(metrics, f"C08_rotation_H_{index}", np.linalg.norm(transformed.magnetic_ecef - rotation @ original.magnetic_ecef) / np.linalg.norm(original.magnetic_ecef), 1e-10)
            compare(metrics, f"C08_rotation_k_{index}", np.linalg.norm(transformed.k_ecef_rad_km - rotation @ original.k_ecef_rad_km) / np.linalg.norm(original.k_ecef_rad_km), 1e-10)
            compare(metrics, f"C08_E_local_roundtrip_{index}", np.linalg.norm(rotation.T @ transformed.electric_ecef - transformed.electric_suv) / np.linalg.norm(transformed.electric_suv), 1e-10)
            compare(metrics, f"C08_H_local_roundtrip_{index}", np.linalg.norm(rotation.T @ transformed.magnetic_ecef - transformed.magnetic_suv) / np.linalg.norm(transformed.magnetic_suv), 1e-10)
        local_position = np.array([1000.1, -250.5, 310.7])
        back = position_local_m(position_ecef_km(local_position, POSITION, rotation), POSITION, rotation)
        compare(metrics, "C08_position_roundtrip_m", np.linalg.norm(back - local_position), 1e-6)
        phase_in = np.diag(np.exp(1j * np.array([.2, -.7, .5, -.4])))
        phase_out = np.diag(np.exp(1j * np.array([.9, .1, -.3, -.8])))
        incident = np.array([1., .2j, .3, -.5j])
        transformed_s = phase_out.conj().T @ result.S @ phase_in
        original_output = result.S @ incident
        transformed_output = phase_out @ transformed_s @ phase_in.conj().T @ incident
        compare(metrics, "C08_rephasing_fields", np.linalg.norm(transformed_output - original_output) / np.linalg.norm(original_output), 1e-8)
        rephased_output = transformed_s @ phase_in.conj().T @ incident
        for name, port_modes, start in (
            ("left", result.left_ports.backward, 0),
            ("right", result.right_ports.forward, 2),
        ):
            for field_name in ("electric_ecef", "magnetic_ecef"):
                physical = sum(original_output[start + i] * getattr(port, field_name)
                               for i, port in enumerate(port_modes))
                physical_rephased = sum(rephased_output[start + i] * phase_out[start + i, start + i]
                                        * getattr(port, field_name) for i, port in enumerate(port_modes))
                compare(metrics, f"C08_rephasing_{name}_{field_name}",
                        np.linalg.norm(physical_rephased - physical) / np.linalg.norm(physical), 1e-8)

        resource_check()
        failed = [name for name, item in metrics.items() if not item["pass"]]
        status = "PASS" if not failed else "FAIL"
        write_result(metrics, arrays, status, None if not failed else ", ".join(failed))
        print(f"S1 {status}: {len(metrics) - len(failed)}/{len(metrics)} checks")
        for name in failed:
            print(f"FAIL {name}: {metrics[name]}")
        return 0 if not failed else 1
    except Exception as exc:
        write_result(metrics, arrays, "STOPPED", str(exc))
        print(f"S1 STOPPED: {exc}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
