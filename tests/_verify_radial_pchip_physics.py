"""Physical radial PCHIP checks; independent fullpath diagnostics are not a P4 gate."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.io import loadmat
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF, launch_state_ah3d_ecef, spatial_gradient_ah3d_ecef,
    trace_ray_ah3d_ecef,
)
from _a4_c1m_p31_fullpath import RadialFullPath


def rejected(call):
    try:
        call()
    except ValueError:
        return True
    return False


def main(run):
    out = run / 'radial_pchip_validation'
    out.mkdir(parents=True, exist_ok=False)
    grid_path = ROOT / 'output2/a4_pharlap/c1/c1m/p3_20261009_v1/experiments/grids/smooth_s5.mat'
    grid = loadmat(grid_path, simplify_cells=True)
    def provider(kind='radial_pchip', density=None):
        return GeocentricGridNeProvider3D(
            grid['Ne_m3'] if density is None else density,
            grid['latitude_deg'], grid['longitude_deg'], grid['height_km'],
            6371., interpolation=kind)
    p = provider()
    pp = p._radial_spline
    roots = pp.derivative().roots(extrapolate=False)
    q = np.unique(np.r_[p.height_km, roots[np.isfinite(roots)], np.linspace(0, 460, 10001)])
    values = pp(q)
    checks = {
        'nonnegative_over_segments': bool(np.all(values >= 0)),
        'input_nodes_preserved': bool(np.allclose(pp(p.height_km), p.Ne_m3[0, 0], rtol=1e-14, atol=1e-8)),
        'vacuum_entrance_value_and_slope_zero': bool(pp(0) == 0 and pp.derivative()(0) == 0),
        'vacuum_extension': p.get_Ne_geocentric(30, 120, -1) == 0,
        'nonzero_top_not_extended': rejected(lambda: p.get_Ne_geocentric(30, 120, 461)),
        'old_hermite_still_rejects_negative_support': rejected(lambda: provider('radial_hermite')),
        'PH_compatibility_not_silently_expanded': rejected(lambda: RefractiveIndexAH3DECEF(
            p, earth_radius_km=6371., freq_MHz=10., compatibility_mode='pharlap_radial_b0')),
    }
    bad = grid['Ne_m3'].copy(); bad[0, 0, -1] += 1.
    checks['horizontal_variation_rejected'] = rejected(lambda: provider(density=bad))
    bad = grid['Ne_m3'].copy(); bad[0, 0, 0] = -1.
    checks['negative_input_rejected'] = rejected(lambda: provider(density=bad))
    model = RefractiveIndexAH3DECEF(p, earth_radius_km=6371., freq_MHz=10., gyro_MHz=0., collision_s=0., spatial_fd_base_km=.125)
    launch = dict(tx_lat_deg=30., tx_lon_deg=120., tx_height_km=101., elevation_deg=22., bearing_deg=0.)
    state = launch_state_ah3d_ecef(model, **launch)
    exact = model.k0**2 * model.plasma_coefficient / 1e14 * float(pp.derivative()(101)) * state[:3] / np.linalg.norm(state[:3])
    gradient, _ = spatial_gradient_ah3d_ecef(state, model)
    relative = float(np.linalg.norm(gradient - exact) / np.linalg.norm(exact))
    checks['Hamiltonian_FD_matches_same_polynomial_gradient'] = relative < 2e-3
    np.savez_compressed(out / 'medium.npz', height=p.height_km, coefficients=pp.c, query=q, density=values)
    cases = []
    for angle in (22., 39.5):
        for projection in (True, False):
            model = RefractiveIndexAH3DECEF(p, earth_radius_km=6371., freq_MHz=10., gyro_MHz=0., collision_s=0., spatial_fd_base_km=.125, integration_profile='reference')
            launch.update(tx_height_km=0., elevation_deg=angle)
            initial = launch_state_ah3d_ecef(model, **launch)
            ref = RadialFullPath({'h':p.height_km, 'coefficients':pp.c}, initial[:3], initial[3:6]/model.k0, model.plasma_coefficient, 10.)
            ray = trace_ray_ah3d_ecef(model, **launch, h_stop_km=450., max_path_km=3000., apply_dispersion_projection=projection)
            s = np.asarray(ray['trajectory_geom_path_km']); xyz = np.asarray(ray['trajectory_ecef_km'])
            mask = s <= ref.s_end
            truth = ref.query(s[mask])
            error = np.linalg.norm(xyz[mask] - truth['xyz'], axis=1)*1000
            metrics = {'native_rms_m':float(np.sqrt(np.mean(error**2))), 'native_max_m':float(error.max()),
                       'landing_m':float(np.linalg.norm(xyz[-1]-ref.endpoint)*1000),
                       'group_rel':float(abs(ray['group_path_km']-ref.group_end)/ref.group_end)}
            passed = (ray['termination_status'] == 'ground' and metrics['native_rms_m'] <= 2 and
                      metrics['native_max_m'] <= 10 and metrics['landing_m'] <= 10 and metrics['group_rel'] <= 2e-6)
            name = f'e{angle:g}_projection{int(projection)}'
            np.savez_compressed(out/(name+'.npz'), s=s, xyz=xyz, wavevector=ray['wavevector_ecef_km^-1'], group=ray['trajectory_group_path_km'], reference_s=s[mask], reference_xyz=truth['xyz'])
            cases.append(dict(case_id=name, projection=projection, metrics=metrics, diagnostic_passed=bool(passed), termination=ray['termination_status'], coefficient=model.plasma_coefficient, reference_s_end_km=ref.s_end, AH_s_end_km=float(s[-1])))
            print(json.dumps(cases[-1]), flush=True)
    result = dict(checks=checks, gradient_relative_error=relative, cases=cases,
                  engineering_and_diagnostics_passed=bool(all(checks.values()) and all(c['diagnostic_passed'] for c in cases)),
                  formal_physical_acceptance_passed=False, recommend_p4=False,
                  limitation='Native-node diagnostics only; full convergence, reconstruction and FD/grid budgets remain required.',
                  input_sha256=hashlib.sha256(grid_path.read_bytes()).hexdigest(), source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    if not result['engineering_and_diagnostics_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--run', required=True, type=Path)
    main(parser.parse_args().run.resolve())
