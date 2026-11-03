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
from _a4_c1m_p31_fullpath import (RadialFullPath, LocalRadialPath,
    polynomial_precision_report, cartesian_control, analytic_control_report)
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline


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


PHYSICAL_LIMITS = {'rms_m': 2., 'max_m': 10., 'landing_m': 10., 'group_rel': 2e-6}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def physics_identity(run):
    manifest = json.loads((run/'manifest.json').read_text(encoding='utf-8'))
    for record in manifest['source_records']:
        source = ROOT/record['path']
        if hashlib.sha256(source.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError('Execution source changed after freeze: '+record['path'])
    for record in manifest['input_records']:
        if hashlib.sha256((run/record['path']).read_bytes()).hexdigest() != record['sha256']:
            raise ValueError('Frozen input changed: '+record['path'])
    return manifest


def make_provider(path):
    with np.load(path) as data:
        h, ne = data['h'], data['Ne']
        lat, lon = data['lat'], data['lon']
        if ne.ndim == 1:
            ne = np.broadcast_to(ne, (len(lat), len(lon), len(h))).copy()
        provider = GeocentricGridNeProvider3D(ne, lat, lon, h, 6371., interpolation='radial_pchip')
        if not np.array_equal(provider._radial_spline.c, data['coefficients']):
            raise ValueError('Frozen and reconstructed PCHIP coefficients differ')
        field = {'h': h.copy(), 'coefficients': data['coefficients'].copy()}
    return provider, field


def physics_trace(case, provider, profile, coefficient, fd=.125, projection=True):
    from config import UPDATE2
    from models.ray_tracer_ah_3d_ecef import integrate_fixed_path_ah3d_ecef
    key = '_p31_frozen'
    previous = UPDATE2['AH3D_INTEGRATION'].get(key)
    UPDATE2['AH3D_INTEGRATION'][key] = dict(profile)
    try:
        model = RefractiveIndexAH3DECEF(provider, earth_radius_km=6371.,
            freq_MHz=case['frequency_MHz'], gyro_MHz=0., collision_s=0.,
            compatibility_mode=None, plasma_coefficient=coefficient, spatial_fd_base_km=fd)
        launch = dict(tx_lat_deg=30., tx_lon_deg=120., tx_height_km=case['origin_height_km'],
                      elevation_deg=case['elevation_deg'], bearing_deg=0.)
        initial = launch_state_ah3d_ecef(model, **launch)
        local = case['scope'] == 'first_100_km'
        ray = integrate_fixed_path_ah3d_ecef(model, **launch, profile=key,
            path_length_km=100. if local else case['max_path_km'],
            h_stop_km=case['h_stop_km'], stop_on_physical_events=not local,
            apply_dispersion_projection=projection)
        state = {'s': np.asarray(ray['trajectory_geom_path_km']),
                 'xyz': np.asarray(ray['trajectory_ecef_km']),
                 'p': np.asarray(ray['wavevector_ecef_km^-1'])/model.k0,
                 'group': np.asarray(ray['trajectory_group_path_km']),
                 'tangent': np.asarray(ray['wavevector_ecef_km^-1'])/np.linalg.norm(ray['wavevector_ecef_km^-1'], axis=1)[:, None],
                 'initial': initial, 'termination': ray['termination_status']}
        return state, model
    finally:
        if previous is None:
            UPDATE2['AH3D_INTEGRATION'].pop(key)
        else:
            UPDATE2['AH3D_INTEGRATION'][key] = previous


def save_trace(path, trace):
    np.savez_compressed(path, **trace)


def reference_for(case, field, trace, coefficient, eps=2e-10):
    constructor = LocalRadialPath if case['scope'] == 'first_100_km' else RadialFullPath
    return constructor(field, trace['xyz'][0], trace['p'][0], coefficient, case['frequency_MHz'], eps=eps)


def splines(trace):
    return (CubicHermiteSpline(trace['s'], trace['xyz'], trace['tangent']),
            CubicHermiteSpline(trace['s'], trace['group'], 1/np.linalg.norm(trace['p'], axis=1)))


def metrics(trace, reference, q):
    xyz, group = splines(trace)
    truth = reference.query(q)
    d = np.linalg.norm(xyz(q)-truth['xyz'], axis=1)*1000
    return {'rms_m': float(np.sqrt(np.mean(d*d))), 'max_m': float(d.max()),
            'landing_m': float(np.linalg.norm(trace['xyz'][-1]-reference.endpoint)*1000),
            'group_rel': float(abs(trace['group'][-1]-reference.group_end)/reference.group_end)}


def cross_metrics(a, b):
    end = min(float(a['s'][-1]), float(b['s'][-1]))
    q = np.unique(np.r_[np.arange(0., end, .5), end])
    axyz, ag = splines(a); bxyz, bg = splines(b)
    d = np.linalg.norm(axyz(q)-bxyz(q), axis=1)*1000
    return {'rms_m': float(np.sqrt(np.mean(d*d))), 'max_m': float(d.max()),
            'landing_m': float(np.linalg.norm(a['xyz'][-1]-b['xyz'][-1])*1000),
            'group_rel': float(abs(a['group'][-1]-b['group'][-1])/b['group'][-1])}


def medium_and_state(case, provider, field, trace, model):
    pp = provider._radial_spline
    critical = pp.derivative().roots(extrapolate=False)
    heights = np.unique(np.r_[field['h'], critical[np.isfinite(critical)]])
    h = np.linalg.norm(trace['xyz'], axis=1)-6371.
    ne = pp(np.clip(h, field['h'][0], field['h'][-1]))
    residual = np.sum(trace['p']**2, axis=1)-(1-model.plasma_coefficient*ne/(model.freq_MHz*1e6)**2)
    angular = np.cross(trace['xyz'], trace['p'])
    relative = np.linalg.norm(angular-angular[0], axis=1)/np.linalg.norm(angular[0])
    normal = angular[0]/np.linalg.norm(angular[0])
    plane_m = float(np.max(np.abs(trace['xyz']@normal))*1000)
    local = case['scope'] == 'first_100_km'
    coverage = (abs(trace['s'][-1]-100.) < 1e-8 if local else trace['termination'] == 'ground')
    checks = {'nonnegative_segments': bool(np.min(pp(heights)) >= 0),
              'node_reproduction': bool(np.allclose(pp(field['h']), provider.Ne_m3[0, 0], rtol=1e-14, atol=1e-8)),
              'smooth_ground_entry': bool(local or (pp(0.) == 0 and pp.derivative()(0.) == 0)),
              'coverage_termination': bool(coverage),
              'dispersion_residual': bool(np.max(np.abs(residual)) <= 1e-10),
              'initial_dispersion': bool(abs(residual[0]) <= 1e-10),
              'plane': bool(plane_m <= .2),
              'radial_invariant': bool(np.max(relative) <= 2e-6),
              'monotonic_arc_group': bool(np.all(np.diff(trace['s']) > 0) and np.all(np.diff(trace['group']) > 0))}
    return {'checks': checks, 'passed': all(checks.values()),
            'residual_threshold_source': 'production initialization/projection tolerance 1e-10',
            'residual_max': float(np.max(np.abs(residual))), 'plane_max_m': plane_m,
            'angular_invariant_relative_max': float(np.max(relative)),
            'boundary_density': [float(pp(field['h'][i])) for i in (0, -1)],
            'boundary_slope': [float(pp.derivative()(field['h'][i])) for i in (0, -1)],
            'segment_density_min': float(np.min(pp(heights))),
            'initial_p': trace['p'][0].tolist(), 'termination': trace['termination']}


def physics_case(run, case, manifest):
    registration = json.loads((run/'reference_registration.json').read_text(encoding='utf-8'))
    registered = next(row for row in registration['cases'] if row['case_id'] == case['case_id'])
    if not registered['reference_registered']:
        raise ValueError('Independent reference topology not registered before AH')
    folder = run/'cases'/case['case_id']
    if folder.exists():
        raise FileExistsError('Refuse case overwrite: '+str(folder))
    folder.mkdir(parents=True)
    provider, field = make_provider(run/case['medium_path'])
    traces = []
    for i, profile in enumerate(manifest['integration_profiles'][:3]):
        trace, model = physics_trace(case, provider, profile, manifest['coefficient'],
                                     projection=manifest['projection'])
        save_trace(folder/f'ah_l{i+1}.npz', trace); traces.append(trace)
    final = traces[-1]
    state = medium_and_state(case, provider, field, final, model)
    write_json(folder/'medium_checks.json', state)
    ref = reference_for(case, field, final, manifest['coefficient'])
    if not np.array_equal(final['xyz'][0], registered['initial_xyz']) or not np.array_equal(final['p'][0], registered['initial_p']):
        raise ValueError('AH initial state differs from independent reference registration')
    if case['scope'] != 'first_100_km' and abs(ref.s_end-final['s'][-1])*1000 > 10.:
        raise ValueError('Natural path coverage differs by more than 10 m')
    precision = polynomial_precision_report(ref)
    write_json(folder/'reference_precision.json', precision)
    end = min(float(final['s'][-1]), ref.s_end)
    iq = np.unique(np.r_[np.arange(0., end, .25), end])
    native = final['s'][final['s'] <= end]
    mids = (native[:-1]+native[1:])/2
    q = np.unique(np.r_[iq, native, mids])
    truth = ref.query(q)
    np.savez_compressed(folder/'reference.npz', **truth)
    control = cartesian_control(ref)
    control_end = min(float(control.t[-1]), ref.s_end)
    control_q = np.unique(np.r_[np.linspace(0., control_end, 201),
                               ref.half[0] if hasattr(ref, 'half') else control_end/2])
    control_truth = ref.query(control_q)
    control_error = np.linalg.norm(control.sol(control_q)[:3].T-control_truth['xyz'], axis=1)*1000
    control_group = float(np.max(abs(control.sol(control_q)[6]-control_truth['group']))/ref.group_end)
    control_landing = float(np.linalg.norm(control.y[:3, -1]-ref.endpoint)*1000)
    control_total_group = float(abs(control.y[6, -1]-ref.group_end)/ref.group_end)
    precision['cartesian_control_max_m'] = float(control_error.max())
    precision['cartesian_control_group_rel'] = control_group
    precision['cartesian_natural_endpoint_m'] = control_landing
    precision['cartesian_natural_group_rel'] = control_total_group
    precision['position_error_estimate_m'] += max(float(control_error.max()), control_landing)
    precision['group_error_estimate_rel'] += max(control_group, control_total_group)
    write_json(folder/'reference_precision.json', precision)
    np.savez_compressed(folder/'cartesian_control.npz', s=control_q, state=control.sol(control_q),
                        reference_xyz=control_truth['xyz'], reference_group=control_truth['group'],
                        natural_s_end=float(control.t[-1]), natural_endpoint=control.y[:3, -1],
                        natural_group_end=float(control.y[6, -1]),
                        termination='fixed_geometric_path_completed' if case['scope'] == 'first_100_km' else 'ground')
    xyz, group = splines(final)
    ix = np.searchsorted(q, native)
    native_error = np.linalg.norm(final['xyz'][:len(native)]-truth['xyz'][ix], axis=1)*1000
    refxyz = CubicHermiteSpline(native, truth['xyz'][ix], truth['tangent'][ix])
    refgroup = CubicHermiteSpline(native, truth['group'][ix], 1/np.linalg.norm(truth['p'][ix], axis=1))
    mx = np.searchsorted(q, mids)
    reconstruction = {'position_m': float(np.max(np.linalg.norm(refxyz(mids)-truth['xyz'][mx], axis=1))*1000),
                      'group_rel': float(np.max(abs(refgroup(mids)-truth['group'][mx]))/ref.group_end)}
    samples = []
    for ds in (1., .5, .25):
        # All three grids use a common independently evaluated reference array.
        query = np.unique(np.r_[np.arange(0., end, ds), end])
        index = np.searchsorted(q, query)
        error = np.linalg.norm(xyz(query)-truth['xyz'][index], axis=1)*1000
        samples.append({'step_km': ds, 'rms_m': float(np.sqrt(np.mean(error**2))), 'max_m': float(error.max())})
    sampling = {k: max(abs(samples[i][k]-samples[i-1][k]) for i in (1, 2)) for k in ('rms_m', 'max_m')}
    sampling.update(landing_m=0., group_rel=0.)
    metrics_final = dict(samples[-1]); metrics_final.pop('step_km')
    metrics_final.update(landing_m=float(np.linalg.norm(final['xyz'][-1]-ref.endpoint)*1000),
                         group_rel=float(abs(final['group'][-1]-ref.group_end)/ref.group_end))
    propagation = {k: float(native_error.max()) for k in ('rms_m', 'max_m')}
    propagation.update(landing_m=metrics_final['landing_m'], group_rel=metrics_final['group_rel'])
    rec = {k: reconstruction['position_m'] for k in ('rms_m', 'max_m')}
    rec.update(landing_m=0., group_rel=reconstruction['group_rel'])
    reference_error = {k: precision['position_error_estimate_m'] for k in ('rms_m', 'max_m', 'landing_m')}
    reference_error['group_rel'] = precision['group_error_estimate_rel']
    total = {k: propagation[k]+rec[k]+sampling[k]+reference_error[k] for k in PHYSICAL_LIMITS}
    increments = [cross_metrics(traces[i-1], traces[i]) for i in (1, 2)]
    monotonic = all(increments[1][k] <= increments[0][k] for k in PHYSICAL_LIMITS)
    bounded = all(increments[-1][k] <= v for k, v in PHYSICAL_LIMITS.items())
    strict = (monotonic and bounded and all(metrics_final[k] <= v for k, v in PHYSICAL_LIMITS.items())
              and all(metrics_final[k] <= 3*(increments[-1][k]+rec[k]+sampling[k]+reference_error[k]) for k in PHYSICAL_LIMITS))
    budget_passed = (all(total[k] <= v for k, v in PHYSICAL_LIMITS.items())
                    and all(part[k] <= .1*v for part in (rec, sampling, reference_error) for k, v in PHYSICAL_LIMITS.items()))
    fd_traces = []
    for fd in (.25, .125, .0625):
        if fd == .125:
            trace = final
        else:
            trace, _ = physics_trace(case, provider, manifest['integration_profiles'][2], manifest['coefficient'],
                                      fd=fd, projection=manifest['projection'])
            save_trace(folder/f'ah_fd{fd:g}.npz', trace)
        fd_traces.append(trace)
    fd_increments = [cross_metrics(fd_traces[i-1], fd_traces[i]) for i in (1, 2)]
    fd_passed = (all(fd_increments[1][k] <= .2*v for k, v in PHYSICAL_LIMITS.items())
                 and all(fd_increments[1][k] <= fd_increments[0][k] for k in PHYSICAL_LIMITS))
    # A nonmonotone FD platform needs derivative evidence, not a small-error exemption.
    derivative_rows = []
    for fd in (.25, .125, .0625):
        fd_model = RefractiveIndexAH3DECEF(provider, earth_radius_km=6371., freq_MHz=case['frequency_MHz'],
            gyro_MHz=0., collision_s=0., plasma_coefficient=manifest['coefficient'], spatial_fd_base_km=fd)
        errs = []
        for index in np.linspace(0, len(final['s'])-1, min(21, len(final['s'])), dtype=int):
            x = final['xyz'][index]; height = np.linalg.norm(x)-6371.
            if height < field['h'][0] or height > field['h'][-1]:
                continue
            actual, _ = spatial_gradient_ah3d_ecef(np.r_[x, final['p'][index]*fd_model.k0, final['group'][index]], fd_model)
            exact = fd_model.k0**2*fd_model.plasma_coefficient/(fd_model.freq_MHz*1e6)**2*float(provider._radial_spline.derivative()(height))*x/np.linalg.norm(x)
            errs.append(float(np.linalg.norm(actual-exact)/fd_model.k0**2))
        derivative_rows.append({'fd_km': fd, 'normalized_gradient_error_max': max(errs)})
    np.savetxt(folder/'point_metrics.csv', np.c_[q, np.linalg.norm(xyz(q)-truth['xyz'], axis=1)*1000,
               group(q)-truth['group']], delimiter=',', header='s_km,position_error_m,group_error_km', comments='')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(q, np.linalg.norm(xyz(q)-truth['xyz'], axis=1)*1000)
    axes[0].set(xlabel='Absolute geometric path (km)', ylabel='Position error (m)')
    axes[1].plot(final['s'], np.linalg.norm(final['xyz'], axis=1)-6371.)
    axes[1].set(xlabel='Geometric path (km)', ylabel='Height (km)')
    figure.tight_layout(); figure.savefig(folder/'trajectory_errors.png', dpi=140); plt.close(figure)
    result = {'case_id': case['case_id'], 'scope': case['scope'], 'metrics': metrics_final,
              'medium_initial_state_passed': state['passed'], 'strict_convergence_passed': bool(strict),
              'reference_certification_passed': bool(budget_passed), 'reconstruction_budget_passed': bool(budget_passed),
              'passed': bool(state['passed'] and (strict or budget_passed)),
              'acceptance_route': 'strict' if strict else 'independent_reference' if budget_passed else 'failed',
              'propagation_envelope': propagation, 'reconstruction': rec, 'sampling': sampling,
              'reference_error': reference_error, 'total_error_estimate': total,
              'sampling_levels': samples, 'integration_increments': increments,
              'fd_increments': fd_increments, 'fd_passed': bool(fd_passed), 'fd_derivative_controls': derivative_rows,
              'reference_natural_s_end_km': ref.s_end, 'reference_natural_group_km': ref.group_end,
              'reference_endpoint': ref.endpoint.tolist(), 'reference_topology': 'local_upward' if case['scope'] == 'first_100_km' else 'single_turn_ground',
              'reference_error_is_interval_proof': False}
    write_json(folder/'result.json', result)
    print(json.dumps({'case_id': case['case_id'], 'passed': result['passed'], 'fd_passed': result['fd_passed'], 'metrics': metrics_final}), flush=True)
    return result


def physics_execute(run, only=None):
    manifest = physics_identity(run)
    registry = json.loads((run/'case_registry.json').read_text(encoding='utf-8'))
    records = []
    for case in registry['cases']:
        if only and case['case_id'] != only:
            continue
        result_path = run/'cases'/case['case_id']/'result.json'
        if result_path.exists():
            # Raw data hash verification is performed by the independent verifier.
            raise FileExistsError('Use a new subrun; execution never silently reuses results')
        try:
            records.append(physics_case(run, case, manifest))
        except Exception as error:
            folder = run/'cases'/case['case_id']; folder.mkdir(parents=True, exist_ok=True)
            record = {'case_id': case['case_id'], 'passed': False, 'status': 'tool_or_physical_failure',
                      'exception_type': type(error).__name__, 'reason': str(error)}
            write_json(folder/'failure.json', record); records.append(record)
            print(json.dumps(record), flush=True)
    write_json(run/('execution_'+only+'.json' if only else 'execution.json'), records)
    return records


def physics_register(run):
    manifest = physics_identity(run)
    registry = json.loads((run/'case_registry.json').read_text(encoding='utf-8'))
    if (run/'reference_registration.json').exists() or (run/'cases').exists():
        raise FileExistsError('Reference registration must precede every AH execution')
    rows = []
    for case in registry['cases']:
        provider, field = make_provider(run/case['medium_path'])
        model = RefractiveIndexAH3DECEF(provider, earth_radius_km=6371., freq_MHz=case['frequency_MHz'],
            gyro_MHz=0., collision_s=0., compatibility_mode=None, plasma_coefficient=manifest['coefficient'],
            spatial_fd_base_km=.125)
        initial = launch_state_ah3d_ecef(model, tx_lat_deg=30., tx_lon_deg=120.,
            tx_height_km=case['origin_height_km'], elevation_deg=case['elevation_deg'], bearing_deg=0.)
        trace = {'xyz': initial[None, :3], 'p': initial[None, 3:6]/model.k0}
        try:
            reference = reference_for(case, field, trace, manifest['coefficient'])
            rows.append({'case_id': case['case_id'], 'reference_registered': True,
                'initial_xyz': trace['xyz'][0].tolist(), 'initial_p': trace['p'][0].tolist(),
                'topology': 'local_upward' if case['scope'] == 'first_100_km' else 'single_turn_ground',
                'natural_s_end': reference.s_end, 'natural_group_end': reference.group_end,
                'natural_endpoint': reference.endpoint.tolist(),
                'apex_km': None if case['scope'] == 'first_100_km' else reference.rt-6371.})
        except Exception as error:
            rows.append({'case_id': case['case_id'], 'reference_registered': False, 'reason': str(error)})
    report = {'passed': all(r['reference_registered'] for r in rows), 'cases': rows,
              'AH_execution_started': False, 'purpose': 'topology registered after frozen inputs, before AH'}
    write_json(run/'reference_registration.json', report)
    print(json.dumps({'reference_registered': sum(r['reference_registered'] for r in rows), 'case_count': len(rows)}), flush=True)
    if not report['passed']:
        raise ValueError('Reference preregistration has unsupported topology')
    return report


def physics_controls(run):
    out = run/'reference_controls'
    out.mkdir(parents=True, exist_ok=False)
    report = analytic_control_report(out)
    checks = {'straight': report['checks']['straight_reconstruction'],
              'curved_cartesian': all(report['cartesian_agreement'][k] <= .1*v for k, v in PHYSICAL_LIMITS.items()),
              'turn_neighborhood': report['checks']['turn_point_finite'],
              'early_termination': report['checks']['early_termination_rejected'],
              'truncated_prefix': report['checks']['truncated_common_prefix_rejected'],
              'natural_endpoint': report['checks']['natural_endpoint_mismatch_detected'],
              'natural_group': report['checks']['natural_group_mismatch_detected']}
    from config import K_FP
    coefficient = K_FP**2
    h = np.array([0., 460.])
    f = {'h': h, 'coefficients': np.zeros((4, 1))}
    local = LocalRadialPath(f, [6471., 0., 0.], [.5, np.sqrt(.75), 0.], coefficient, 10.)
    truth = local.query([0., 10., 100.])
    expected = np.array([6471., 0., 0.])+np.array([0., 10., 100.])[:, None]*np.array([.5, np.sqrt(.75), 0.])
    checks['local_vacuum'] = bool(np.max(abs(truth['xyz']-expected)) < 1e-9 and abs(local.group_end-100.) < 1e-9)
    nonzero = {'h': h, 'coefficients': np.array([[0.], [0.], [1e8], [1e11]])}
    n = np.sqrt(1-coefficient*1.1e11/1e14)
    p = np.array([.5, np.sqrt(.75), 0.])*n
    checks['wrong_initial_state'] = rejected(lambda: LocalRadialPath(nonzero, [6471., 0., 0.], p*1.001, coefficient, 10.))
    checks['wrong_frequency_units'] = rejected(lambda: LocalRadialPath(nonzero, [6471., 0., 0.], p, coefficient, 10e6))
    increments = [{'rms_m': .01, 'max_m': .02, 'landing_m': .01, 'group_rel': 1e-10},
                  {'rms_m': .02, 'max_m': .03, 'landing_m': .02, 'group_rel': 2e-10}]
    checks['nonmonotonic_not_strict'] = not all(increments[1][k] <= increments[0][k] for k in PHYSICAL_LIMITS)
    report['physical_checks'] = checks
    report['physical_passed'] = bool(all(checks.values()))
    report['new_physical_limits'] = PHYSICAL_LIMITS
    report['local_precision'] = polynomial_precision_report(local)
    write_json(out/'physical_controls.json', report)
    if not report['physical_passed']:
        raise ValueError('Reference positive/negative controls failed')
    return report


def continuous_control(case, kind, coefficient, gradient=None, max_step=.25):
    """Independent Cartesian integration of the frozen node-generating formula."""
    provider, _ = make_provider(case['_medium_absolute_path'])
    model = RefractiveIndexAH3DECEF(provider, earth_radius_km=6371., freq_MHz=case['frequency_MHz'],
                                  gyro_MHz=0., collision_s=0., plasma_coefficient=coefficient)
    initial = launch_state_ah3d_ecef(model, tx_lat_deg=30., tx_lon_deg=120.,
        tx_height_km=case['origin_height_km'], elevation_deg=case['elevation_deg'], bearing_deg=0.)
    alpha = coefficient/(case['frequency_MHz']*1e6)**2
    def density(h):
        if kind == 'affine':
            return 1e11+gradient*(h-100.), gradient
        u = (h-350.)/100.
        ne = 607635874822.0494*np.exp(.5*(1-u-np.exp(-u)))
        derivative = ne*.5*(np.exp(-u)-1)/100.
        if kind == 'smooth':
            x = float(np.clip((h-20.)/60., 0., 1.))
            taper = x**3*(10+x*(-15+6*x))
            slope = 30*x*x*(1-x)**2/60. if 0 < x < 1 else 0.
            derivative = derivative*taper+ne*slope
            ne *= taper
        return float(ne), float(derivative)
    def rhs(s, y):
        r = np.linalg.norm(y[:3]); ne, dne = density(r-6371.)
        n = np.sqrt(1-alpha*ne)
        return np.r_[y[3:6]/n, -alpha*dne*y[:3]/(2*n*r), 1/n]
    local = case['scope'] == 'first_100_km'
    def ground(s, y):
        return np.linalg.norm(y[:3])-6371. if s > 1e-7 else 1e-7
    ground.terminal = True; ground.direction = -1
    ode = solve_ivp(rhs, (0., 100. if local else case['max_path_km']),
        np.r_[initial[:3], initial[3:6]/model.k0, 0.], method='DOP853',
        rtol=2e-12, atol=2e-14, max_step=max_step, dense_output=True,
        events=None if local else ground)
    if not ode.success or (not local and not len(ode.t_events[0])):
        raise ValueError('Continuous generator reference did not complete prescribed/natural path')
    return ode


def physics_grid(run):
    manifest = physics_identity(run)
    registry = json.loads((run/'case_registry.json').read_text(encoding='utf-8'))
    out = run/'grid_controls'; out.mkdir(exist_ok=False)
    records = []
    for case in registry['cases']:
        target = out/case['case_id']; target.mkdir()
        medium = registry['media'][case['medium_id']]['continuous_background']
        kind = medium['kind']
        analytic_case = dict(case, _medium_absolute_path=run/case['medium_path'])
        try:
            analytic = continuous_control(analytic_case, kind, manifest['coefficient'], medium.get('gradient_m3_per_km'))
            refined = continuous_control(analytic_case, kind, manifest['coefficient'], medium.get('gradient_m3_per_km'), .125)
            q = np.linspace(0., min(analytic.t[-1], refined.t[-1]), 201)
            analytic_error = float(np.max(np.linalg.norm(analytic.sol(q)[:3].T-refined.sol(q)[:3].T, axis=1))*1000)
            analytic_group_error = float(abs(analytic.y[6, -1]-refined.y[6, -1])/refined.y[6, -1])
            np.savez_compressed(target/'continuous_reference.npz', s=q, state=refined.sol(q),
                coarse_state=analytic.sol(q), natural_s_end=float(refined.t[-1]),
                natural_endpoint=refined.y[:3, -1], natural_group_end=float(refined.y[6, -1]))
            levels = []
            for variant in case['grid_variants']:
                path = variant['path']
                provider, field = make_provider(run/path)
                trace, model = physics_trace(case, provider, manifest['integration_profiles'][2],
                                              manifest['coefficient'], projection=manifest['projection'])
                spacing = float(field['h'][1]-field['h'][0])
                save_trace(target/f'ah_d{spacing:g}.npz', trace)
                reference = reference_for(case, field, trace, manifest['coefficient'])
                end = min(reference.s_end, trace['s'][-1], refined.t[-1])
                q = np.unique(np.r_[np.arange(0., end, .5), end])
                actual = metrics(trace, reference, q)
                truth = reference.query(q)
                analytic_xyz = refined.sol(q)[:3].T
                difference = np.linalg.norm(truth['xyz']-analytic_xyz, axis=1)*1000
                grid_error = {'rms_m': float(np.sqrt(np.mean(difference**2))), 'max_m': float(difference.max()),
                    'landing_m': float(np.linalg.norm(reference.endpoint-refined.y[:3, -1])*1000),
                    'group_rel': float(abs(reference.group_end-refined.y[6, -1])/refined.y[6, -1])}
                np.savez_compressed(target/f'grid_reference_d{spacing:g}.npz', s=q,
                                    xyz=truth['xyz'], group=truth['group'], continuous_xyz=analytic_xyz)
                levels.append({'spacing_km': spacing, 'medium_path': path, 'AH_to_actual_medium': actual,
                    'actual_medium_to_continuous_generator': grid_error,
                    'node_and_gradient_reproduction_max': None if kind != 'affine' else
                        float(np.max(abs(provider._radial_spline.derivative()(field['h'])-medium['gradient_m3_per_km']))),
                    'passed': bool(all(actual[k] <= v and grid_error[k] <= v for k, v in PHYSICAL_LIMITS.items()))})
            baseline = next(row for row in levels if row['medium_path'] == case['medium_path'])
            record = {'case_id': case['case_id'], 'levels': levels, 'passed': all(r['passed'] for r in levels),
                'baseline_passed': baseline['passed'], 'baseline_spacing_km': baseline['spacing_km'],
                'continuous_control_refinement_max_m': analytic_error,
                'continuous_control_refinement_group_rel': analytic_group_error,
                'continuous_control_precision_passed': bool(analytic_error <= .2 and analytic_group_error <= 2e-7)}
            record['passed'] = bool(record['passed'] and record['continuous_control_precision_passed'])
        except Exception as error:
            record = {'case_id': case['case_id'], 'passed': False, 'reason': str(error), 'exception_type': type(error).__name__}
        write_json(target/'result.json', record); records.append(record)
        print(json.dumps({'grid_case': case['case_id'], 'passed': record['passed'], 'reason': record.get('reason')}), flush=True)
    projection = []
    template = next(c for c in registry['cases'] if c['medium_id'] == 'smooth_s5')
    provider, _ = make_provider(run/template['medium_path'])
    for angle in (22., 39.5):
        case = dict(template, elevation_deg=angle)
        traces = []
        for enabled in (True, False):
            trace, _ = physics_trace(case, provider, manifest['integration_profiles'][2],
                                      manifest['coefficient'], projection=enabled)
            save_trace(out/f'projection_e{angle:g}_{int(enabled)}.npz', trace); traces.append(trace)
        response = cross_metrics(*traces)
        projection.append({'elevation_deg': angle, 'response': response,
                           'passed': all(response[k] <= .2*v for k, v in PHYSICAL_LIMITS.items())})
    report = {'cases': records, 'projection_controls': projection,
              'passed': bool(all(r['passed'] for r in records) and all(r['passed'] for r in projection)),
              'scope': 'separate AH integration error and actual PCHIP medium discretization; frozen baseline retained'}
    write_json(out/'summary.json', report)
    return report


def physics_preconditions(run):
    """Diagnose unsupported registered topology without executing a production ray."""
    import mpmath as mp
    from numpy.polynomial import Polynomial
    manifest = physics_identity(run)
    registry = json.loads((run/'case_registry.json').read_text(encoding='utf-8'))
    registration = json.loads((run/'reference_registration.json').read_text(encoding='utf-8'))
    out = run/'precondition_checks'; out.mkdir(exist_ok=False)
    records = []
    for row in registration['cases']:
        if row['reference_registered']:
            continue
        case = next(c for c in registry['cases'] if c['case_id'] == row['case_id'])
        provider, field = make_provider(run/case['medium_path'])
        model = RefractiveIndexAH3DECEF(provider, earth_radius_km=6371., freq_MHz=case['frequency_MHz'],
                                      gyro_MHz=0., collision_s=0., plasma_coefficient=manifest['coefficient'])
        initial = launch_state_ah3d_ecef(model, tx_lat_deg=30., tx_lon_deg=120.,
            tx_height_km=case['origin_height_km'], elevation_deg=case['elevation_deg'], bearing_deg=0.)
        x, p = initial[:3], initial[3:6]/model.k0
        r0 = np.linalg.norm(x); er = x/r0
        J = r0*np.linalg.norm(p-(p@er)*er)
        alpha = manifest['coefficient']/(case['frequency_MHz']*1e6)**2
        E = float(p@p+alpha*provider._radial_spline(r0-6371.))
        rows = []
        with mp.workdps(80):
            exact = lambda v: mp.mpf(float(v).as_integer_ratio()[0])/float(v).as_integer_ratio()[1]
            for i, (a, b) in enumerate(zip(field['h'][:-1], field['h'][1:])):
                radius = Polynomial([6371.+a, 1.])
                ne = Polynomial(field['coefficients'][:, i][::-1])
                q = (Polynomial([E])-alpha*ne)*radius**2-J**2
                candidates = [0., b-a]+[float(v.real) for v in q.deriv().roots()
                    if abs(v.imag) < 1e-8 and 0 < v.real < b-a]
                for t in candidates:
                    coeff = [exact(v) for v in field['coefficients'][:, i]]
                    tm, rm = exact(t), exact(6371.)+exact(a)+exact(t)
                    nem = ((coeff[0]*tm+coeff[1])*tm+coeff[2])*tm+coeff[3]
                    value = (exact(E)-exact(alpha)*nem)*rm**2-exact(J)**2
                    rows.append({'height_km': float(a+t), 'Q_km2_80dps': mp.nstr(value, 80),
                                 'radial_p_squared': float(value/rm**2), 'segment': i})
        minimum = min(rows, key=lambda v: v['radial_p_squared'])
        def rhs(s, y):
            r = np.linalg.norm(y[:3]); h = r-6371.
            ne = float(provider._radial_spline(h, extrapolate=True))
            dne = float(provider._radial_spline.derivative()(h, extrapolate=True))
            n = np.sqrt(E-alpha*ne)
            return np.r_[y[3:6]/n, -alpha*dne*y[:3]/(2*n*r), 1/n]
        def upper(s, y):
            return np.linalg.norm(y[:3])-6371.-case['h_stop_km']
        upper.terminal = True; upper.direction = 1
        controls = []
        for step in (.5, .25):
            ode = solve_ivp(rhs, (0., case['max_path_km']), np.r_[x, p, 0.], method='DOP853',
                rtol=2e-12, atol=2e-14, max_step=step, events=upper, dense_output=True)
            controls.append(ode)
            np.savez_compressed(out/f"{case['case_id']}_cartesian_step{step:g}.npz", s=ode.t, state=ode.y,
                                upper_event=ode.t_events[0])
        last = controls[-1]
        record = {'case_id': case['case_id'], 'reason': row['reason'],
            'minimum_radial_p_squared': minimum['radial_p_squared'], 'minimum': minimum,
            'all_segment_stationary_evaluations': rows,
            'method': 'every segment endpoint and real stationary root; exact binary64 conversion at 80 dps; numerical evidence, no interval claim',
            'cartesian_upper_reached': bool(last.success and len(last.t_events[0]) == 1),
            'cartesian_endpoint_height_km': float(np.linalg.norm(last.y[:3, -1])-6371.),
            'cartesian_endpoint_radial_p': float(last.y[3:6, -1]@last.y[:3, -1]/np.linalg.norm(last.y[:3, -1])),
            'cartesian_s_end_km': float(last.t[-1]), 'cartesian_group_end_km': float(last.y[6, -1]),
            'cartesian_refinement_endpoint_m': float(np.linalg.norm(controls[0].y[:3, -1]-last.y[:3, -1])*1000),
            'cartesian_refinement_group_rel': float(abs(controls[0].y[6, -1]-last.y[6, -1])/last.y[6, -1]),
            'AH_executed': False, 'classification': 'outgoing_upper_boundary' if len(last.t_events[0]) else 'unresolved'}
        records.append(record)
        print(json.dumps({k: record[k] for k in ('case_id', 'minimum_radial_p_squared', 'cartesian_upper_reached', 'classification')}), flush=True)
    report = {'blocking_cases': records, 'passed': not records,
              'formal_matrix_executed': False, 'recommend_p4_physics': False}
    write_json(out/'preconditions.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--stage', choices=('diagnostic', 'controls', 'register', 'preconditions', 'execute', 'grid'), default='diagnostic')
    parser.add_argument('--case')
    args = parser.parse_args()
    if args.stage == 'diagnostic':
        main(args.run.resolve())
    elif args.stage == 'controls':
        physics_controls(args.run.resolve())
    elif args.stage == 'register':
        physics_register(args.run.resolve())
    elif args.stage == 'grid':
        physics_grid(args.run.resolve())
    elif args.stage == 'preconditions':
        physics_preconditions(args.run.resolve())
    else:
        physics_execute(args.run.resolve(), args.case)
