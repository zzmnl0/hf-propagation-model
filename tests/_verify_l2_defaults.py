"""Current formal L2 routing and bounded real short-ray regression."""
from __future__ import annotations

import argparse
import ast
from contextlib import ExitStack
import json
import math
import os
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import config as cfg
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
import models.ray_tracer_ah_3d_ecef as ah
import models.unified_hybrid_pipeline as unified

OUT = None  # The formal runner supplies its task-specific child directory.


def rejected(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError('invalid profile/settings were not rejected')


def provider(values):
    grid = np.broadcast_to(np.asarray(values, dtype=float), (2, 2, 3)).copy()
    return GeocentricGridNeProvider3D(grid, [20., 40.], [110., 130.],
                                      [80., 130., 180.], 6371.,
                                      interpolation='radial_hermite')


def model(ne, freq=12., mode='O', gyro=0., profile=None):
    return ah.RefractiveIndexAH3DECEF(
        ne, earth_radius_km=6371., freq_MHz=freq, wave_mode=mode,
        gyro_MHz=gyro, b_ecef_hat=(0., 0., 1.), integration_profile=profile)


def request(ne):
    return unified.UnifiedPropagationRequest(
        scenario='l2_regression', provider=ne,
        coordinate_model=ah.COORDINATE_MODEL, earth_radius_km=6371.,
        tx_geocentric=(30., 120., 0.), rx_geocentric=(30., 120., 0.),
        x_km=np.array([0., 1.]), z_km=np.array([80., 130.]),
        rt_config={**cfg.RT_3D, 'beta_min': 35., 'beta_max': 35., 'n_fan': 1},
        p2p_config={**cfg.P2P_3D, 'candidate_nearest': 1}, make_pd=False)


def verify_cli_ast():
    tree = ast.parse((ROOT/'main2.py').read_text(encoding='utf-8'))
    found = {}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == 'add_argument' and node.args
                and isinstance(node.args[0], ast.Constant)
                and node.args[0].value in ('--p2p-profile', '--p2p-final-profile')):
            found[node.args[0].value] = {k.arg: ast.literal_eval(k.value)
                                        for k in node.keywords if k.arg in ('choices', 'default', 'help')}
    assert set(found) == {'--p2p-profile', '--p2p-final-profile'}
    for entry in found.values():
        assert entry['choices'] == ('L2',) and entry['default'] is None and 'L2' in entry['help']
    return found


def route_checks(ne):
    req = request(ne)
    assert cfg.UPDATE2['AH3D_INTEGRATION_PROFILE'] == 'L2'
    assert cfg.P2P_3D['solve_profile'] == cfg.P2P_3D['final_retrace_profile'] == 'L2'
    assert cfg.UPDATE2['AH3D_INTEGRATION']['L2'] == cfg.AH3D_L2_SETTINGS
    direct = model(ne)
    factory = unified._ah_model(req, 'O')
    assert direct.integration_profile == factory.integration_profile == 'L2'
    assert direct.spatial_fd_base_km == factory.spatial_fd_base_km == 0.03125
    assert direct.plasma_coefficient == factory.plasma_coefficient == cfg.AH3D_PLASMA_COEFFICIENT_SI
    meta = ne.get_metadata()
    with patch.object(ne, 'get_metadata', return_value={**meta, 'components': [
            {'enabled': True, 'minimum_structure_scale_km': .6}]}):
        assert math.isclose(unified._ah_model(req, 'O').spatial_fd_base_km, .01)

    for old in ('smoke', 'research', 'reference', 'wgs84_reference', 'unknown', ''):
        rejected(lambda old=old: unified._ah_model(req, 'O', old))
        for key in ('solve_profile', 'final_retrace_profile'):
            bad = request(ne)
            bad.p2p_config[key] = old
            rejected(lambda bad=bad: unified.run_unified_pipeline(bad))
    with patch.dict(cfg.UPDATE2, {'AH3D_INTEGRATION_PROFILE': 'smoke'}):
        rejected(lambda: unified.run_unified_pipeline(request(ne)))
    with patch.dict(cfg.UPDATE2['AH3D_INTEGRATION']['L2'], {'rtol': 1e-5}):
        rejected(lambda: unified.run_unified_pipeline(request(ne)))
    with patch.dict(cfg.UPDATE2, {'AH3D_SPATIAL_FD_BASE_KM': .5}):
        rejected(lambda: unified.run_unified_pipeline(request(ne)))

    calls = []
    landing = np.array([-2760., 4780., 3180.])
    def fake_ray(m, **kw):
        calls.append({'leaf': 'fan', 'mode': m.wave_mode, 'profile': kw['profile'],
                      'coefficient': m.plasma_coefficient})
        return {'termination_status': 'ground', 'wave_mode': m.wave_mode,
                'elevation_tx_deg': kw['elevation_deg'], 'bearing_tx_deg': kw['bearing_deg'],
                'trajectory_ecef_km': np.array([landing + 1., landing]),
                'apex_height_km': 100., 'group_path_km': 300., 'tau_ms': 1.}
    def fake_p2p(m, **kw):
        assert kw['seeds']
        calls.append({'leaf': 'p2p', 'mode': m.wave_mode, 'profile': kw['profile'],
                      'coefficient': m.plasma_coefficient})
        return [{'termination_status': 'ground', 'wave_mode': m.wave_mode,
                 'elevation_tx_deg': 35., 'bearing_tx_deg': 0.,
                 'trajectory_ecef_km': np.array([landing + 1., landing]),
                 'apex_height_km': 100., 'group_path_km': 300., 'tau_ms': 1.}]
    with ExitStack() as stack:
        stack.enter_context(patch.dict(cfg.UPDATE2, {'RAY_FAN_INCLUDE_X': True}))
        stack.enter_context(patch.object(ne, 'sample_link_plane', return_value=np.zeros((2, 2))))
        stack.enter_context(patch.object(unified, 'trace_ray_ah3d_ecef', fake_ray))
        stack.enter_context(patch.object(unified, 'find_seeded_rays_ah3d_ecef', fake_p2p))
        stack.enter_context(patch.object(unified, '_nearest_ground_seeds', return_value=[{'seed_id': 'S'}]))
        for name in ('_add_power', '_add_scatter_diagnostics', '_apply_p3_energy_closure',
                     '_annotate_mode_causes', '_tag_common_semantics'):
            stack.enter_context(patch.object(unified, name, lambda *args, **kwargs: None))
        for coefficient in (None, 81.0):
            req = request(ne)
            req.plasma_coefficient = coefficient
            before = len(calls)
            result = unified.run_unified_pipeline(req)
            observed = calls[before:]
            assert len(observed) == 4 and {c['leaf'] for c in observed} == {'fan', 'p2p'}
            assert {c['mode'] for c in observed} == {'O', 'X'}
            assert all(c['profile'] == 'L2' for c in observed)
            expected = cfg.AH3D_PLASMA_COEFFICIENT_SI if coefficient is None else coefficient
            assert all(c['coefficient'] == expected for c in observed)
            assert len(result.modes) == 2
            assert all(m['p2p_trace_profile'] == m['p2p_optimization_trace_profile'] == 'L2'
                       and m['p2p_final_retrace_applied'] is False
                       and m['p2p_final_retrace_reason'] == 'same_profile_solution_reused'
                       for m in result.modes)
    return calls


def settings_checks():
    # Metadata must describe the same default used by the formal entry point.
    no_profile = {k: v for k, v in cfg.UPDATE2.items()
                  if k != 'AH3D_INTEGRATION_PROFILE'}
    with patch.dict(cfg.UPDATE2, no_profile, clear=True):
        assert ('AH3D_PROFILE=' + cfg.AH3D_L2_PROFILE) in unified._feature_gates().split(';')

    # Historical/custom profiles remain explicit diagnostic inputs.
    assert ah._integration_profile('smoke')['atol_k'] == 1e-7
    rejected(lambda: ah._integration_profile('unknown'))
    name = '__l2_invalid_test__'
    canonical = dict(cfg.AH3D_L2_SETTINGS)
    invalid = [dict(canonical, atol_p=x) for x in (0., -1., float('nan'), float('inf'))]
    invalid += [dict(canonical, atol_k=1e-12),
                {k: v for k, v in canonical.items() if k != 'atol_p'},
                dict(canonical, ds_min_km=1., ds_max_km=.5)]
    try:
        for item in invalid:
            cfg.UPDATE2['AH3D_INTEGRATION'][name] = item
            rejected(lambda: ah._integration_profile(name))
    finally:
        cfg.UPDATE2['AH3D_INTEGRATION'].pop(name, None)


def real_short_traces():
    vacuum = provider([0., 0., 0.])
    smooth = provider([8e10, 1e11, 1.2e11])
    launch = dict(tx_lat_deg=30., tx_lon_deg=120., tx_height_km=120.,
                  elevation_deg=30., bearing_deg=25., path_length_km=.25)
    observed = []
    original = ah.DOP853
    def spy(*args, **kwargs):
        observed.append({'rtol': float(kwargs['rtol']), 'atol': np.array(kwargs['atol']),
                         'first_step': float(kwargs['first_step']),
                         'max_step': float(kwargs['max_step'])})
        return original(*args, **kwargs)
    def trace(m, profile=None):
        ray = ah.integrate_fixed_path_ah3d_ecef(m, **launch, profile=profile)
        consumed = observed[-1]
        assert ray['termination_status'] == 'fixed_geometric_path_completed'
        assert consumed['rtol'] == 1e-10 and consumed['max_step'] == .5
        assert consumed['first_step'] == .002
        np.testing.assert_array_equal(consumed['atol'][:3], np.full(3, 1e-8))
        np.testing.assert_array_equal(consumed['atol'][3:6], np.full(3, m.k0*1e-12))
        assert consumed['atol'][6] == 1e-12
        assert math.isclose(ray['integration_settings']['atol_k']/m.k0, 1e-12, rel_tol=1e-15)
        return ray
    rows = []
    frozen = '__l2_frozen_diagnostic__'
    try:
        with patch.object(ah, 'DOP853', side_effect=spy):
            for freq in (8., 12.):
                for mode in ('O', 'X'):
                    m = model(vacuum, freq, mode)
                    ray = trace(m)
                    expected = (np.asarray(ray['trajectory_ecef_km'][0]) +
                                m.launch_direction(30., 120., 30., 25.)*.25)
                    pos_error = float(np.linalg.norm(ray['trajectory_ecef_km'][-1]-expected))
                    group_error = abs(float(ray['trajectory_group_path_km'][-1])-.25)
                    assert pos_error < 1e-8 and group_error < 1e-8
                    rows.append({'case': 'vacuum', 'freq_MHz': freq, 'mode': mode,
                                 'position_error_km': pos_error, 'group_error_km': group_error,
                                 'actual_atol_k': ray['integration_settings']['atol_k']})
            for mode in ('O', 'X'):
                primary = model(smooth, 12., mode, gyro=1.)
                cfg.UPDATE2['AH3D_INTEGRATION'][frozen] = {
                    'rtol': 1e-10, 'atol_r_km': 1e-8, 'atol_k': primary.k0*1e-12,
                    'ds_min_km': .002, 'ds_max_km': .5}
                nominal = trace(primary)
                explicit = trace(model(smooth, 12., mode, gyro=1., profile=frozen), frozen)
                for key in ('trajectory_ecef_km', 'wavevector_ecef_km^-1', 'trajectory_group_path_km'):
                    np.testing.assert_allclose(nominal[key], explicit[key], rtol=0., atol=1e-12)
                rows.append({'case': 'smooth_magnetized_vs_frozen', 'freq_MHz': 12., 'mode': mode,
                             'max_position_delta_km': float(np.max(np.abs(
                                 nominal['trajectory_ecef_km']-explicit['trajectory_ecef_km']))),
                             'canonical_atol_p': nominal['integration_settings']['atol_p'],
                             'legacy_explicit_atol_k': explicit['integration_settings']['atol_k']})
    finally:
        cfg.UPDATE2['AH3D_INTEGRATION'].pop(frozen, None)
    assert len(observed) == 8
    return rows


def summary_check():
    import main2
    saved = main2._OUT
    try:
        main2._OUT = str(OUT)
        main2._save_p2p_summary('l2_test', SimpleNamespace(modes=[], traceability={}, warnings=[]))
        item = json.loads((Path(OUT)/'p2p_unified_l2_test.json').read_text(encoding='utf-8'))
        assert item['solve_profile'] == item['final_retrace_profile'] == 'L2'
    finally:
        main2._OUT = saved


def run():
    if OUT is None:
        raise RuntimeError('output directory is required')
    routes = route_checks(provider([1e11, 1e11, 1e11]))
    settings_checks()
    cli = verify_cli_ast()
    summary_check()
    real = real_short_traces()
    result = {'status': 'PASS', 'scope': 'formal routing and bounded short traces only',
              'routing_checks': routes, 'cli_ast': cli, 'real_trajectory_count': 8,
              'max_path_km': .25, 'real_short_cases': real,
              'limitations': 'no PH, IRI, IGRF, full P2P or full-scene validation'}
    path = Path(OUT)/'l2-defaults-result.json'
    with path.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print('L2_DEFAULTS_PASS routed_stubs=8 real_short_trajectories=8')
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    OUT = args.output.resolve()
    OUT.mkdir(parents=True, exist_ok=False)
    raise SystemExit(run())
