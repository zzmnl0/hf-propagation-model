"""Radial PCHIP boundary/FD engineering diagnostics; no formal P4 decision."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.io import loadmat
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
from models.ray_tracer_ah_3d_ecef import (
    RefractiveIndexAH3DECEF, SpatialDerivativeError,
    launch_state_ah3d_ecef, spatial_gradient_ah3d_ecef,
)


def main(out):
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve()
    shutil.copy2(source, out / source.name)
    grid_path = ROOT / 'output2/a4_pharlap/c1/c1m/p3_20261009_v1/experiments/grids/smooth_s5.mat'
    grid = loadmat(grid_path, simplify_cells=True)
    provider = GeocentricGridNeProvider3D(
        grid['Ne_m3'], grid['latitude_deg'], grid['longitude_deg'],
        grid['height_km'], 6371., interpolation='radial_pchip')
    pp = provider._radial_spline
    boundary_checks = {}
    toy = GeocentricGridNeProvider3D(
        np.broadcast_to([0., 1., 2.], (2, 2, 3)).copy(),
        [29., 31.], [119., 121.], [0., 1., 2.], 6371., interpolation='radial_pchip')
    boundary_checks['zero_value_nonzero_slope_is_not_C1_vacuum'] = bool(
        toy.get_Ne_geocentric(30., 120., -1.) == 0.
        and toy._radial_spline(0.) == 0.
        and toy._radial_spline.derivative()(0.) != 0.)
    for name, lat, lon in (
        ('latitude_low', provider.latitude_deg[0] - .01, 120.),
        ('latitude_high', provider.latitude_deg[-1] + .01, 120.),
        ('longitude_low', 30., provider.longitude_deg[0] - .01),
        ('longitude_high', 30., provider.longitude_deg[-1] + .01),
    ):
        try:
            provider.get_Ne_geocentric(lat, lon, 100.)
        except ValueError:
            boundary_checks[name] = True
        else:
            boundary_checks[name] = False
    heights = (0., .03125, .125, .5, 19.5, 19.875, 20., 20.125,
               20.5, 99.875, 100., 100.125, 101.875, 102., 102.125)
    rows = []
    for height in heights:
        for step in (.25, .125, .0625):
            model = RefractiveIndexAH3DECEF(
                provider, earth_radius_km=6371., freq_MHz=10., gyro_MHz=0.,
                collision_s=0., spatial_fd_base_km=step)
            state = launch_state_ah3d_ecef(
                model, tx_lat_deg=30., tx_lon_deg=120., tx_height_km=height,
                elevation_deg=22., bearing_deg=0.)
            actual_height = float(np.linalg.norm(state[:3]) - model.earth_radius_km)
            derivative = float(pp.derivative()(actual_height)) if actual_height >= 0 else 0.
            exact = (model.k0**2 * model.plasma_coefficient / 1e14
                     * derivative * state[:3] / np.linalg.norm(state[:3]))
            limit = (model.spatial_fd_stability_rtol * np.linalg.norm(exact)
                     + model.k0**2 * model.spatial_fd_stability_atol_n2_per_km)
            row = dict(height_km=height, actual_height_km=actual_height,
                       fd_base_km=step, analytic_normalized_norm=float(np.linalg.norm(exact) / model.k0**2),
                       engineering_normalized_limit=float(limit / model.k0**2),
                       normalized_absolute_error=None, error_to_existing_tolerance=None,
                       platform_stability_ratio=None, selected_steps_km=None,
                       engineering_consistent=False, status='fd_rejected', failure=None)
            try:
                gradient, diagnostics = spatial_gradient_ah3d_ecef(state, model)
            except SpatialDerivativeError as exc:
                row['failure'] = str(exc)
            else:
                error = float(np.linalg.norm(gradient - exact))
                row.update(normalized_absolute_error=error / model.k0**2,
                           error_to_existing_tolerance=error / limit,
                           platform_stability_ratio=diagnostics['stability_ratio'],
                           selected_steps_km=list(diagnostics['selected_steps_km']),
                           engineering_consistent=bool(error <= limit), status='evaluated')
            rows.append(row)
    convergence = []
    for height in heights:
        samples = [r for r in rows if r['height_km'] == height]
        errors = [r['normalized_absolute_error'] for r in samples]
        convergence.append(dict(height_km=height, normalized_errors_coarse_to_fine=errors,
                                nonincreasing_error=(None if any(e is None for e in errors)
                                                     else bool(errors[2] <= errors[1] <= errors[0]))))
    with (out / 'fd_checks.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = dict(boundary_checks=boundary_checks, rows=rows, convergence=convergence,
                  engineering_consistency_passed=bool(all(boundary_checks.values())
                                                     and all(r['engineering_consistent'] for r in rows)),
                  formal_physical_acceptance_passed=False, recommend_p4=False,
                  scope='Boundary rejection and local same-PCHIP FD checks using existing model tolerance; convergence is recorded separately, with no new formal P4 threshold.',
                  provider_metadata=provider.get_metadata(),
                  input_sha256=hashlib.sha256(grid_path.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (out / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(dict(engineering_consistency_passed=result['engineering_consistency_passed'],
                          evaluated=len(rows), inconsistent=sum(not r['engineering_consistent'] for r in rows),
                          formal_physical_acceptance_passed=False)), flush=True)
    if not result['engineering_consistency_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True, type=Path)
    main(parser.parse_args().out.resolve())
