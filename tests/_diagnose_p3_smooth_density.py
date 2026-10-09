"""Read-only native PH diagnostics for the frozen P3 smooth-entry rejection."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline

import _run_a4_c1m_attribution as common


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    directory = run / 'experiments'
    manifest_path = directory / 'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    grid_path = directory / 'grids/smooth_s5.mat'
    grid = common.mat(grid_path)
    height = np.asarray(grid['height_km']).reshape(-1)
    density = np.asarray(grid['ne_profile_m3']).reshape(-1)
    spline = CubicHermiteSpline(height, density, np.gradient(density, height, edge_order=2))
    out = run / 'production_diagnostics/ph_native_smooth_entry'
    out.mkdir(parents=True, exist_ok=False)
    rows, summaries, sources = [], [], []
    for job in manifest['jobs']:
        if job['role'] not in ('R2', 'R4') or job['level'] != 3:
            continue
        path = directory / 'pharlap_raw' / (job['job_id'] + '.mat')
        data = common.mat(path)
        assert data['job_sha256'] == job['job_sha256']
        assert data['grid_sha256'] == job['grid_sha256']
        native = data['ray_path_data']
        h = common.vector(native, 'height')
        ne = common.vector(native, 'electron_density') * 1e6
        n = common.vector(native, 'refractive_index')
        p = np.column_stack([common.vector(data['ray_state_vec'], k) for k in ('dir_x', 'dir_y', 'dir_z')])
        assert len(h) == len(p) == len(ne) == len(n)
        signed = spline(h)
        selected = np.flatnonzero((h >= 18.) & (h <= 24.))
        negative = selected[signed[selected] < 0.]
        scale = manifest['coefficient'] / (job['frequency_MHz'] * 1e6) ** 2
        for i in selected:
            rows.append(dict(job_id=job['job_id'], node=int(i), height_km=float(h[i]),
                             signed_Hermite_Ne_m3=float(signed[i]), PH_output_Ne_m3=float(ne[i]),
                             clipped_Hermite_Ne_m3=float(max(0., signed[i])),
                             output_n2=float(n[i]**2), state_p2=float(p[i] @ p[i]),
                             state_residual_output_Ne=float(p[i] @ p[i] - (1-scale*ne[i])),
                             n2_residual_output_Ne=float(n[i]**2 - (1-scale*ne[i]))))
        summaries.append(dict(job_id=job['job_id'], negative_support_nodes=int(len(negative)),
                              zero_output_at_negative_support=bool(len(negative) and np.all(ne[negative] == 0.)),
                              max_clipped_Ne_difference_m3=float(np.max(np.abs(ne[selected]-np.maximum(0., signed[selected])))),
                              max_state_residual=float(np.max(np.abs(np.sum(p[selected]**2, axis=1)-(1-scale*ne[selected]))))))
        sources.append(dict(path=str(path), sha256=common.sha(path)))
    assert len(summaries) == 6, 'Expected four R2 and two R4 finest PH rays'
    with (out / 'native_points.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = dict(diagnostic_only=True, production_changed=False, cases=summaries,
                  interpretation='Nonnegative reported Ne is evidence about output semantics, not proof of internal gradient clipping or state reset.',
                  sources=sources, grid_sha256=common.sha(grid_path), manifest_sha256=common.sha(manifest_path),
                  source_sha256=common.sha(Path(__file__)))
    (out / 'summary.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    (out / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
