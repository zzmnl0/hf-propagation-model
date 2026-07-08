"""Stage-1 P0 compact 3-D parameter sweep."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

import config as cfg
from models.ne_provider import GridNeProvider, ParameterizedIonosphere3DProvider
from models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output"
WORD_DIR = ROOT / "word"


def synthetic_background():
    x = np.linspace(-50.0, 900.0, 96)
    z = np.linspace(70.0, 420.0, 71)
    zz = z[None, :]
    ne = 9.0e11 * np.exp(-((zz - 250.0) / 85.0) ** 2)
    return x, z, np.repeat(ne, x.size, axis=0)


def default_cases():
    return [
        {"case_id": "bg_8MHz", "freq_MHz": 8.0, "tid_amp": 0.0, "bubble_delta": 0.0},
        {"case_id": "tid_10MHz", "freq_MHz": 10.0, "tid_amp": 0.08, "bubble_delta": 0.0},
        {"case_id": "bubble_12MHz", "freq_MHz": 12.0, "tid_amp": 0.0, "bubble_delta": 0.35},
    ]


def build_provider(case):
    x, z, ne = synthetic_background()
    base = GridNeProvider(ne, x, z)
    tid = {"enable": float(case.get("tid_amp", 0.0)) > 0.0, "n_components": 1, "amplitude_list": [float(case.get("tid_amp", 0.0))], "az_deg_list": [30.0], "period_s_list": [2400.0], "lambda_h_km_list": [300.0], "z0_km": 250.0, "sigma_z_km": 80.0}
    bubble = {"enable": float(case.get("bubble_delta", 0.0)) > 0.0, "delta_max": float(case.get("bubble_delta", 0.0)), "x0_km": 450.0, "y0_km": 0.0, "z0_km": 260.0, "sigma_perp_km": 45.0, "sigma_para_km": 520.0}
    provider = ParameterizedIonosphere3DProvider(base, tid_params=tid, bubble_params=bubble, geomag={**cfg.GEOMAG, "link_bearing_deg": cfg.LINK_BEARING_DEG}, seed=11)
    return x, z, provider


def run_sweep(cases=None, output_dir: Path | None = None):
    out_dir = OUT_DIR if output_dir is None else Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    rt = {**cfg.RT_3D, "ds_km": 2.0, "ds_min_km": 0.5, "ds_max_km": 3.0, "max_steps": 900, "z_stop_km": 420.0}
    for case in (default_cases() if cases is None else cases):
        x, z, provider = build_provider(case)
        ne2d = provider.get_Ne_grid(x, z)
        n3d = RefractiveIndex3D(ne2d, x, z, float(case["freq_MHz"]), earth_flat=False, ne_provider=provider)
        ray = trace_single_ray_3d((0.0, 0.0), 28.0, 0.0, n3d, float(case["freq_MHz"]), rt)
        summary = provider.perturbation_summary([0.0, 450.0, 850.0], [-40.0, 0.0, 40.0], [180.0, 250.0, 320.0])
        rows.append({"case_id": case["case_id"], "freq_MHz": float(case["freq_MHz"]), "tid_amp": float(case.get("tid_amp", 0.0)), "bubble_delta": float(case.get("bubble_delta", 0.0)), "x_land_km": ray["x_land_km"], "y_land_km": ray["y_land_km"], "tau_ms": ray["tau_ms"], "h_reflect_km": ray["h_reflect_km"], "max_relative_perturbation": summary["max_relative_perturbation"]})
    summary_path = out_dir / "stage1_sweep_3d_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    WORD_DIR.mkdir(parents=True, exist_ok=True)
    (WORD_DIR / "stage1_traceability_matrix.md").write_text(build_traceability(), encoding="utf-8")
    return rows, summary_path


def build_traceability():
    return "\n".join([
        "# stage1 P0 traceability matrix",
        "",
        "| ID | Requirement | Code | Verification | Status |",
        "| --- | --- | --- | --- | --- |",
        "| P0-1 | Default feature gates remain off | config.py | tests/_verify_update2_all.py | implemented |",
        "| P0-2 | 3-D NeProvider and parameterized ionosphere | models/ne_provider.py; models/parameterized_ionosphere_3d.py | tests/_verify_ne_provider_3d.py | implemented |",
        "| P0-3 | AH3D Hamiltonian tracing diagnostics | models/ray_tracer_ah_3d.py | tests/_verify_AH3D_physics.py | implemented |",
        "| P0-4 | Reproducible compact 3-D parameter sweep | analysis/parameter_sweep_3d.py | tests/_verify_sweep_3d_reproducibility.py | implemented |",
        "",
    ])


def main():
    rows, path = run_sweep()
    print("wrote {}".format(path))
    print("cases={}".format(len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
