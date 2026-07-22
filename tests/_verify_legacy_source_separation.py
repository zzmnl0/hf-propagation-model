"""Static acceptance for production/legacy source separation."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ARCHIVED_MODELS = {
    "hybrid_model.py",
    "hybrid_ray_pe.py",
    "local_frame.py",
    "ne_provider.py",
    "p2p_mode_audit_3d.py",
    "parameterized_ionosphere_3d.py",
    "pe_propagator.py",
    "point_to_point.py",
    "point_to_point_3d.py",
    "point_to_point_ah_3d.py",
    "radar_clutter.py",
    "ray_tracer.py",
    "ray_tracer_3d.py",
    "ray_tracer_ah_3d.py",
    "ray_tracer_ah_3d_v2.py",
    "tube_tracer.py",
    "tube_tracer_3d.py",
}


def main() -> int:
    production_models = ROOT / "models"
    legacy_models = ROOT / "legacy" / "models"
    production_files = {path.name for path in production_models.glob("*.py")}
    legacy_files = {path.name for path in legacy_models.glob("*.py")}

    checks = {
        "all_archived_models_present": ARCHIVED_MODELS <= legacy_files,
        "no_archived_model_in_production": not (ARCHIVED_MODELS & production_files),
        "no_duplicate_model_stems": not (
            (production_files - {"__init__.py"})
            & (legacy_files - {"__init__.py"})
        ),
        "legacy_analysis_moved": (
            (ROOT / "legacy" / "analysis" / "parameter_sweep_3d.py").is_file()
            and not (ROOT / "analysis" / "parameter_sweep_3d.py").exists()
        ),
        "legacy_plot_moved": (
            (ROOT / "legacy" / "viz" / "plot_ne_background.py").is_file()
            and not (ROOT / "viz" / "plot_ne_background.py").exists()
        ),
        "legacy_benchmark_moved": (
            (ROOT / "legacy" / "_bench.py").is_file()
            and not (ROOT / "_bench.py").exists()
        ),
    }

    violations = []
    formal_files = [ROOT / "main.py", ROOT / "main2.py"]
    formal_files.extend(production_models.glob("*.py"))
    for path in formal_files:
        text = path.read_text(encoding="utf-8")
        if "legacy." in text:
            violations.append(str(path.relative_to(ROOT)))
    checks["production_has_no_legacy_import"] = not violations

    failures = [name for name, passed in checks.items() if not passed]
    for name, passed in checks.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {name}")
    if violations:
        print("legacy import violations: " + ", ".join(violations))
    print(f"{len(failures)} fail(s)")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
