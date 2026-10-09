"""Production-only dependency guard and small current-Provider checks."""
from __future__ import annotations

import ast
import importlib.abc
import hashlib
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def legacy_module(name):
    return name == "legacy" or name.startswith("legacy.") or name == "code.legacy" or name.startswith("code.legacy.")


def legacy_path(value):
    if not isinstance(value, (str, bytes, os.PathLike)):
        return False
    try:
        return Path(os.fsdecode(value)).resolve().is_relative_to(ROOT / "legacy")
    except (OSError, ValueError):
        return False


def static_violations(paths):
    violations = []
    for path in paths:
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                modules = [node.module or ""]
                modules.extend(f"{node.module}.{alias.name}" for alias in node.names if node.module == "code")
            if any(legacy_module(name) for name in modules):
                violations.append(f"{path}:{node.lineno}: legacy import")
            # Only executable string constants, not docstrings/comments or solver labels.
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                value = node.value
                literal_path = value.replace("\\", "/").lower()
                if "\n" not in value and ("/legacy/" in literal_path or literal_path.startswith("legacy/") or value.startswith("legacy.")):
                    violations.append(f"{path}:{node.lineno}: legacy path/module constant")
    return violations


class _NoLegacy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if legacy_module(fullname):
            raise RuntimeError(f"production boundary blocked import: {fullname}")
        return None


def install_guard(iri_driver=None):
    """Run only in a disposable test process; audit hooks cannot be removed."""
    sys.meta_path.insert(0, _NoLegacy())
    allowed_calls = []
    if iri_driver is not None:
        iri_driver = Path(iri_driver).resolve()
        driver_sha = hashlib.sha256(iri_driver.read_bytes()).hexdigest()

    def allowed_iri(args):
        import subprocess
        executable, command, cwd, _env = args
        if iri_driver is None or (executable and Path(executable).resolve() != iri_driver):
            return False
        if Path(cwd or os.getcwd()).resolve() != ROOT:
            return False
        if isinstance(command, str):
            prefix = subprocess.list2cmdline([str(iri_driver)]) + " "
            suffix = " " + subprocess.list2cmdline([str(iri_driver.parent / "data")])
            if not command.startswith(prefix) or not command.endswith(suffix):
                return False
            numbers = command[len(prefix):-len(suffix)].split()
        else:
            numbers = list(command[1:-1])
            if str(command[0]) != str(iri_driver) or str(command[-1]) != str(iri_driver.parent / "data"):
                return False
        try:
            valid = len(numbers) == 11 and all(math.isfinite(float(value)) for value in numbers)
        except (TypeError, ValueError):
            valid = False
        if not valid or hashlib.sha256(iri_driver.read_bytes()).hexdigest() != driver_sha:
            return False
        allowed_calls.append({"executable": str(iri_driver), "sha256": driver_sha,
                              "numeric_arguments": numbers, "cwd": str(ROOT)})
        return True

    def audit(event, args):
        if event == "subprocess.Popen" and allowed_iri(args):
            return
        if event in ("subprocess.Popen", "os.system", "os.posix_spawn", "os.spawn", "os.exec", "os.startfile", "os.startfile/2"):
            raise RuntimeError(f"production boundary blocked child process: {event}")
        if event == "import" and legacy_module(str(args[0])):
            raise RuntimeError("production boundary blocked legacy import")
        if event == "open" and legacy_path(args[0]):
            raise RuntimeError("production boundary blocked legacy file")
        if event == "exec" and legacy_path(getattr(args[0], "co_filename", "")):
            raise RuntimeError("production boundary blocked legacy execution")
        if event == "compile" and len(args) > 1 and legacy_path(args[1]):
            raise RuntimeError("production boundary blocked legacy compilation")

    sys.addaudithook(audit)
    check_loaded()
    return allowed_calls


def check_loaded():
    bad = [name for name, module in tuple(sys.modules.items())
           if legacy_module(name)
           or legacy_path(getattr(module, "__file__", ""))]
    bad.extend(str(path) for path in sys.path if legacy_path(path))
    if bad:
        raise RuntimeError(f"production boundary loaded legacy: {bad}")
    return True


def run():
    import copy
    import importlib
    import subprocess
    import tempfile
    import numpy as np
    import config as cfg
    import main
    import main2
    from models.ne_provider_geocentric import GeocentricGridNeProvider3D, VacuumGeocentricNeProvider3D
    from models.unified_hybrid_pipeline import reject_archived_physics, run_unified_pipeline, UnifiedPropagationRequest

    checks = {}
    with tempfile.TemporaryDirectory(dir=ROOT / "tests") as folder:
        fixture = Path(folder) / "probe.py"
        fixture.write_text("from legacy.models import example\n", encoding="utf-8")
        checks["static_import_rejected"] = bool(static_violations([fixture]))
        fixture.write_text("import importlib\nimportlib.import_module('legacy.models')\n", encoding="utf-8")
        checks["static_dynamic_module_rejected"] = bool(static_violations([fixture]))
        fixture.write_text("solver = 'iso3d_legacy'\n", encoding="utf-8")
        checks["rejection_label_allowed"] = not static_violations([fixture])

    install_guard()
    probes = {
        "dynamic_import_blocked": lambda: importlib.import_module("legacy.models"),
        "legacy_file_blocked": lambda: open(ROOT / "legacy" / "__init__.py"),
        "path_execution_blocked": lambda: exec(compile("pass", str(ROOT / "legacy" / "probe.py"), "exec")),
        "child_process_blocked": lambda: subprocess.run([sys.executable, "-c", "pass"], check=True),
    }
    for name, action in probes.items():
        try:
            action()
        except RuntimeError as exc:
            checks[name] = "production boundary" in str(exc)
        else:
            checks[name] = False

    sys.path.append(str(ROOT / "legacy"))
    try:
        check_loaded()
    except RuntimeError:
        checks["legacy_search_path_rejected"] = True
    else:
        checks["legacy_search_path_rejected"] = False
    finally:
        sys.path.pop()

    saved_absorption, saved_update = copy.deepcopy(cfg.ABSORPTION), copy.deepcopy(cfg.UPDATE2)
    checks["obsolete_absorption_mode_removed"] = "mode" not in cfg.UPDATE2["ABSORPTION_COMPLEX"]
    try:
        UnifiedPropagationRequest(scenario="boundary", provider=VacuumGeocentricNeProvider3D(cfg.RE_KM),
                                  pe_irregularity_params={})
    except TypeError:
        checks["obsolete_pe_payload_rejected"] = True
    else:
        checks["obsolete_pe_payload_rejected"] = False
    try:
        for label, target, key in (
                ("d_layer", cfg.ABSORPTION, "enable"),
                ("pe_enable", cfg.UPDATE2["PE_COUPLING"], "enable"),
                ("pe_power", cfg.UPDATE2["PE_COUPLING"], "apply_to_power")):
            old = target.get(key)
            target[key] = True
            for entry_name, action in (
                    ("shared", reject_archived_physics),
                    ("main", main._enable_main_reference_physics),
                    ("main2", main2._enable_reference_physics),
                    ("power", lambda: main._add_power_3d([])),
                    ("pipeline", lambda: run_unified_pipeline(UnifiedPropagationRequest(
                        scenario="boundary", provider=VacuumGeocentricNeProvider3D(cfg.RE_KM),
                        make_ray_fan=False, make_p2p=False, make_pd=False))),
                    ("request", lambda: main2._request_for("tid"))):
                before = copy.deepcopy((cfg.ABSORPTION, cfg.UPDATE2))
                try:
                    action()
                except RuntimeError:
                    checks[f"{label}_{entry_name}_rejected"] = True
                else:
                    checks[f"{label}_{entry_name}_rejected"] = False
                checks[f"{label}_{entry_name}_not_overwritten"] = repr(before) == repr((cfg.ABSORPTION, cfg.UPDATE2))
            target[key] = old
    finally:
        cfg.ABSORPTION.clear()
        cfg.ABSORPTION.update(saved_absorption)
        cfg.UPDATE2.clear()
        cfg.UPDATE2.update(saved_update)
    checks["configuration_restored"] = repr((saved_absorption, saved_update)) == repr((cfg.ABSORPTION, cfg.UPDATE2))
    mode = {"Pr_W": 2.0, "tau_ms": 3.0}
    result = main._add_ray_pe_reference_3d([mode])
    checks["default_pe_noop_metadata"] = (result[0] is mode and mode["Pr_W"] == 2.0
        and mode["tau_ms"] == 3.0 and mode["pe_solver_dimensionality"] == "off"
        and mode["pe_apply_to_power"] is False
        and mode["pe_validity"] == "archived_local_2p5d_not_in_production_chain")
    try:
        run_unified_pipeline(UnifiedPropagationRequest(
            scenario="boundary", provider=VacuumGeocentricNeProvider3D(cfg.RE_KM),
            apply_ray_pe_reference=True, make_ray_fan=False, make_p2p=False, make_pd=False))
    except RuntimeError:
        checks["request_pe_explicitly_rejected"] = True
    else:
        checks["request_pe_explicitly_rejected"] = False

    lat, lon, height = np.array([29., 31.]), np.array([119., 122., 124.]), np.array([80., 100., 130., 180.])
    density = 1e10 + 1e8 * lat[:, None, None] + 2e8 * lon[None, :, None] + 3e7 * height[None, None, :]
    provider = GeocentricGridNeProvider3D(density, lat, lon, height, 6371.)
    expected = 1e10 + 1e8 * 30. + 2e8 * 120. + 3e7 * 110.
    checks["provider_axis_order_affine"] = bool(abs(provider.get_Ne_geocentric(30., 120., 110.) / expected - 1) < 1e-14)
    checks["provider_domain_fill"] = bool(provider.get_Ne_geocentric(20., 120., 110.) == 0.)
    for name, values, axis in [("provider_shape_rejected", density.transpose(1, 0, 2), lat),
                               ("provider_axis_rejected", density, lat[::-1]),
                               ("provider_negative_rejected", -density, lat)]:
        try:
            GeocentricGridNeProvider3D(values, axis, lon, height, 6371.)
        except ValueError:
            checks[name] = True
        else:
            checks[name] = False
    checks["loaded_modules_clean"] = check_loaded()
    for name, passed in checks.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {name}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.path.insert(0, str(ROOT))
    raise SystemExit(run())
