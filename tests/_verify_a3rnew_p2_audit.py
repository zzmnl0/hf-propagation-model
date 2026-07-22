"""Control-flow verification for the A3Rnew-P2 landing-map audit."""

from pathlib import Path
import tempfile
import json
import sys
import os

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from legacy.models.p2p_mode_audit_3d import adaptive_landing_audit


class _Provider:
    def get_metadata(self):
        return {"provider_hash": "test-provider", "config_hash": "test-config"}


class _AH:
    wave_mode = "O"
    geomag = {"fH_MHz": 0.0}
    n3d = _Provider()


def main():
    tmp = Path(__file__).resolve().parents[1] / "output2" / (
        "_p2_audit_test_{}".format(os.getpid()))
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        params = {
            "audit_beta_bounds_deg": (2.0, 12.0),
            "audit_az_bounds_deg": (-180.0, 180.0),
            "audit_beta_step_deg": 10.0,
            "audit_az_step_deg": 90.0,
            "audit_max_refinement": 1,
            "audit_target_near_km": 100.0,
            "audit_fold_condition": 1.0e4,
            "audit_checkpoint_dir": tmp,
            "audit_resume": True,
            "audit_progress": False,
        }

        def trace(beta, az):
            # A smooth landing map with a deliberate fold in beta and a
            # periodic azimuth coordinate; this is only a control-flow field.
            a = float(az)
            y = ((a + 180.0) % 360.0) - 180.0
            return {
                "x_land_km": 50.0 + 20.0 * (float(beta) - 7.0) ** 2,
                "y_land_km": y,
                "h_reflect_km": 350.0,
                "tau_ms": 10.0,
                "turning_point_count": 1,
                "branch_status": "ground",
            }

        seeds, diag = adaptive_landing_audit(
            (0.0, 0.0, 0.0), (50.0, 0.0, 0.0), _AH(), 10.0, params, 0.0, trace)
        assert diag["full_domain"] and diag["periodic_azimuth"]
        assert diag["node_count"] >= 12
        assert seeds and any("landing_map_fold" in s["audit_triggers"] for s in seeds)
        metadata = list(Path(tmp).rglob("audit.json"))[0]
        saved = json.loads(metadata.read_text(encoding="utf-8"))
        saved["frozen_hash"] = "changed"
        metadata.write_text(json.dumps(saved), encoding="utf-8")
        try:
            adaptive_landing_audit(
                (0.0, 0.0, 0.0), (50.0, 0.0, 0.0), _AH(), 10.0, params, 0.0, trace)
        except RuntimeError as exc:
            assert "hash mismatch" in str(exc)
        else:
            raise AssertionError("checkpoint mismatch was not rejected")
    finally:
        for item in sorted(tmp.rglob("*"), reverse=True):
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                item.rmdir()
        tmp.rmdir()
    print("A3Rnew-P2 audit control flow: 0 fail(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
