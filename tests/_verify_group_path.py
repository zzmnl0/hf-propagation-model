"""
Stage A verification: geometric path vs cold-plasma group path.

Checks:
  1. 2-D single ray returns geom_path_km and group_path_km = integral(ds/n).
  2. 3-D no-flat az=0 matches 2-D tau for the same launch angle.
  3. tau_ms is always computed from group_path_km.
  4. O/X post-processing writes direct AH group paths and matching tau_ms.

Run:
    conda run -n pytorch_cpu python tests/_verify_group_path.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import config as cfg
from utils import group_delay_ms
from models.ionosphere_model import IonosphereModel
from legacy.models.ray_tracer import RefractiveIndex, trace_single_ray
from legacy.models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from main import _apply_OX_split_3d


def _ok(label, condition):
    print("  {} {}".format(label, "OK" if condition else "FAIL"))
    return 0 if condition else 1


def run():
    n_fail = 0

    print("=" * 55)
    print("  Stage A verification: group path and delay")
    print("=" * 55)

    iono = IonosphereModel()
    Ne, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)

    nm_2d = RefractiveIndex(Ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)
    nm_3d_noflat = RefractiveIndex3D(
        Ne, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ, earth_flat=False)

    beta = 30.0
    ray2d = trace_single_ray(cfg.TX_POS, beta, nm_2d, cfg.FREQ_MHZ)
    ray3d = trace_single_ray_3d(
        cfg.TX_POS, beta, 0.0, nm_3d_noflat, cfg.FREQ_MHZ)

    print("\n[1] 2-D group path uses integral(ds/n)")
    print("  geom={:.3f} km  group={:.3f} km  tau={:.6f} ms".format(
        ray2d["geom_path_km"], ray2d["group_path_km"], ray2d["tau_ms"]))
    n_fail += _ok("group_path >= geom_path:",
                  ray2d["group_path_km"] >= ray2d["geom_path_km"])
    n_fail += _ok("tau == group_path/c:",
                  abs(ray2d["tau_ms"] -
                      group_delay_ms(ray2d["group_path_km"])) < 1e-9)

    print("\n[2] 3-D no-flat az=0 matches 2-D single ray")
    diff_tau = abs(ray3d["tau_ms"] - ray2d["tau_ms"])
    diff_gp = abs(ray3d["group_path_km"] - ray2d["group_path_km"])
    print("  tau2d={:.6f} ms  tau3d={:.6f} ms  diff={:.6f} ms".format(
        ray2d["tau_ms"], ray3d["tau_ms"], diff_tau))
    print("  gp2d={:.3f} km  gp3d={:.3f} km  diff={:.3f} km".format(
        ray2d["group_path_km"], ray3d["group_path_km"], diff_gp))
    n_fail += _ok("tau diff < 0.1 ms:", diff_tau < 0.1)
    n_fail += _ok("group-path diff < 30 km:", diff_gp < 30.0)
    n_fail += _ok("3-D tau == group_path/c:",
                  abs(ray3d["tau_ms"] -
                      group_delay_ms(ray3d["group_path_km"])) < 1e-9)

    print("\n[3] O/X direct AH group paths")
    modes_ox = _apply_OX_split_3d([ray3d], Ne, nm_3d_noflat)
    n_fail += _ok("O/X modes produced:", len(modes_ox) >= 1)
    for mode in modes_ox:
        tau_from_gp = group_delay_ms(mode["group_path_km"])
        ok_tau = abs(mode["tau_ms"] - tau_from_gp) < 1e-9
        ok_gp = np.isfinite(mode["group_path_km"]) and mode["group_path_km"] > 0.0
        print("  {} gp={:.3f} km tau={:.6f} ms ratio={:.6f} {}".format(
            mode.get("label", "?"), mode["group_path_km"], mode["tau_ms"],
            mode.get("OX_gp_ratio", float("nan")),
            "OK" if (ok_tau and ok_gp) else "FAIL"))
        if not (ok_tau and ok_gp):
            n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
