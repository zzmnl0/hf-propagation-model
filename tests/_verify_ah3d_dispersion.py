"""A1 checks for physical AH roots, cutoff states, and group velocity."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.dispersion_ah import C_KM_S, appleton_hartree_roots, dispersion_D
from legacy.models.ray_tracer_3d import RefractiveIndex3D
from legacy.models.ray_tracer_ah_3d import RefractiveIndexAH3D, group_velocity_ah3d


def check(label, value):
    print(" {} {}".format(label, "OK" if value else "FAIL"))
    return 0 if value else 1


def model(Ne, mode="O", fH=0.0):
    x = np.array([0.0, 10.0])
    z = np.array([100.0, 110.0])
    n3d = RefractiveIndex3D(np.full((2, 2), Ne), x, z, cfg.FREQ_MHZ,
                            earth_flat=False)
    return RefractiveIndexAH3D(n3d, cfg.FREQ_MHZ, mode,
                               {**cfg.GEOMAG, "fH_MHz": fH},
                               use_dispersion_v2=True)


def run():
    failures = 0
    omega = 2.0 * np.pi * cfg.FREQ_MHZ * 1e6
    roots = appleton_hartree_roots(omega, 0.0, [1.0, 0.0, 0.2],
                                   [0.0, 0.0, 1.0], 1.2)
    failures += check("vacuum_roots", all(abs(v - 1.0) < 1e-12 for v in roots.values()))

    vacuum = model(0.0)
    k = np.array([omega / C_KM_S, 0.0, 0.0])
    vg = group_velocity_ah3d(5.0, 0.0, 105.0, k, vacuum)
    failures += check("vacuum_group_velocity", abs(np.linalg.norm(vg) / C_KM_S - 1.0) < 1e-7)

    Ne = 0.35 * (cfg.FREQ_MHZ * 1e6 / cfg.K_FP) ** 2
    iso = model(Ne)
    n = np.sqrt(iso.n2_raw(5.0, 0.0, 105.0, k))
    k_iso = np.array([n * omega / C_KM_S, 0.0, 0.0])
    vg_iso = group_velocity_ah3d(5.0, 0.0, 105.0, k_iso, iso)
    failures += check("isotropic_group_velocity", abs(np.linalg.norm(vg_iso) / (C_KM_S * n) - 1.0) < 2e-6)

    def D(kv, w):
        root = appleton_hartree_roots(w, Ne, kv, iso.b_hat, iso.fH_MHz)["O"]
        return dispersion_D(kv, w, root).real

    hk, hw = 1e-4, 100.0
    dk_fd = np.array([(D(k_iso + np.eye(3)[i] * hk, omega)
                       - D(k_iso - np.eye(3)[i] * hk, omega)) / (2.0 * hk)
                      for i in range(3)])
    dw_fd = (D(k_iso, omega + hw) - D(k_iso, omega - hw)) / (2.0 * hw)
    vg_fd = -dk_fd / dw_fd
    failures += check("group_velocity_high_precision_fd",
                      np.linalg.norm(vg_iso - vg_fd) / C_KM_S < 2e-5)

    cutoff = (cfg.FREQ_MHZ * 1e6 / cfg.K_FP) ** 2
    near = model(0.995 * cutoff)
    evan = model(1.005 * cutoff)
    failures += check("cutoff_states",
                      near.branch_status(5.0, 0.0, 105.0, k) == "near_cutoff"
                      and evan.branch_status(5.0, 0.0, 105.0, k) == "evanescent")
    failures += check("no_n2_floor", evan.n2(5.0, 0.0, 105.0, k) < 0.0)
    print("A1 dispersion: {} fail(s)".format(failures))
    return failures


if __name__ == "__main__":
    raise SystemExit(1 if run() else 0)
