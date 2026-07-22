"""
P1 verification: default 3-D P2P seeding is independent of display ray fan.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import config as cfg
import legacy.models.point_to_point_3d as p2p3d


class DummyN3D:
    pass


def _fake_shoot(tx, rx, n3d, freq_MHz, beta_init, az_init=0.0, **kwargs):
    beta = float(beta_init)
    return {
        "x_land_km": float(rx[0]),
        "y_land_km": 0.0,
        "z_land_km": 0.0,
        "h_reflect_km": 120.0 + beta,
        "tau_ms": 1.0 + beta / 100.0,
        "group_path_km": 300.0 + beta,
        "beta_deg_tx": beta,
        "az_deg_tx": float(az_init),
        "beta_deg_rx": beta,
        "az_deg_rx": float(az_init),
        "elevation_tx_deg": beta,
        "elevation_rx_deg": beta,
        "trajectory_3d": [np.array([0.0, 0.0, 0.0, 1.0, 0.0, 1.0]),
                          np.array([float(rx[0]), 0.0, 0.0, 1.0, 0.0, -1.0])],
    }


def _fan(values):
    return [{"x_land_km": 100.0 + i, "beta_deg_tx": float(v)} for i, v in enumerate(values)]


def _betas(modes):
    return [round(float(m["beta_deg_tx"]), 6) for m in modes]


def _ok(label, cond):
    print("  {} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


def run():
    old = p2p3d._shoot_3d
    p2p3d._shoot_3d = _fake_shoot
    try:
        base = {**cfg.P2P_3D, "n_init": 4, "seed_source": "internal_scan"}
        m_a = p2p3d.find_all_rays_3d((0.0, 0.0), (100.0, 0.0), DummyN3D(), cfg.FREQ_MHZ,
                                      p2p_params=base, rays_fan=_fan([10, 20, 30, 40]))
        m_b = p2p3d.find_all_rays_3d((0.0, 0.0), (100.0, 0.0), DummyN3D(), cfg.FREQ_MHZ,
                                      p2p_params=base, rays_fan=_fan([70, 75, 80, 85]))
        warm = {**cfg.P2P_3D, "n_init": 3, "seed_source": "fan_warm_start"}
        m_w = p2p3d.find_all_rays_3d((0.0, 0.0), (100.0, 0.0), DummyN3D(), cfg.FREQ_MHZ,
                                      p2p_params=warm, rays_fan=_fan([70, 75, 80, 85]))
        n_fail = 0
        n_fail += _ok("internal_scan ignores rays_fan:", _betas(m_a) == _betas(m_b))
        n_fail += _ok("fan_warm_start uses fan seeds:", _betas(m_w) == [70.0, 75.0, 80.0])
        n_fail += _ok("seed source recorded:", all(m.get("p2p_seed_source") == "internal_scan" for m in m_a)
                      and all(m.get("p2p_seed_source") == "fan_warm_start" for m in m_w))
        n_fail += _ok("cause fields recorded:", all("cause_layer" in m and "cause_launch" in m
                      and "cause_mechanism" in m and "cause_summary" in m for m in m_a + m_w))
    finally:
        p2p3d._shoot_3d = old
    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == "__main__":
    sys.exit(run())
