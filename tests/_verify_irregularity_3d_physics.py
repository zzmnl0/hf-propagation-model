"""Verify P1 3-D ionospheric component physics semantics."""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
from models.irregularity_detector import diagnose_irregularity_field
from models.ne_provider import ExternalNeProvider3D, GridNeProvider, ParameterizedIonosphere3DProvider
from models.parameterized_ionosphere_3d import ParameterizedIonosphere3DFactory


def _ok(label, cond):
    print(" {} {}".format(label, "OK" if cond else "FAIL"))
    return 0 if cond else 1


def _base():
    x = np.linspace(0.0, 400.0, 33)
    z = np.linspace(90.0, 380.0, 30)
    zz = z[None, :]
    ne = 1.0e10 + 1.2e12 * np.exp(-0.5 * ((zz - 280.0) / 70.0) ** 2)
    ne = np.repeat(ne, x.size, axis=0)
    return GridNeProvider(ne, x, z)


def _check_components_and_factory():
    base = _base()
    factory = ParameterizedIonosphere3DFactory(base, geomag=cfg.GEOMAG, seed=11)
    provider, meta = factory.build([
        {
            "type": "tid",
            "amplitude_list": [0.05],
            "az_deg_list": [25.0],
            "period_s_list": [1200.0],
            "lambda_h_km_list": [200.0],
            "z0_km": 260.0,
            "sigma_z_km": 50.0,
        },
        {
            "type": "bubble",
            "delta_max": 0.5,
            "x0_km": 180.0,
            "y0_km": 0.0,
            "z0_km": 280.0,
            "sigma_perp_km": 35.0,
            "sigma_para_km": 500.0,
        },
    ])
    vals = provider.sample_grid3d([120.0, 180.0], [-20.0, 0.0], [240.0, 280.0])
    comps = meta.get("components", [])
    ok = vals.shape == (2, 2, 2) and np.all(vals >= 0.0)
    ok = ok and len(comps) == 2 and any(c["component_type"] == "TID3D" for c in comps)
    ok = ok and any(c["component_type"] == "field_aligned_bubble" for c in comps)
    return _ok("component_factory_metadata", ok)


def _check_tid_phase_motion():
    base = _base()
    provider = ParameterizedIonosphere3DProvider(
        base,
        tid_params={
            "enable": True,
            "n_components": 1,
            "amplitude_list": [0.1],
            "az_deg_list": [0.0],
            "period_s_list": [1000.0],
            "lambda_h_km_list": [100.0],
            "z0_km": 260.0,
            "sigma_z_km": 80.0,
        },
        geomag=cfg.GEOMAG,
    )
    a = provider.get_Ne((0.0, 0.0, 260.0), t=0.0)
    b = provider.get_Ne((50.0, 0.0, 260.0), t=500.0)
    ok = abs(a - b) / max(a, 1.0) < 0.03
    return _ok("tid_phase_speed_reference", ok)


def _check_bubble_field_aligned():
    base = _base()
    p = {
        "enable": True,
        "delta_max": 0.6,
        "x0_km": 200.0,
        "y0_km": 0.0,
        "z0_km": 280.0,
        "sigma_perp_km": 30.0,
        "sigma_para_km": 600.0,
    }
    provider = ParameterizedIonosphere3DProvider(base, bubble_params=p, geomag=cfg.GEOMAG)
    b_hat = np.asarray(provider.get_metadata()["components"][0]["b_hat"])
    center = np.array([200.0, 0.0, 280.0])
    along = center + b_hat * 250.0
    across = center + np.array([0.0, 250.0, 0.0])
    bg_center = base.get_Ne(center)
    dep_along = 1.0 - provider.get_Ne(along) / max(base.get_Ne(along), 1.0)
    dep_across = 1.0 - provider.get_Ne(across) / max(base.get_Ne(across), 1.0)
    dep_center = 1.0 - provider.get_Ne(center) / max(bg_center, 1.0)
    ok = dep_center > 0.4 and dep_along > dep_across
    return _ok("bubble_field_aligned_extension", ok)


def _check_spreadf_seed_and_spectrum():
    base = _base()
    p = {
        "enable": True,
        "Cs": 0.03,
        "p": 3.0,
        "L0_km": 80.0,
        "l0_m": 1000.0,
        "h_screen_km": 280.0,
        "sigma_z_km": 70.0,
        "anisotropy_ratio": 20.0,
        "n_modes": 128,
    }
    a = ParameterizedIonosphere3DProvider(base, spreadf_params=p, geomag=cfg.GEOMAG, seed=21)
    b = ParameterizedIonosphere3DProvider(base, spreadf_params=p, geomag=cfg.GEOMAG, seed=21)
    x = np.linspace(0.0, 400.0, 64)
    vals_a = np.asarray([a.get_Ne((xx, 0.0, 280.0)) for xx in x])
    vals_b = np.asarray([b.get_Ne((xx, 0.0, 280.0)) for xx in x])
    diag = diagnose_irregularity_field(a, np.linspace(0.0, 400.0, 24), [-20.0, 0.0, 20.0], [240.0, 280.0, 320.0])
    ok = np.allclose(vals_a, vals_b) and diag["spectral_slope_x"] > 0.5
    return _ok("spreadf_seed_spectral_diagnostic", ok)


def _check_es_and_external_npz():
    base = _base()
    es = ParameterizedIonosphere3DProvider(
        base,
        es_params={**cfg.ES_BG, **cfg.ES, "enable": True, "foEs_MHz": 5.0},
        geomag=cfg.GEOMAG,
        seed=5,
    )
    near = es.get_Ne((100.0, 0.0, cfg.ES["h_Es_km"]))
    far = es.get_Ne((100.0, 0.0, cfg.ES["h_Es_km"] + 5.0))
    meta = es.get_metadata()["components"][0]
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "ne.npz")
        x = np.array([0.0, 10.0])
        y = np.array([-5.0, 5.0])
        z = np.array([100.0, 110.0])
        ne_cm3 = np.ones((2, 2, 2)) * 1.0e5
        np.savez(path, Ne=ne_cm3, x_km=x, y_km=y, z_km=z, units="cm^-3")
        ext = ExternalNeProvider3D(path)
        val = ext.get_Ne((5.0, 0.0, 105.0))
    ok = near > far and meta["component_type"] == "Es3D" and abs(val - 1.0e11) < 1.0e6
    return _ok("es_external_npz_semantics", ok)


def run():
    print("=" * 55)
    print(" P1 3-D irregularity physics verification")
    print("=" * 55)
    n_fail = 0
    n_fail += _check_components_and_factory()
    n_fail += _check_tid_phase_motion()
    n_fail += _check_bubble_field_aligned()
    n_fail += _check_spreadf_seed_and_spectrum()
    n_fail += _check_es_and_external_npz()
    print("irregularity_3d_physics: {} fail(s)".format(n_fail))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(run())
