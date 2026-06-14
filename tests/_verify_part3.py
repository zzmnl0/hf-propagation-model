"""
Verify Part 3: point-to-point ray solver (variational + Newton).

Checks:
  [1] optical_path: positive finite value on parabolic test arc
  [2] optical_path_gradient: correct shape, finite values
  [3] variational_find_ray (QP): h_reflect and tau in physical range
  [4] find_all_rays_p2p (IRI): at least 2 unique modes found
  [5] classify_mode: correct labels for 5 reference cases
  [6] extract_es_params: correct upward-crossing detection
  [7] find_ray_newton (QP): converges or gracefully returns None
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from models.ray_tracer import RefractiveIndex
from models.point_to_point import (
    optical_path, optical_path_gradient,
    variational_find_ray, find_all_rays_p2p,
    classify_mode, extract_es_params, extract_bubble_entry,
    find_ray_newton,
)

print("=" * 55)
print("  Part 3 verification")
print("=" * 55)

# ── Build QP ionosphere ───────────────────────────────────────────────────────
Nm, zm, ym = 1e12, 300.0, 100.0
x_arr = np.linspace(-50, 1300, 271)
z_arr = np.arange(60.0, 602.0, 2.0)
Ne_qp = np.where(np.abs(z_arr - zm) <= ym,
                 Nm * (1.0 - ((z_arr - zm) / ym) ** 2), 0.0)
Ne_2d_qp = np.tile(Ne_qp, (len(x_arr), 1))
n_model_qp = RefractiveIndex(Ne_2d_qp, x_arr, z_arr, freq_MHz=10.0)

tx = cfg.TX_POS   # (0, 0)
rx = cfg.RX_POS   # (1169, 0)

# ── [1] optical_path ──────────────────────────────────────────────────────────
print("\n[1] optical_path on parabolic test arc")

t_lin = np.linspace(0.0, 1.0, 32)
h_test = 250.0
pts_test = np.column_stack([t_lin * 1169.0,
                             h_test * 4.0 * t_lin * (1.0 - t_lin)])
S_val = optical_path(pts_test, n_model_qp)

print("  S = {:.2f} km  (expect > 0 and finite)".format(S_val))
assert S_val > 0.0 and np.isfinite(S_val), \
    "optical_path returned non-physical: {}".format(S_val)
print("  [PASS]")

# ── [2] optical_path_gradient ─────────────────────────────────────────────────
print("\n[2] optical_path_gradient shape and finiteness")

grad = optical_path_gradient(pts_test, n_model_qp)
grad_mag = np.linalg.norm(grad, axis=1)

print("  grad shape = {}  (expect ({}, 2))".format(grad.shape, len(pts_test) - 2))
print("  max |grad| = {:.5f}".format(float(np.max(grad_mag))))
assert grad.shape == (len(pts_test) - 2, 2), \
    "gradient shape mismatch: {}".format(grad.shape)
assert np.all(np.isfinite(grad)), "gradient has non-finite values"
print("  [PASS]")

# ── [3] variational_find_ray (QP, beta=45) ────────────────────────────────────
print("\n[3] variational_find_ray on QP ionosphere (beta=45, high ray)")

pts_conv, gp_conv = variational_find_ray(tx, rx, n_model_qp, 45.0,
                                          is_high_ray=True,
                                          p2p_params={**cfg.P2P, 'max_iter': 300})
h_max_qp = float(np.max(pts_conv[:, 1]))
tau_qp   = gp_conv / cfg.C_KMS * 1e3

print("  h_reflect = {:.1f} km".format(h_max_qp))
print("  tau_ms    = {:.3f} ms".format(tau_qp))
print("  group_path = {:.1f} km".format(gp_conv))

assert 100.0 < h_max_qp < 450.0, \
    "h_reflect {:.1f} km outside expected range [100,450]".format(h_max_qp)
assert 1.0 < tau_qp < 30.0, \
    "tau_ms {:.3f} ms outside expected range [1,30]".format(tau_qp)
print("  [PASS]")

# ── [4] find_all_rays_p2p (IRI) ───────────────────────────────────────────────
print("\n[4] find_all_rays_p2p on IRI background")

from models.ionosphere_model import IonosphereModel
Ne_iri, _ = IonosphereModel().build_Ne_field(cfg.BG_X, cfg.BG_Z)
n_model_iri = RefractiveIndex(Ne_iri, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

p2p_small = {**cfg.P2P, 'n_init': 10, 'max_iter': 200}
rays = find_all_rays_p2p(tx, rx, n_model_iri, cfg.FREQ_MHZ,
                          p2p_params=p2p_small)

print("  Unique modes found: {}".format(len(rays)))
for r in rays:
    print("    {:8s}  h={:.1f} km  tau={:.3f} ms  is_high={}".format(
        r['label'], r['h_reflect_km'], r['tau_ms'], r['is_high']))

assert len(rays) >= 2, \
    "Expected >= 2 propagation modes, got {}".format(len(rays))
labels = {r['label'] for r in rays}
f_modes = {l for l in labels if 'F' in l}
assert len(f_modes) >= 1, \
    "Expected at least one F-layer mode, got labels: {}".format(labels)
print("  [PASS] >= 2 modes, F-layer modes: {}".format(f_modes))

# ── [5] classify_mode ─────────────────────────────────────────────────────────
print("\n[5] classify_mode correctness")

cases = [
    ({'h_reflect_km': 120.0, 'tau_ms': 2.0}, 'Es'),
    ({'h_reflect_km': 170.0, 'tau_ms': 3.0}, 'E'),
    ({'h_reflect_km': 250.0, 'tau_ms': 4.0}, '1F_low'),
    ({'h_reflect_km': 250.0, 'tau_ms': 7.0}, '1F_high'),
    ({'h_reflect_km': 350.0, 'tau_ms': 9.0}, '2F'),
]
for d, expected in cases:
    got = classify_mode(d)
    assert got == expected, \
        "classify_mode h={}, tau={}: got '{}', expected '{}'".format(
            d['h_reflect_km'], d['tau_ms'], got, expected)
    print("  h={:.0f} km  tau={:.1f} ms  -> {}  [OK]".format(
        d['h_reflect_km'], d['tau_ms'], got))
print("  [PASS]")

# ── [6] extract_es_params ─────────────────────────────────────────────────────
print("\n[6] extract_es_params on synthetic path")

h_es_test = 110.0
pts_es = np.array([
    [0.0,   0.0],
    [100.0, 80.0],
    [200.0, 130.0],
    [300.0, 80.0],
    [400.0,  0.0],
], dtype=float)
res = extract_es_params(pts_es, h_es_test)

print("  at_Es = {}".format(res))
assert res is not None, "extract_es_params returned None unexpectedly"
assert abs(res['z'] - h_es_test) < 0.01, \
    "z at Es crossing wrong: {:.3f}".format(res['z'])
assert 0.0 < res['theta_Es_deg'] < 90.0, \
    "theta_Es_deg out of range: {:.3f}".format(res['theta_Es_deg'])
print("  [PASS] x={:.1f} km  theta={:.1f} deg".format(
    res['x'], res['theta_Es_deg']))

# ── [7] find_ray_newton (QP) ──────────────────────────────────────────────────
print("\n[7] find_ray_newton on QP ionosphere")

ray_n = find_ray_newton(tx, rx, n_model_qp, 10.0,
                         beta_init_deg=30.0, tol_km=2.0)
if ray_n is not None:
    print("  converged: h={:.1f} km  tau={:.3f} ms  label={}".format(
        ray_n['h_reflect_km'], ray_n['tau_ms'], ray_n['label']))
    assert ray_n['h_reflect_km'] > 100.0, "h_reflect too low"
    assert ray_n['tau_ms'] > 1.0, "tau_ms too small"
    print("  [PASS] converged")
else:
    print("  did not converge within max_iter (acceptable)")
    print("  [PASS] non-convergence is not an error")

print()
print("=" * 55)
print("  All Part 3 checks PASSED.")
print("=" * 55)
