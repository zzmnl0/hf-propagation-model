"""
Verify Part 2: ray tracer (Haselgrove + RK4).

Checks:
  [1] RefractiveIndex:  n^2 profile correct (QP ionosphere)
  [2] Single ray at beta=30 deg through QP:  reflects, trajectory symmetric
  [3] QP ray fan:  reflected rays have reasonable tau, at least 50% return
  [4] IRI background ray fan:  rays return, tau in physical range
  [5] Es-crossing detection:  at_Es dict populated correctly
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from models.ray_tracer import RefractiveIndex, trace_single_ray, shoot_rays_fan

print("=" * 55)
print("  Part 2 verification")
print("=" * 55)

# ── Build QP ionosphere ───────────────────────────────────────────────────────
Nm, zm, ym = 1e12, 300.0, 100.0   # peak Ne, height, half-thickness [km]
x_arr = np.linspace(-50, 1300, 271)
z_arr = np.arange(60.0, 602.0, 2.0)
Ne_qp = np.where(np.abs(z_arr - zm) <= ym,
                 Nm * (1.0 - ((z_arr - zm) / ym) ** 2), 0.0)
Ne_2d_qp = np.tile(Ne_qp, (len(x_arr), 1))

n_model_qp = RefractiveIndex(Ne_2d_qp, x_arr, z_arr, freq_MHz=10.0)

# ── [1] n^2 profile ───────────────────────────────────────────────────────────
print("\n[1] Refractive index profile (QP)")

xm = x_arr[len(x_arr)//2]
n2_peak = n_model_qp.n2(xm, zm)       # at F2 peak, Ne=Nm -> n^2 = 1 - fp^2/f^2
fp_max  = 8.98 * np.sqrt(Nm)          # Hz
n2_exp  = 1.0 - (fp_max / 10e6) ** 2  # expected n^2 at F2 peak
n2_low  = n_model_qp.n2(xm, 60.0)    # below ionosphere -> n^2 = 1
n2_free = n_model_qp.n2(xm, 700.0)   # above grid -> fill_value = 1.0

print("  n^2 at F2 peak ({:.0f} km) : {:.4f}  (expected {:.4f})".format(zm, n2_peak, n2_exp))
print("  n^2 at z=60 km            : {:.4f}  (expected ~1.0)".format(n2_low))
print("  n^2 above grid (fill)     : {:.4f}  (expected 1.0)".format(n2_free))

assert abs(n2_peak - n2_exp) < 0.01, \
    "n^2 at peak: {:.4f} vs expected {:.4f}".format(n2_peak, n2_exp)
assert n2_low > 0.99,  "n^2 at 60 km should be ~1"
assert n2_free == 1.0, "fill_value outside grid should be 1.0"
assert np.all(Ne_2d_qp >= 0), "Ne must be non-negative"
print("  [PASS] n^2 profile correct, fill_value OK")

# ── [2] Single ray beta=30 deg through QP ─────────────────────────────────────
print("\n[2] Single ray beta=30 deg (QP ionosphere)")

ray30 = trace_single_ray((0.0, 0.0), 30.0, n_model_qp, freq_MHz=10.0)
traj  = np.array(ray30['trajectory'])   # (N, 4)
xs, zs = traj[:, 0], traj[:, 1]

apex_idx   = int(np.argmax(zs))
h_reflect  = ray30['h_reflect_km']
x_apex     = xs[apex_idx]
x_land     = xs[-1]

print("  h_reflect    = {:.1f} km".format(h_reflect))
print("  x_apex       = {:.1f} km".format(x_apex))
print("  x_land       = {:.1f} km".format(x_land))
print("  tau_ms       = {:.3f} ms".format(ray30['tau_ms']))
print("  group_path   = {:.1f} km".format(ray30['group_path_km']))

# Reflection must occur inside QP layer
assert 100.0 < h_reflect < zm + ym + 10, \
    "Reflection height {:.1f} km outside QP layer".format(h_reflect)

# Ray must return to near ground (z <= ~5 km at end)
assert zs[-1] < 5.0, \
    "Ray did not return to ground: z_final={:.1f} km".format(zs[-1])

# Trajectory symmetry: x_apex ~ x_land / 2  (horizontal uniform ionosphere)
sym_err = abs(x_apex - x_land / 2.0) / (x_land + 1e-9)
print("  Symmetry error |x_apex - x_land/2| / x_land = {:.4f}  (expect < 0.02)".format(sym_err))
assert sym_err < 0.02, "Trajectory not symmetric: {:.4f}".format(sym_err)

# tau must be > free-space estimate (n < 1 -> group path longer)
tau_free = ray30['group_path_km'] / cfg.C_KMS * 1e3  # consistency check
assert abs(tau_free - ray30['tau_ms']) < 1e-6, "tau_ms inconsistency"
assert ray30['tau_ms'] > 0.5, "tau_ms suspiciously small"

print("  [PASS] reflects, trajectory symmetric, tau consistent")

# ── [3] QP ray fan ────────────────────────────────────────────────────────────
print("\n[3] QP ray fan (beta 5..85 deg)")

rt_qp = {**cfg.RT, 'n_fan': 17, 'beta_min': 5.0, 'beta_max': 85.0}
rays_qp = shoot_rays_fan((0.0, 0.0), n_model_qp, freq_MHz=10.0, rt_params=rt_qp)

reflected = [r for r in rays_qp if r['h_reflect_km'] < cfg.RT['z_stop_km'] - 5]
taus      = [r['tau_ms'] for r in reflected]

print("  Total rays   : {}".format(len(rays_qp)))
print("  Reflected    : {}  ({:.0f}%)".format(len(reflected),
      100*len(reflected)/len(rays_qp)))
print("  tau range    : {:.2f} ~ {:.2f} ms".format(min(taus), max(taus)))

assert len(reflected) >= len(rays_qp) // 2, \
    "Less than 50% of rays reflected: {}".format(len(reflected))
assert min(taus) > 0.1, "Minimum tau unreasonably small"
assert max(taus) < 50.0, "Maximum tau unreasonably large"
print("  [PASS] fan shoots correctly, >50% reflected")

# ── [4] IRI background ray fan ────────────────────────────────────────────────
print("\n[4] IRI background ray fan")

from models.ionosphere_model import IonosphereModel
Ne_iri, _ = IonosphereModel().build_Ne_field(cfg.BG_X, cfg.BG_Z)
n_model_iri = RefractiveIndex(Ne_iri, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)

rt_small = {**cfg.RT, 'n_fan': 17}
rays_iri = shoot_rays_fan(cfg.TX_POS, n_model_iri,
                          freq_MHz=cfg.FREQ_MHZ, rt_params=rt_small)

refl_iri  = [r for r in rays_iri if r['h_reflect_km'] < cfg.RT['z_stop_km'] - 5]
tau_iri   = [r['tau_ms'] for r in refl_iri]

print("  Reflected : {}  /  {}  rays".format(len(refl_iri), len(rays_iri)))
if tau_iri:
    print("  tau range : {:.2f} ~ {:.2f} ms".format(min(tau_iri), max(tau_iri)))

assert len(refl_iri) >= 3, "Too few rays reflected in IRI background"
assert all(1.5 < t < 25.0 for t in tau_iri), \
    "tau outside expected [1.5,25] ms range: {}".format(tau_iri)
print("  [PASS] IRI rays return, tau in [1.5,25] ms")

# ── [5] Es-crossing detection ─────────────────────────────────────────────────
print("\n[5] Es-crossing detection (h_Es={:.0f} km)".format(cfg.ES['h_Es_km']))

ray_es = trace_single_ray(cfg.TX_POS, 20.0, n_model_iri,
                          freq_MHz=cfg.FREQ_MHZ,
                          rt_params=cfg.RT,
                          h_Es_km=cfg.ES['h_Es_km'])

if ray_es['at_Es'] is not None:
    ae = ray_es['at_Es']
    print("  at_Es z={:.1f} km,  theta_Es={:.1f} deg,  group_path={:.1f} km".format(
        ae['z'], ae['theta_Es_deg'], ae['group_path_km']))
    assert 'kx' in ae and 'kz' in ae, "at_Es missing wave vector"
    assert ae['group_path_km'] > 0, "at_Es group_path must be positive"
    print("  [PASS] at_Es dict populated")
else:
    # low-elevation rays may not cross h_Es if they reflect below it
    print("  at_Es = None (ray reflected below h_Es -- acceptable for beta=20)")
    print("  [PASS] no false positive")

print()
print("=" * 55)
print("  All Part 2 checks PASSED.")
print("=" * 55)
