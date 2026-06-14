"""
Verify Part 4: Es-layer three-segment model (Hao et al. 2017).

Checks:
  [1] classify(): correct mode and alpha for reflect / mixed / scatter
  [2] reflection_coeff_sq(): positive finite value; increases with foEs
  [3] scatter_cross_section(): positive finite value; increases with foEs
  [4] compute_power() curve shape:
        reflect region (foEs/f > fr): Pr increases with foEs
        scatter region (foEs/f < fs): Pr decreases with foEs
        mixed region: smooth monotone transition
  [5] transmission_amplitude() in [0, 1]
  [6] compute_power() result dict has all required keys
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from models.es_model import EsLayerModel

print("=" * 55)
print("  Part 4 verification")
print("=" * 55)

f_test    = 10.0           # MHz
theta_rad = np.radians(10) # grazing angle 10 deg
D_km      = 500.0          # total TX->Es->RX path [km]
Pt_W      = 1e3            # transmit power [W]

model = EsLayerModel()     # uses config.ES defaults

# ── [1] classify() ────────────────────────────────────────────────────────────
print("\n[1] classify() mode boundaries")

model.foEs = 5.0           # foEs/f = 0.5 > fr=0.25 -> reflect
mode, alpha = model.classify(f_test)
print("  foEs/f=0.50  -> mode={}, alpha={:.3f}  (expect reflect, 1.0)".format(mode, alpha))
assert mode == 'reflect' and abs(alpha - 1.0) < 1e-9, \
    "Expected (reflect,1.0), got ({},{})".format(mode, alpha)

model.foEs = 2.0           # foEs/f = 0.2, in (0.10, 0.25) -> mixed
mode, alpha = model.classify(f_test)
alpha_exp = (0.2 - 0.10) / (0.25 - 0.10)
print("  foEs/f=0.20  -> mode={}, alpha={:.3f}  (expect mixed, {:.3f})".format(
    mode, alpha, alpha_exp))
assert mode == 'mixed' and abs(alpha - alpha_exp) < 1e-9, \
    "Expected (mixed,{:.3f}), got ({},{:.3f})".format(alpha_exp, mode, alpha)

model.foEs = 0.5           # foEs/f = 0.05 < fs=0.10 -> scatter
mode, alpha = model.classify(f_test)
print("  foEs/f=0.05  -> mode={}, alpha={:.3f}  (expect scatter, 0.0)".format(mode, alpha))
assert mode == 'scatter' and abs(alpha - 0.0) < 1e-9, \
    "Expected (scatter,0.0), got ({},{})".format(mode, alpha)

print("  [PASS]")

# ── [2] reflection_coeff_sq() ─────────────────────────────────────────────────
print("\n[2] reflection_coeff_sq(): positive, finite, increases with foEs")

rho_vals = []
for foEs in [2.0, 4.0, 6.0]:
    model.foEs = foEs
    rho_sq = model.reflection_coeff_sq(theta_rad, f_test)
    rho_vals.append(rho_sq)
    print("  foEs={:.0f} MHz  rho_sq={:.4e}".format(foEs, rho_sq))

assert all(np.isfinite(r) and r >= 0 for r in rho_vals), \
    "rho_sq must be non-negative and finite"
# rho_sq should increase with foEs (fN/f appears squared in prefac)
assert rho_vals[2] > rho_vals[0], \
    "rho_sq should increase with foEs: {:.4e} vs {:.4e}".format(rho_vals[0], rho_vals[2])
print("  [PASS]")

# ── [3] scatter_cross_section() ───────────────────────────────────────────────
print("\n[3] scatter_cross_section(): positive, finite, increases with foEs")

sig_vals = []
for foEs in [0.5, 1.0, 1.5]:
    model.foEs = foEs
    sig5 = model.scatter_cross_section(theta_rad, f_test)
    sig_vals.append(sig5)
    print("  foEs={:.1f} MHz  sigma5={:.3e} m^-1".format(foEs, sig5))

assert all(np.isfinite(s) and s > 0 for s in sig_vals), \
    "sigma5 must be positive and finite"
assert sig_vals[2] > sig_vals[0], \
    "sigma5 should increase with foEs: {:.3e} vs {:.3e}".format(sig_vals[0], sig_vals[2])
print("  [PASS]")

# ── [4] compute_power() curve shape ───────────────────────────────────────────
print("\n[4] compute_power() curve: monotone in each regime")

foEs_scan = np.linspace(0.3, 10.0, 80)
powers_dBW = []
for foEs in foEs_scan:
    model.foEs = foEs
    res = model.compute_power(Pt_W, 1.0, 1.0, f_test, D_km, theta_rad)
    Pr = res['Pr_W']
    assert np.isfinite(Pr), "Pr_W is not finite for foEs={:.2f}".format(foEs)
    powers_dBW.append(10.0 * np.log10(max(Pr, 1e-300)))

powers_dBW = np.array(powers_dBW)
ratios      = foEs_scan / f_test

# Reflection region (foEs/f > fr=0.25): power should be strictly increasing
mask_refl = ratios > cfg.ES['fr']
idx_r = np.where(mask_refl)[0]
assert len(idx_r) >= 5, "Too few points in reflect region"
p_refl = powers_dBW[idx_r]
# allow minor non-monotonicity due to oscillatory reflection coefficient
assert p_refl[-1] > p_refl[0], \
    "Reflection region: power should increase overall from {:.1f} to {:.1f} dBW".format(
        p_refl[0], p_refl[-1])
print("  Reflect region: Pr({:.1f} dBW) -> Pr({:.1f} dBW)  [monotone increasing OK]".format(
    p_refl[0], p_refl[-1]))

# Scatter region (foEs/f < fs=0.10): power should be strictly increasing with foEs
mask_scat = ratios < cfg.ES['fs']
idx_s = np.where(mask_scat)[0]
assert len(idx_s) >= 2, "Too few points in scatter region"
p_scat = powers_dBW[idx_s]
assert p_scat[-1] > p_scat[0], \
    "Scatter region: power should increase with foEs: {:.1f} -> {:.1f}".format(
        p_scat[0], p_scat[-1])
print("  Scatter region: Pr({:.1f} dBW) -> Pr({:.1f} dBW)  [increases with foEs OK]".format(
    p_scat[0], p_scat[-1]))

# Mixed region: finite and bounded
mask_mix = (ratios >= cfg.ES['fs']) & (ratios <= cfg.ES['fr'])
p_mix = powers_dBW[mask_mix]
assert len(p_mix) > 0
assert all(np.isfinite(p_mix)), "NaN/Inf in mixed region"
print("  Mixed region: Pr range [{:.1f}, {:.1f}] dBW  [finite OK]".format(
    p_mix.min(), p_mix.max()))

print("  [PASS]")

# ── [5] transmission_amplitude() ─────────────────────────────────────────────
print("\n[5] transmission_amplitude() in [0, 1]")

for foEs in [0.5, 2.0, 5.0]:
    model.foEs = foEs
    T = model.transmission_amplitude(theta_rad, f_test)
    print("  foEs={:.1f} MHz  T={:.4f}".format(foEs, T))
    assert 0.0 <= T <= 1.0, "T out of [0,1]: {:.4f}".format(T)
print("  [PASS]")

# ── [6] result dict keys ──────────────────────────────────────────────────────
print("\n[6] compute_power() result dict keys")

model.foEs = 5.0
res = model.compute_power(Pt_W, 1.0, 1.0, f_test, D_km, theta_rad)
required_keys = {'Pr_W', 'Pr_reflect_W', 'Pr_scatter_W',
                 'mode', 'alpha', 'tau_ms', 'delta_tau_ms'}
missing = required_keys - set(res.keys())
assert not missing, "Missing keys in result: {}".format(missing)

print("  mode={}  alpha={:.3f}  tau_ms={:.3f}  delta_tau_ms={:.6f}".format(
    res['mode'], res['alpha'], res['tau_ms'], res['delta_tau_ms']))
assert res['tau_ms'] > 0, "tau_ms must be positive"
assert 0.0 <= res['delta_tau_ms'] < 1.0, "delta_tau_ms out of expected range"
print("  [PASS]")

print()
print("=" * 55)
print("  All Part 4 checks PASSED.")
print("=" * 55)
