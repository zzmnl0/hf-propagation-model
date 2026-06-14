"""Verify Part 1: IonosphereModel - IRI background, TID, Es, plasma bubble."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from models.ionosphere_model import IonosphereModel

print("=" * 55)
print("  Part 1 verification")
print("=" * 55)

# ── 1. IRI background ─────────────────────────────────────────────────────────
print("\n[1] IRI background profile")
iono_base = IonosphereModel()
Ne_2d, n_2d = iono_base.build_Ne_field(cfg.BG_X, cfg.BG_Z)

mid = len(cfg.BG_X) // 2
Ne_col   = Ne_2d[mid, :]
peak_idx = Ne_col.argmax()
peak_Ne  = Ne_col[peak_idx]
peak_alt = cfg.BG_Z[peak_idx]
foF2     = cfg.K_FP * np.sqrt(peak_Ne) / 1e6

print("  Ne_2d shape : {}  (expected ({}, {}))".format(
    Ne_2d.shape, len(cfg.BG_X), len(cfg.BG_Z)))
print("  F2 peak     : Ne = {:.2e} m^-3   alt = {:.1f} km".format(peak_Ne, peak_alt))
print("  foF2 ~ {:.2f} MHz  (expect ~7-9 MHz)".format(foF2))

assert Ne_2d.shape == (len(cfg.BG_X), len(cfg.BG_Z)), "Shape mismatch"
assert 200 < peak_alt < 400, "F2 peak alt {} km outside 200-400 km".format(peak_alt)
assert 1e11 < peak_Ne < 1e13, "F2 peak Ne {:.2e} outside expected range".format(peak_Ne)
assert np.all(Ne_2d >= 0), "Negative Ne values"
print("  [PASS] shape, F2 peak altitude/density, non-negative")

# ── 2. Horizontal homogeneity ─────────────────────────────────────────────────
print("\n[2] Horizontal homogeneity (IRI only, no TID/bubble)")
col_diff = np.max(np.abs(Ne_2d - Ne_2d[0:1, :]))
print("  Max column difference: {:.2e} m^-3  (expect 0)".format(col_diff))
assert col_diff < 1.0, "IRI baseline should be horizontally uniform"
print("  [PASS] all columns identical")

# ── 3. TID perturbation ───────────────────────────────────────────────────────
print("\n[3] TID perturbation (enable=True, amplitude=0.10)")
iono_tid = IonosphereModel(tid_params={**cfg.TID, 'enable': True})
Ne_tid, _ = iono_tid.build_Ne_field(cfg.BG_X, cfg.BG_Z, t=0.0)

Ne_bg_col  = Ne_2d[mid, :]
Ne_tid_col = Ne_tid[mid, :]

# Use same 1%-of-peak threshold that _add_tid uses for normalisation.
# E-layer altitudes (below the threshold) can have >10% relative perturbation
# in absolute terms while still satisfying the F-region 10% target.
Ne_thresh = Ne_bg_col.max() * 0.01
mask = Ne_bg_col > Ne_thresh                       # matches f_mask in _add_tid
dNe_rel = np.abs(Ne_tid_col[mask] - Ne_bg_col[mask]) / Ne_bg_col[mask]
max_dNe_rel = dNe_rel.max()
print("  Max |dNe/Ne0| in F+E (>1% peak): {:.4f}  (expect <= 0.10)".format(max_dNe_rel))

xvar = Ne_tid[:, peak_idx].std() / Ne_tid[:, peak_idx].mean()
print("  x-variation std/mean @ F-peak alt: {:.4f}  (expect > 0.01)".format(xvar))

assert np.all(Ne_tid >= 0), "Negative Ne in TID field"
assert 0.05 < max_dNe_rel <= 0.10 + 1e-9, \
    "TID amplitude {:.4f} outside (0.05, 0.10]".format(max_dNe_rel)
assert xvar > 0.01, "TID should introduce x-variation"
print("  [PASS] amplitude ~10%, x-variation present, non-negative")

# ── 4. Es layer ───────────────────────────────────────────────────────────────
print("\n[4] Es layer (foEs=5 MHz, h=110 km, dh=115 m)")
iono_es = IonosphereModel(es_params={**cfg.ES, 'enable': True})
Ne_es, _ = iono_es.build_Ne_field(cfg.BG_X, cfg.BG_Z)

h_idx      = np.argmin(np.abs(cfg.BG_Z - cfg.ES['h_Es_km']))
Nmax_exp   = (cfg.ES['foEs_MHz'] * 1e6 / cfg.K_FP) ** 2
Ne_bg_h    = Ne_2d[0, h_idx]
Ne_es_h    = Ne_es[0, h_idx]
Es_contrib = Ne_es_h - Ne_bg_h

print("  z = {:.0f} km:  IRI {:.2e}  ->  with Es {:.2e} m^-3".format(
    cfg.BG_Z[h_idx], Ne_bg_h, Ne_es_h))
print("  Es contribution: {:.2e}  (Nmax expected: {:.2e})".format(Es_contrib, Nmax_exp))

es_xvar = Ne_es[:, h_idx].std()
print("  Es x-variation std: {:.2e}  (expect ~0)".format(es_xvar))

assert abs(Es_contrib - Nmax_exp) / Nmax_exp < 0.01, \
    "Es peak {:.2e} != Nmax {:.2e}".format(Es_contrib, Nmax_exp)
assert es_xvar < 1.0, "Es should be uniform in x"
assert np.all(Ne_es >= 0), "Negative Ne in Es field"
if h_idx > 0:
    assert abs(Ne_es[0, h_idx-1] - Ne_2d[0, h_idx-1]) < 1.0, \
        "No Es contribution at h_Es-2km"
if h_idx < len(cfg.BG_Z) - 1:
    assert abs(Ne_es[0, h_idx+1] - Ne_2d[0, h_idx+1]) < 1.0, \
        "No Es contribution at h_Es+2km"
print("  [PASS] Es peak = Nmax, uniform in x, confined to +/-dh window")

# ── 5. Plasma bubble ──────────────────────────────────────────────────────────
print("\n[5] Plasma bubble (delta_max=0.6, x0=600 km, z0=350 km)")
iono_bub = IonosphereModel(bubble_params={**cfg.BUBBLE, 'enable': True})
Ne_bub, _ = iono_bub.build_Ne_field(cfg.BG_X, cfg.BG_Z)

xi = np.argmin(np.abs(cfg.BG_X - cfg.BUBBLE['x0_km']))
zi = np.argmin(np.abs(cfg.BG_Z - cfg.BUBBLE['z0_km']))
Ne_ref    = Ne_2d[xi, zi]
Ne_centre = Ne_bub[xi, zi]
depletion = 1.0 - Ne_centre / Ne_ref if Ne_ref > 0 else 0.0

print("  Bubble centre ({:.0f}, {:.0f}) km: IRI {:.2e} -> {:.2e} m^-3".format(
    cfg.BG_X[xi], cfg.BG_Z[zi], Ne_ref, Ne_centre))
print("  Depletion: {:.3f}  (expect ~{:.2f})".format(depletion, cfg.BUBBLE['delta_max']))

assert np.all(Ne_bub >= 0), "Negative Ne in bubble field"
assert abs(depletion - cfg.BUBBLE['delta_max']) < 0.05, \
    "Bubble depletion {:.3f} != {:.2f}".format(depletion, cfg.BUBBLE['delta_max'])
xi_far = np.argmin(np.abs(cfg.BG_X - 0.0))
assert abs(Ne_bub[xi_far, zi] - Ne_2d[xi_far, zi]) / (Ne_2d[xi_far, zi] + 1) < 0.01, \
    "Bubble should not reach TX position"
print("  [PASS] depletion ~60%, unaffected at TX, non-negative")

# ── 6. Refractive index ───────────────────────────────────────────────────────
print("\n[6] Refractive index n (IRI baseline)")
n_peak = n_2d[mid, peak_idx]
n_low  = n_2d[mid, 0]
print("  n at F2 peak ({:.0f} km): {:.5f}  (expect < 1)".format(peak_alt, n_peak))
print("  n at z = {:.0f} km:        {:.6f}  (expect ~1)".format(cfg.BG_Z[0], n_low))
assert n_peak < 1.0, "n at F2 peak should be < 1"
assert n_low  > 0.999, "n at low altitude should be ~1"
print("  [PASS] refractive index physically correct")

print()
print("=" * 55)
print("  All Part 1 checks PASSED.")
print("=" * 55)
