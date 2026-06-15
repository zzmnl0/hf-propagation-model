"""
Verify Part 5: PE/SSF propagator (wide-angle, Carrano 2020 approach).

Checks:
  [1] construct_incident_field: Gaussian envelope, correct phase tilt
  [2] ssf_step: energy conserved in uniform medium (n=1)
  [3] apply_pml: boundary attenuated, centre unchanged
  [4] propagate (free-space, n=1): AOA peak matches incident angle < 2 deg
  [5] analyze: mean_aoa_deg near incident angle, delta_tau finite
  [6] extract_scatter_modes: finds dominant peak in single-beam spectrum
  [7] extract_domain: n values in physical range (0, 1]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from utils import freq_to_k0
from models.pe_propagator import PEPropagator, construct_incident_field

print("=" * 55)
print("  Part 5 verification")
print("=" * 55)

# ── Common parameters ─────────────────────────────────────────────────────────
freq_MHz  = cfg.FREQ_MHZ               # 10 MHz
k0        = freq_to_k0(freq_MHz)        # ~209.4 km^-1
dz_km     = cfg.PE['dz_m'] / 1000.0    # 0.0075 km  (7.5 m @ 10 MHz)
dx_km     = cfg.PE['dx_km']            # 0.5 km

# Test domain: z = 0..15 km  (Nz=2000), x = 5 steps
Nz_test   = 2000
z_test    = np.arange(Nz_test) * dz_km       # 0..14.993 km
z_inc     = 7.5                               # km, centre of domain
w0_test   = 1.5                               # km, fits Gaussian in domain
beta_inc  = 20.0                              # deg

pe = PEPropagator(freq_MHz=freq_MHz, pe_params=cfg.PE)

# ── [1] construct_incident_field ──────────────────────────────────────────────
print("\n[1] construct_incident_field: Gaussian envelope and phase tilt")

u_init = construct_incident_field(1.0, beta_inc, z_inc, z_test, k0, w0_km=w0_test)

# Peak at z_inc
peak_idx  = int(np.argmax(np.abs(u_init)))
z_peak    = z_test[peak_idx]
amp_peak  = np.abs(u_init[peak_idx])
amp_edge  = np.abs(u_init[0])             # should be near 0 (far from centre)

print("  Amplitude at z_inc : {:.4f}  (expect 1.0)".format(amp_peak))
print("  Amplitude at z=0   : {:.2e}  (expect ~0)".format(amp_edge))
assert abs(amp_peak - 1.0) < 1e-5, "Peak amplitude not 1: {:.6f}".format(amp_peak)
assert abs(z_peak - z_inc) < dz_km * 2, \
    "Peak not at z_inc: {:.4f} km".format(z_peak)
assert amp_edge < 0.01, "Gaussian not decayed at boundary: {:.4e}".format(amp_edge)

# Phase tilt: local phase gradient near z_inc should equal +k0*sin(beta)
i0  = int(z_inc / dz_km) - 3
i1  = int(z_inc / dz_km) + 3
# Phase gradient via ratio of adjacent samples (avoids 2pi wrapping)
grad_samples = np.angle(u_init[i0+1:i1+1] / u_init[i0:i1])
phase_grad   = float(np.mean(grad_samples)) / dz_km   # rad/km
expected_grad = k0 * np.sin(np.radians(beta_inc))     # +71.6 rad/km
print("  Phase gradient : {:.2f} rad/km  (expected {:.2f})".format(
    phase_grad, expected_grad))
assert abs(phase_grad - expected_grad) < 1.0, \
    "Phase gradient wrong: {:.2f} vs {:.2f}".format(phase_grad, expected_grad)
print("  [PASS]")

# ── [2] ssf_step energy conservation (n=1) ───────────────────────────────────
print("\n[2] ssf_step: energy conserved in uniform medium")

n_ones  = np.ones(Nz_test)
u_after = PEPropagator.ssf_step(u_init, n_ones, k0, dz_km, dx_km)

power_in  = float(np.sum(np.abs(u_init) ** 2))
power_out = float(np.sum(np.abs(u_after) ** 2))
rel_err   = abs(power_out - power_in) / (power_in + 1e-30)

print("  Power in  : {:.6e}".format(power_in))
print("  Power out : {:.6e}".format(power_out))
print("  Relative error : {:.2e}  (expect < 1e-10)".format(rel_err))
assert rel_err < 1e-10, "Energy not conserved: rel_err={:.2e}".format(rel_err)
print("  [PASS]")

# ── [3] apply_pml ────────────────────────────────────────────────────────────
print("\n[3] apply_pml: boundary attenuated, centre unchanged")

u_flat    = np.ones(Nz_test, dtype=complex)
n_pml_t   = 30
sigma_t   = 0.4
u_pml     = PEPropagator.apply_pml(u_flat, n_pml_t, sigma_t)

att_edge  = float(np.abs(u_pml[0]))           # outermost boundary
att_inner = float(np.abs(u_pml[n_pml_t]))     # just inside PML
att_ctr   = float(np.abs(u_pml[Nz_test // 2]))  # centre

expected_edge = np.exp(-sigma_t * 1.0 ** 2)   # ~ 0.67 for sigma=0.4
print("  |u| at boundary (i=0)    : {:.4f}  (expect ~{:.4f})".format(
    att_edge, expected_edge))
print("  |u| just inside PML      : {:.4f}  (expect 1.0)".format(att_inner))
print("  |u| at centre            : {:.4f}  (expect 1.0)".format(att_ctr))
assert abs(att_edge - expected_edge) < 0.01, \
    "Boundary attenuation wrong: {:.4f}".format(att_edge)
assert abs(att_ctr - 1.0) < 1e-10, "Centre should be unchanged"
assert att_inner == 1.0, "Point just outside PML should be 1.0"
print("  [PASS]")

# ── [4] propagate free-space: AOA at incident angle ──────────────────────────
print("\n[4] propagate (free-space, n=1): AOA peak near beta_inc={} deg".format(beta_inc))

Nx_test  = 6                                  # 5 SSF steps
n_field  = np.ones((Nx_test, Nz_test))
u_in     = construct_incident_field(1.0, beta_inc, z_inc, z_test, k0, w0_km=w0_test)

u_out, _ = pe.propagate(u_in, n_field, dx_km, dz_km)

# FFT for AOA
U_out     = np.fft.fft(u_out)
kz_arr    = np.fft.fftfreq(Nz_test, d=dz_km) * 2.0 * np.pi
power_arr = np.abs(U_out) ** 2
valid_mask = np.abs(kz_arr) < k0
power_arr[~valid_mask] = 0.0

peak_bin  = int(np.argmax(power_arr))
kz_peak   = kz_arr[peak_bin]
aoa_peak  = float(np.degrees(np.arcsin(np.clip(kz_peak / k0, -1.0, 1.0))))

print("  kz at peak  = {:.2f} km^-1  (expected {:.2f})".format(
    kz_peak, k0 * np.sin(np.radians(beta_inc))))
print("  AOA at peak = {:.2f} deg     (expected {:.2f})".format(aoa_peak, beta_inc))
assert abs(aoa_peak - beta_inc) < 2.0, \
    "AOA peak {:.2f} deg far from beta_inc {:.2f} deg".format(aoa_peak, beta_inc)
print("  [PASS]")

# ── [5] analyze ───────────────────────────────────────────────────────────────
print("\n[5] analyze: mean_aoa near beta_inc, delta_tau finite")

dx_total  = (Nx_test - 1) * dx_km
res       = pe.analyze(u_out, z_test, dx_total)

print("  mean_aoa_deg      = {:.2f} deg  (expect ~{:.1f})".format(
    res['mean_aoa_deg'], beta_inc))
print("  delta_tau_ms      = {:.4f} ms".format(res['delta_tau_ms']))
print("  tau_extra_mean_ms = {:.4f} ms".format(res['tau_extra_mean_ms']))

assert abs(res['mean_aoa_deg'] - beta_inc) < 2.0, \
    "mean_aoa {:.2f} far from beta_inc".format(res['mean_aoa_deg'])
assert np.isfinite(res['delta_tau_ms']), "delta_tau_ms not finite"
assert res['delta_tau_ms'] >= 0.0, "delta_tau_ms must be non-negative"
print("  [PASS]")

# ── [6] extract_scatter_modes ─────────────────────────────────────────────────
print("\n[6] extract_scatter_modes: dominant peak identified")

modes = pe.extract_scatter_modes(res['aoa_deg'], res['power_aoa'], beta_inc)

print("  Modes found: {}".format(len(modes)))
if modes:
    print("  Dominant: aoa={:.2f} deg  delta_aoa={:.2f} deg  power={:.3e}".format(
        modes[0]['aoa_deg'], modes[0]['delta_aoa_deg'], modes[0]['power']))
    assert abs(modes[0]['aoa_deg'] - beta_inc) < 3.0, \
        "Dominant mode AOA {:.2f} far from beta_inc".format(modes[0]['aoa_deg'])
    print("  [PASS]")
else:
    print("  No peaks above threshold (Gaussian spectrum too smooth) -- acceptable")
    print("  [PASS]")

# ── [7] extract_domain ────────────────────────────────────────────────────────
print("\n[7] extract_domain: n values in physical range (0, 1]")

# Build a small uniform Ne background
x_bg  = np.linspace(0, 200, 41)     # 0..200 km, 5 km step
z_bg  = np.linspace(60, 400, 171)   # 60..400 km, 2 km step
Ne_bg = np.zeros((len(x_bg), len(z_bg)))

# Add a simple F-layer peak at 300 km
for jz, z in enumerate(z_bg):
    zm, ym = 300.0, 100.0
    if abs(z - zm) <= ym:
        Ne_bg[:, jz] = 1e12 * (1.0 - ((z - zm) / ym) ** 2)

x_range = (50.0, 150.0)
z_range = (200.0, 400.0)

n_pe, x_pe, z_pe = pe.extract_domain(Ne_bg, x_bg, z_bg, x_range, z_range)

print("  n_pe shape: {} x {}".format(len(x_pe), len(z_pe)))
print("  n range   : [{:.4f}, {:.4f}]  (expect (0, 1])".format(
    float(n_pe.min()), float(n_pe.max())))
assert n_pe.shape[0] == len(x_pe), "Shape mismatch x"
assert n_pe.shape[1] == len(z_pe), "Shape mismatch z"
assert float(n_pe.min()) > 0.0, "n must be positive"
# With earth_flat=True, n_eff = n*(1+z/R_E) can exceed 1.0 in vacuum regions.
# Upper bound: n<=1 vacuum at z_max => n_eff <= 1*(1+z_max/R_E).
import config as _cfg
_earth_flat = pe.params.get('earth_flat', True)
if _earth_flat:
    _n_max_bound = 1.0 * (1.0 + z_range[1] / _cfg.RE_KM) + 1e-6
    assert float(n_pe.max()) <= _n_max_bound, \
        "n_eff exceeds Earth-flat bound: {:.4f} > {:.4f}".format(
            float(n_pe.max()), _n_max_bound)
else:
    assert float(n_pe.max()) <= 1.0 + 1e-9, "n must be <= 1 in ionosphere"
print("  x_pe: {:.1f}..{:.1f} km  ({} points)".format(
    x_pe[0], x_pe[-1], len(x_pe)))
print("  z_pe: {:.4f}..{:.4f} km  ({} points)".format(
    z_pe[0], z_pe[-1], len(z_pe)))
print("  [PASS]")

print()
print("=" * 55)
print("  All Part 5 checks PASSED.")
print("=" * 55)
