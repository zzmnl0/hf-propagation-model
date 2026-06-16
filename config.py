"""
Shared configuration for the HF shortwave propagation hybrid model.
Test link : TX @ 30°N 120°E, RX @ 1169 km, 10 MHz
Coordinates: 2-D Cartesian  x = horizontal distance from TX [km]
                             z = height above Earth surface  [km]
Earth-flattening approximation is applied (valid for paths < ~1000 km).
"""
from datetime import datetime
import numpy as np

# ── Physical constants ────────────────────────────────────────────────────────
C_MS   = 2.998e8      # Speed of light  [m/s]
C_KMS  = 2.998e5      # Speed of light  [km/s]
RE_KM  = 6371.0       # Earth mean radius [km]
E_C    = 1.602e-19    # Electron charge [C]
M_E    = 9.109e-31    # Electron mass   [kg]
EPS0   = 8.854e-12    # Vacuum permittivity [F/m]
K_FP   = 8.980        # fp [Hz] = K_FP * sqrt(Ne [m⁻³])

# ── Test link ─────────────────────────────────────────────────────────────────
TX_LAT    = 30.0           # Transmitter latitude  [°N]
TX_LON    = 120.0          # Transmitter longitude [°E]
TX_POS    = (0.0, 0.0)     # TX in 2-D grid (x=0, z=0) [km, km]
RX_RANGE  = 1169.0         # TX–RX horizontal distance  [km]
RX_POS    = (1169.0, 0.0)  # RX in 2-D grid [km, km]
FREQ_MHZ  = 10.0           # Operating frequency [MHz]
PT_W      = 1000.0         # Transmit power [W]
GT        = 1.0            # TX gain (isotropic, linear)
GR        = 1.0            # RX gain (isotropic, linear)

# ── Background grid (IRI + TID density field, coarse) ────────────────────────
# Ray tracer interpolates from this grid.
# Margins on both sides capture ray paths that swing beyond TX/RX.
BG_X_MIN = -50.0     # [km]
BG_X_MAX = 1300.0    # [km]
BG_DX    = 5.0       # [km]
BG_Z_MIN = 60.0      # [km]  IRI lower limit
BG_Z_MAX = 600.0     # [km]
BG_DZ    = 2.0       # [km]

BG_X = np.arange(BG_X_MIN, BG_X_MAX + BG_DX, BG_DX)
BG_Z = np.arange(BG_Z_MIN, BG_Z_MAX + BG_DZ, BG_DZ)

# ── IRI background ────────────────────────────────────────────────────────────
IRI_DT  = datetime(2020, 6, 1, 12, 0)  # Local noon, June, moderate solar activity
IRI_LAT = 32.5    # Mid-path latitude  [°N]  (TX@30°N, path ~1169 km northward → mid ~32.5°N)
IRI_LON = 120.0   # Mid-path longitude [°E]

# ── TID parameters  (Hooke 1968 / Koval 2018 MSTID) ──────────────────────────
TID = {
    'enable':         False,
    'lambda_h_km':    300.0,              # Horizontal wavelength [km]
    'T_s':            2400.0,             # Period [s] = 40 min (typical MSTID)
    'amplitude':      0.10,               # Peak δNe/Ne₀  (0–1)
    'u_para_ms':      50.0,               # Neutral wind ∥ B  [m/s]
    'I_dip_deg':      50.0,               # Geomagnetic dip angle [deg]  (mid-latitude)
    'H_km':           60.0,               # Chapman scale height [km]
    't_s':            0.0,                # Snapshot time [s]
    'omega_b_rad_s':  2.0 * np.pi / 1200.0,  # Brunt-Väisälä freq [rad/s] (20-min period)
}

# ── Es layer parameters  (Hao et al. 2017) ───────────────────────────────────
ES = {
    'enable':     False,
    'foEs_MHz':   5.0,    # Es plasma frequency  [MHz]
    'h_Es_km':    110.0,  # Es center height     [km]
    'delta_h_m':  115.0,  # Half-thickness       [m]   (Hao 2017 best fit)
    'n_exp':      5,      # Density profile exponent
    'L1_m':       300.0,  # Horizontal irregularity scale (∥ drift) [m]
    'L2_m':       300.0,  # Horizontal irregularity scale (⊥ drift) [m]
    'L3_m':       30.0,   # Vertical irregularity scale   [m]
    'delta_N_N':  0.3,    # Relative density fluctuation ΔN/N
    'fr':         0.25,   # Reflection threshold  foEs/f
    'fs':         0.10,   # Scattering threshold  foEs/f
}

# ── Plasma bubble parameters (analytic Gaussian depletion) ───────────────────
BUBBLE = {
    'enable':    False,
    'delta_max': 0.6,    # Maximum density depletion fraction  (0–1)
    'x0_km':     600.0,  # Bubble centre x  [km]  (~mid-path)
    'z0_km':     350.0,  # Bubble centre height [km]
    'Lx_km':     100.0,  # Horizontal half-width  [km]
    'Lz_km':     150.0,  # Vertical   half-extent [km]
}

# ── Ray tracer  (Part 2) ──────────────────────────────────────────────────────
RT = {
    'ds_km':       0.5,   # RK4 nominal step      [km]
    'ds_min_km':   0.1,   # Adaptive step minimum [km]
    'ds_max_km':   2.0,   # Adaptive step maximum [km]
    'beta_min':    5.0,   # Fan: minimum elevation [deg]
    'beta_max':    85.0,  # Fan: maximum elevation [deg]
    'n_fan':       33,    # Fan: number of rays
    'z_stop_km':   550.0, # Abort ray if z exceeds this [km]
    'max_steps':   4000,  # Maximum RK4 steps per ray
}

# ── Point-to-point solver  (Part 3, Nosikov 2020) ────────────────────────────
P2P = {
    'n_init':      18,    # Number of initial-guess paths (elevation scan)
    'n_ctrl':      30,    # Interior control points per path (endpoints fixed)
    'alpha_km':    0.5,   # Gradient-descent step size  [km/iter]
    'k_spring':    0.1,   # Smoothness spring coefficient
    'max_iter':    500,   # Maximum iterations
    'tol_km':      0.05,  # Convergence criterion: ||dS/dr||_perp  [km]
    'clust_h_km':  10.0,  # Deduplication height tolerance   [km]
    'clust_tau_ms':0.05,  # Deduplication delay  tolerance   [ms]
    'n_workers':   1,     # Parallel workers: 0=auto-detect, 1=sequential
                          # Set >1 only in scripts guarded by if __name__=='__main__'
}

# ── PE / SSF  (Part 5, Carrano 2020) ─────────────────────────────────────────
# λ @ 10 MHz = 30 m   →  Nyquist: dz ≤ λ/2 = 15 m
#                         recommended:   dz = λ/4 = 7.5 m
_LAM_M = C_MS / (FREQ_MHZ * 1e6)   # free-space wavelength [m]

PE = {
    'dx_km':          0.5,          # Propagation step       [km]
    'dz_m':           _LAM_M / 4,  # Vertical sampling      [m]  = 7.5 m @10 MHz
    'n_pml':          60,           # PML layer thickness    [grid points]
    'sigma_pml':      0.4,          # PML max attenuation coefficient
    'w0_km':          20.0,         # RT->PE Gaussian beam waist [km]
    'z_pe_min_km':    200.0,        # PE domain lower bound  [km]
    'z_pe_max_km':    500.0,        # PE domain upper bound  [km]
    'store_history':  False,        # Store all x-slices (debug; uses ~GB RAM)
    'min_power_frac': 0.01,         # Mode extraction: min power / peak power
    'earth_flat':     True,         # Apply n_eff = n*(1+z/R_E) Earth-flattening correction
}

# ── Mode classification  (Part 6) ────────────────────────────────────────────
MODE = {
    'h_Es_km':  140.0,  # h_r < 140 km  -> Es mode
    'h_E_km':   200.0,  # 140 <= h_r < 200 km -> E-layer mode
    # h_r >= 200 km -> F-layer mode (further split by delay into 1F-low / 1F-high / 2F)
}

# ── OTH Radar parameters (Phase 1) ───────────────────────────────────────────
RADAR = {
    'mode':             'monostatic',  # 'monostatic' | 'bistatic'
    'target_range_km':  1169.0,        # TX -> target one-way distance [km]
    'sigma_rcs_m2':     5.0,           # Target RCS [m^2] (aircraft resonance region @10 MHz)
    'two_way':          True,          # Output two-way delay tau_2way_ms
}

# TX -> target bearing (arbitrary 0 deg for testing; update for real scenario)
LINK_BEARING_DEG = 0.0   # [deg], clockwise from North

# ── Geomagnetic parameters (Phase 2) ─────────────────────────────────────────
# Computed via ppigrf (IGRF-14) at IRI mid-path point (32.5N, 120E, 300 km, 2020-01-01).
# Re-run utils.get_geomag(IRI_LAT, IRI_LON, dt=IRI_DT) to update for a new path.
GEOMAG = {
    'fH_MHz':    1.197,   # Gyrofrequency [MHz]  (B=42766 nT, ppigrf IGRF-14)
    'dip_deg':   48.7,    # Magnetic dip angle [deg]  (inclination, positive down)
    'decl_deg':  -5.5,    # Magnetic declination [deg] (negative = westward)
    'enable_OX': False,   # Master switch: True -> O/X mode splitting via AH equation
}
