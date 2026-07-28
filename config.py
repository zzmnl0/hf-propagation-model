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

# ── Formal PHaRLAP-compatible spherical link contract ────────────────────────
EARTH_MODEL = {
    'name': 'sphere',
    'radius_km': RE_KM,
    'coordinate_system': 'geocentric_spherical',
    'height_reference': 'above_spherical_earth',
}
TX_LAT    = 30.0           # Geocentric latitude [deg]
TX_LON    = 120.0          # East longitude [deg]
TX_HEIGHT_KM = 0.0
RX_RANGE  = 1169.0         # Derived great-circle distance [km]
LINK_BEARING_DEG = 0.0     # Clockwise from geocentric North [deg]


def _config_spherical_destination(lat_deg, lon_deg, bearing_deg, distance_km):
    lat, lon, bearing = np.radians([lat_deg, lon_deg, bearing_deg])
    angle = float(distance_km) / RE_KM
    lat2 = np.arcsin(np.sin(lat) * np.cos(angle)
                     + np.cos(lat) * np.sin(angle) * np.cos(bearing))
    lon2 = lon + np.arctan2(
        np.sin(bearing) * np.sin(angle) * np.cos(lat),
        np.cos(angle) - np.sin(lat) * np.sin(lat2))
    return float(np.degrees(lat2)), float((np.degrees(lon2) + 180.0) % 360.0 - 180.0)


RX_LAT, RX_LON = _config_spherical_destination(
    TX_LAT, TX_LON, LINK_BEARING_DEG, RX_RANGE)
RX_HEIGHT_KM = 0.0
TX_GEOCENTRIC = (TX_LAT, TX_LON, TX_HEIGHT_KM)
RX_GEOCENTRIC = (RX_LAT, RX_LON, RX_HEIGHT_KM)

# Archived local-link coordinates are retained only for legacy/ reproduction.
# Formal production entry points must not consume these aliases.
TX_POS = (0.0, 0.0)
RX_POS = (RX_RANGE, 0.0)
FREQ_MHZ  = 10.0           # Operating frequency [MHz]
PT_W      = 1000.0         # Transmit power [W]
GT        = 1.0            # TX gain (isotropic, linear)
GR        = 1.0            # RX gain (isotropic, linear)

# ── Background grid (IRI + TID density field, coarse) ────────────────────────
# Ray tracer interpolates from this grid.
# Margins on both sides capture ray paths that swing beyond TX/RX.
BG_X_MIN = -50.0     # [km]
BG_X_MAX = 2000.0    # [km]
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
    'amplitude':      0.10,               # Peak dNe/Ne0  (0-1)
    'u_para_ms':      50.0,               # Neutral wind || B  [m/s]
    'I_dip_deg':      50.0,               # Geomagnetic dip angle [deg]  (mid-latitude)
    'H_km':           60.0,               # Chapman scale height [km]
    't_s':            0.0,                # Snapshot time [s]
    'omega_b_rad_s':  2.0 * np.pi / 1200.0,  # Brunt-Vaisala freq [rad/s] (20-min period)
    # Phase 5: multi-direction superposition (n_components=1 -> backward compat)
    'n_components':     1,                # Number of TID wave components
    'az_deg_list':      [0.0],            # Propagation azimuth per component [deg from North]
    'amplitude_list':   [0.10],           # Amplitude per component (0-1)
    'period_s_list':    [2400.0],         # Period per component [s]
    'lambda_h_km_list': [300.0],          # Horizontal wavelength per component [km]
    'link_bearing_deg': 0.0,              # Link bearing used to project wave vectors onto 2D plane
}

# ── Multi-direction TID for 3-D scenarios  (Huang & Kelley 1996; Ding et al. 2011) ─
# Superposition of 3 MSTID components from different propagation azimuths.
# Dominant SW component (az=0 relative to link) + two oblique components (+/-40 deg).
# IonosphereModel projects these onto x-z at y=0 for the background Ne snapshot;
# RefractiveIndex3D evaluates the full 3-D phase at every query point.
TID_3D = {
    'enable':           False,
    'n_components':     3,
    'az_deg_list':      [0.0,    40.0,   -40.0],   # [deg] relative to link bearing
    'amplitude_list':   [0.10,   0.05,    0.05],   # dNe/Ne per component
    'period_s_list':    [2400.0, 2700.0, 2100.0],  # [s]   (typical MSTID range)
    'lambda_h_km_list': [300.0,  350.0,  250.0],   # [km]
    'z0_km':            280.0,    # F-region perturbation envelope centre [km]
    'sigma_z_km':       80.0,     # F-region vertical envelope sigma [km]
    'lambda_z_km':      180.0,    # vertical wavelength [km]
    # Shared physical parameters (same as TID)
    'u_para_ms':        50.0,
    'I_dip_deg':        50.0,
    'H_km':             60.0,
    't_s':              0.0,
    'omega_b_rad_s':    2.0 * np.pi / 1200.0,
    'link_bearing_deg': 0.0,
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
    'sigma0_ground_dB': -20.0,         # Ground normalized RCS [dBsm] (sea: -20, land: -25)
}

# ── Flux-tube ray tracer (Phase 4 upgrade, Coleman 1997/1998) ─────────────────
TUBE_TRACER = {
    'delta_beta_deg' : 0.5,    # Fan ray elevation spacing [deg]
    'n_tube_rays'    : 80,     # Max rays in fan (capped at this count)
    'x_tgt_tol_km'  : 80.0,   # Landing-point search window half-width [km]
    'T_pulse_ms'     : 0.5,   # Radar pulse width [ms] (pulse correction ref)
    'newton_tol_km'  : 1.0,   # Newton refinement convergence threshold [km]
    'newton_max_iter': 5,      # Newton max iterations (2-3 usually sufficient)
    'L_cross_km'     : 100.0, # Cross-range beam footprint [km] (2D->3D area)
}

# ── Phase 6 upgrades: 3-D ray tracing + scatter theory ───────────────────────

# SCATTER: Born / Rytov / MPS spectral parameters
# Cs_rel  : relative rms fluctuation amplitude delta_Ne/Ne0 (dimensionless)
#            Same convention as SPREAD_F.Cs (used by spread_f_model.py).
# p       : 3-D von Karman spectral index (typ. 2.5-4.0; F-layer ~3.0)
# L_outer_km : outer irregularity scale (correlation length) [km]
# L_inner_m  : inner scale cutoff [m] (default 0.1 m ~ thermal ion gyro-radius)
# DeltaL_km  : default scatter layer thickness [km] (used when not known from ray)
# z_eff_km   : effective observer distance for Rytov Fresnel filter [km]
#              (set to F2 peak height as default; ray tracer updates per-ray)
SCATTER = {
    'Cs_rel':      0.01,    # 1% rms relative fluctuation (moderate F-layer)
    'p':           3.0,     # spectral index (Rino 1979 F-layer value)
    'L_outer_km':  10.0,    # outer scale [km]
    'L_inner_m':   0.1,     # inner scale [m]
    'DeltaL_km':   100.0,   # default layer thickness [km]
    'z_eff_km':    300.0,   # default observer distance for Fresnel filter [km]
}

# RT_3D: parameters for 3-D ray tracer (extends RT for 3-D tracing)
RT_3D = {
    'ds_km':       0.5,    # RK4 step [km]
    'ds_min_km':   0.1,
    'ds_max_km':   2.0,
    'beta_min':    5.0,    # elevation fan min [deg]
    'beta_max':    85.0,
    'n_fan':       33,
    'z_stop_km':   550.0,
    'max_path_km': 8000.0,
    'max_steps':   4000,
}

# P2P_3D: parameters for 3-D P2P Newton shooter
P2P_3D = {
    'seed_source': 'internal_scan',  # internal_scan | fan_warm_start | explicit_list
    'n_init':      18,     # full-domain coarse landing-map samples
    'tol_km':      5.0,    # Newton convergence tolerance [km] (looser than 2-D)
    'd_beta':      0.05,   # Jacobian step: elevation [deg]
    'd_az':        0.10,   # Jacobian step: azimuth   [deg]
'jacobian':    'central_reference_multiscale',
'jacobian_scales': (0.5, 1.0, 2.0),
'jacobian_min_scale': 0.125,     # refine only when 0.5/1/2 has no plateau
'jacobian_scale_rtol': 1.0e-2,  # adjacent central differences must form a plateau
'jacobian_condition_max': 1.0e10,
'jacobian_rtol': 1.0e-3,       # experimental variational/reference diagnostic only
    'max_iter':    20,     # trust-region residual evaluations per candidate
    'candidate_nearest': 3,
    'clust_h_km':  10.0,   # deduplication height tolerance
    'clust_tau_ms': 0.05,  # deduplication delay tolerance
    'clust_launch_deg': 0.02,
    'clust_trajectory_rms_km': 1.0,
    'launch_low_max_deg': 25.0,   # cause_launch low/medium threshold
    'launch_high_min_deg': 55.0,  # cause_launch medium/high threshold
    'seed_betas_deg': [],         # used only when seed_source='explicit_list'
    # A3Rnew-P2 completeness audit.  The normal scenario drivers keep the
    # established internal scan; frozen validation runs opt in explicitly.
    'audit_mode': 'internal_scan',  # internal_scan | adaptive_2d
    'audit_beta_bounds_deg': (2.0, 88.0),
    'audit_az_bounds_deg': (-180.0, 180.0),
    'audit_beta_step_deg': 5.0,
    'audit_az_step_deg': 15.0,
    'audit_max_refinement': 3,
    'audit_target_near_km': 150.0,
    'audit_fold_condition': 1.0e4,
    'audit_checkpoint_dir': 'output2/a3rnew_p2',
    'audit_resume': True,
    'audit_progress': True,
}

# ── Phase 7 upgrades: 3-D flux-tube, field-aligned bubble, BG Es, phase screen
# ─────────────────────────────────────────────────────────────────────────────

# B4: parameterized field-aligned F-region depletion tube.
# The smooth bottomside gate prevents an F-region depletion from extending
# into the E/D region. This is not a dynamical plasma-bubble plume model.
# sigma_perp_km : half-width perpendicular to B  (~50 km, cross-field)
# sigma_para_km : half-length parallel to B      (~550 km, field-aligned)
BUBBLE_3D = {
    'enable':        False,
    'delta_max':     0.6,        # maximum depletion fraction  (0-1)
    'x0_km':         600.0,      # bubble centre along link    [km]
    'y0_km':         0.0,        # bubble centre cross-path    [km]
    'z0_km':         350.0,      # bubble centre height        [km]
    'sigma_perp_km': 50.0,       # perpendicular-to-B half-width [km]
    'sigma_para_km': 550.0,      # parallel-to-B half-length    [km]
    'bottomside_height_km': 200.0,    # F-region bottomside gate midpoint [km]
    'bottomside_transition_km': 20.0, # smooth gate transition scale [km]
}

# A3: 3-D flux-tube tracer (Coleman 1997/1998 extended to 3-D)
# delta_az_deg : azimuth perturbation step used to compute F_focus_az [deg]
TUBE_3D = {
    'delta_beta_deg': 0.5,       # elevation fan spacing [deg]
    'delta_az_deg':   0.5,       # azimuth perturbation for F_focus_az [deg]
    'n_tube_rays':    80,        # max rays in elevation fan
    'x_tgt_tol_km':  80.0,      # search window half-width around x_tgt [km]
    'T_pulse_ms':    0.5,        # radar pulse width [ms]
    'L_cross_km':    100.0,      # cross-range beam footprint [km]
}

# C4: Booker-Gordon anisotropic Es scatter parameters (Booker 1956)
# L_perp_m : irregularity scale perpendicular to B  (horizontal drift scale)
# L_para_m : irregularity scale along B             (field-aligned elongation)
# p_es     : power-law spectral index for Es (shallower than F-layer p~3)
ES_BG = {
    'L_perp_m':  300.0,          # transverse irregularity scale [m]
    'L_para_m':  10000.0,        # field-aligned scale [m]
    'p_es':      1.5,            # Es spectral index
}

# C3/C5: Rino phase screen (C3) and Knepp MPS (C5)
# N_screen : number of phase screens for MPS (C5 only)
# n_pts    : transverse grid points per dimension (square grid)
# dx_m     : transverse grid spacing [m]
# seed     : random seed for phase screen (None = random)
PHASE_SCREEN = {
    'N_screen': 10,              # MPS screen count
    'n_pts':    512,             # transverse grid size
    'dx_m':     200.0,           # transverse spacing [m]
    'seed':     None,            # random seed
}

# ── Phase 5 upgrades: 3D ionospheric background ──────────────────────────────

# B1: IRI lateral sampling along the link great-circle path
IRI_LATERAL = {
    'enable':      False,   # False -> broadcast single profile (backward compat)
    'spacing_km':  50.0,    # Sample interval along great-circle [km]
}

# B2: position-dependent IGRF via GeomagnField2D
GEOMAG_3D = {
    'enable': False,   # False -> use fixed GEOMAG dict (backward compat)
}

# ── D-layer absorption (Phase 3) ─────────────────────────────────────────────
ABSORPTION = {
    'enable':   False,
    'A0':       500.0,   # [dB*MHz^2]  Pederick & Cervera 2014
    'chi_deg':  60.0,    # solar zenith angle [deg]  (60 = daytime mid-latitude)
}

# ── update2.1 feature gates (default off for backward compatibility) ──────────
UPDATE2 = {
    'AH3D_ENABLE': False,
    'AH3D_SOLVER': 'ecef_dop853',
    'AH3D_INTEGRATION_PROFILE': 'research',
    'AH3D_INTEGRATION': {
        'smoke': {'rtol': 3.0e-5, 'atol_r_km': 1.0e-3, 'atol_k': 1.0e-7,
                  'ds_min_km': 0.05, 'ds_max_km': 5.0},
        'research': {'rtol': 3.0e-7, 'atol_r_km': 1.0e-5, 'atol_k': 1.0e-9,
                     'ds_min_km': 0.01, 'ds_max_km': 2.0},
        'reference': {'rtol': 3.0e-9, 'atol_r_km': 1.0e-7, 'atol_k': 1.0e-11,
                      'ds_min_km': 0.002, 'ds_max_km': 0.5},
    },
    'OX_MODE': 'ah3d',
    'RAY_FAN_INCLUDE_X': False,
    'EXTERNAL_NE_PROVIDER': {
        'enable': False,
        'source': 'internal_parameterized',
        'path': None,
        'coordinate_system': 'geocentric_spherical',
        'axis_order': 'latitude,longitude,height',
        'units': 'm^-3',
        'earth_radius_km': RE_KM,
    },
    'PE_COUPLING': {
        'enable': False,
        'apply_to_power': False,
        'mode': 'local_2p5d',
        'q_boundary': 1.0,
    },
    'SCATTER_PD': {
        'apply_to_pd': False,
        'kernel': 'diagnostic',
    },
    'ABSORPTION_COMPLEX': {
        'enable': False,
        'mode': 'legacy',
        'collision_model': 'd_region_exponential',
        'nu_ref_s': 5.0e6,
        'h_ref_km': 80.0,
        'scale_height_km': 8.0,
        'h_min_km': 50.0,
        'h_max_km': 140.0,
    },
    'RADAR_CLUTTER': {
        'enable': False,
        'apply_to_power': False,
        'surface_type': 'sea',
        'sigma0_dB': None,
        'radial_velocity_ms': 0.0,
        'grazing_angle_deg': None,
        'clutter_class': 'surface_bragg',
    },
}

# ── Spread-F parameters (Phase 3, Rino 1979) ─────────────────────────────────
SPREAD_F = {
    'enable':       False,
    'Cs':           1e-3,   # phase spectral strength
    'p':            3.0,    # power-law spectral index (typical 2.5-4.0)
    'h_screen_km':  300.0,  # phase screen height [km]
    'L0_km':        50.0,   # outer scale [km]
    'l0_m':         100.0,  # inner scale [m]
    'delta_Ne_rel_rms': 0.05,  # relative Ne RMS; distinct from spectral strength Cs
    'sigma_z_km':   40.0,   # F-region irregularity envelope sigma [km]
    'anisotropy_ratio': 20.0,
    'n_modes':      96,
    'seed':         0,
}

# ── Geomagnetic parameters (Phase 2) ─────────────────────────────────────────
# Computed via ppigrf (IGRF-14) at IRI mid-path point (32.5N, 120E, 300 km, 2020-01-01).
# Re-run utils.get_geomag(IRI_LAT, IRI_LON, dt=IRI_DT) to update for a new path.
GEOMAG = {
    'fH_MHz':    1.197,   # Gyrofrequency [MHz]  (B=42766 nT, ppigrf IGRF-14)
    'dip_deg':   48.7,    # Magnetic dip angle [deg]  (inclination, positive down)
    'decl_deg':  -5.5,    # Magnetic declination [deg] (negative = westward)
    'enable_OX': False,   # Master switch: True -> O/X mode splitting via AH equation
}



