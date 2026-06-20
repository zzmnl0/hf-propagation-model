"""
Module C0: Scatter method selector via phase-variance criterion.

For a ray passing through an ionospheric scattering layer, computes the
rms phase variance sigma_phi^2 [rad^2] accumulated through that layer.
This serves as the dimensionless discriminator analogous to Rino's (1979)
strong-scatter parameter U:

  sigma_phi^2 = (r_e * lambda)^2 * delta_Ne_rms^2 * Delta_L * L_outer

where:
  r_e        = 2.818e-15 m   (classical electron radius)
  lambda     = c / f         (free-space wavelength [m])
  delta_Ne_rms = Cs_rel * Ne0  (rms electron density fluctuation [m^-3])
  Delta_L    = scatter layer thickness along ray path [m]
  L_outer    = outer irregularity scale (decorrelation length) [m]

Physical derivation:
  Phase accumulated through path element dz: dphi = r_e * lambda * delta_Ne * dz
  For incoherent addition over Delta_L with correlation length L_outer:
  <phi^2> = r_e^2 * lambda^2 * delta_Ne_rms^2 * Delta_L * L_outer [rad^2]
  Ref: Tatarski (1961) "Wave Propagation in a Turbulent Medium",
       Yeh & Liu (1982) Proc. IEEE 70(4), 324-360.

Selection thresholds (Rino 1979 / Yeh & Liu 1982):
  sigma_phi^2 < 0.1      -> 'Born'   (weak single scatter)
  0.1 <= sigma_phi^2 < 10 -> 'Rytov' (moderate, forward-scatter dominated)
  sigma_phi^2 >= 10       -> 'MPS'   (strong scatter, Knepp 1983)

Typical values at 10 MHz:
  F-layer weak turb  (Cs=0.01, Ne0=5e10, DL=100km, L0=10km): sigma_phi^2 ~ 1.8  -> Rytov
  Bubble strong      (Cs=0.30, Ne0=5e10, DL=100km, L0=50km): sigma_phi^2 ~ 1e4  -> MPS
  Background IRI     (Cs=0.001, Ne0=1e11, DL=50km, L0=1km):  sigma_phi^2 ~ 2e-4 -> Born
"""
import numpy as np

# Physical constants
R_E_M  = 2.818e-15   # classical electron radius [m]
C_MS   = 2.998e8     # speed of light [m/s]


def scatter_phase_variance(Ne0_m3:      float,
                            Cs_rel:      float,
                            DeltaL_km:   float,
                            L_outer_km:  float,
                            freq_MHz:    float) -> float:
    """
    Compute phase variance sigma_phi^2 [rad^2] for a scattering layer.

    Parameters
    ----------
    Ne0_m3     : background electron density at scatter layer [m^-3]
    Cs_rel     : relative rms fluctuation amplitude delta_Ne/Ne0 (0-1)
    DeltaL_km  : layer thickness along ray path [km]
    L_outer_km : outer irregularity scale (correlation length) [km]
    freq_MHz   : radar frequency [MHz]

    Returns
    -------
    sigma_phi_sq : phase variance [rad^2]
    """
    lam_m     = C_MS / (freq_MHz * 1e6)              # free-space wavelength [m]
    delta_Ne  = Cs_rel * max(Ne0_m3, 1.0)            # rms delta_Ne [m^-3]
    DeltaL_m  = DeltaL_km * 1e3                       # layer thickness [m]
    L_m       = L_outer_km * 1e3                      # outer scale [m]
    return (R_E_M * lam_m) ** 2 * delta_Ne**2 * DeltaL_m * L_m


def select_scatter_method(sigma_phi_sq: float) -> str:
    """
    Choose scatter theory based on phase variance [rad^2].

    Returns 'Born', 'Rytov', or 'MPS'.
    """
    if sigma_phi_sq < 0.1:
        return 'Born'
    elif sigma_phi_sq < 10.0:
        return 'Rytov'
    else:
        return 'MPS'


def diagnose_scatter(ray: dict,
                     Ne0_m3:     float,
                     Cs_rel:     float,
                     DeltaL_km:  float,
                     L_outer_km: float,
                     freq_MHz:   float) -> dict:
    """
    Full scatter diagnosis for a ray: compute sigma_phi^2, select method,
    and return diagnostic dict.

    Parameters
    ----------
    ray        : ray dict (from trace_single_ray_3d or trace_single_ray)
    Ne0_m3     : representative background Ne at scatter layer [m^-3]
    Cs_rel     : relative rms fluctuation amplitude
    DeltaL_km  : layer thickness [km]
    L_outer_km : outer scale [km]
    freq_MHz   : frequency [MHz]

    Returns dict with:
      sigma_phi_sq   : [rad^2]
      sigma_phi_rad  : [rad]
      scatter_method : 'Born' | 'Rytov' | 'MPS'
      U_rino_equiv   : sigma_phi_sq (our Rino-U analog)
    """
    sphi2  = scatter_phase_variance(Ne0_m3, Cs_rel, DeltaL_km, L_outer_km, freq_MHz)
    method = select_scatter_method(sphi2)
    return {
        'sigma_phi_sq'   : sphi2,
        'sigma_phi_rad'  : float(np.sqrt(max(sphi2, 0.0))),
        'scatter_method' : method,
        'U_rino_equiv'   : sphi2,
    }
