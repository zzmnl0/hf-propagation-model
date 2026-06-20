"""
Module C1: Born weak-scatter ionospheric backscatter.

Computes the volume backscatter cross-section sigma_v [m^-1] and the
resulting received scatter power for HF radar using the Born (single-scatter)
approximation (Booker & Gordon 1950, Galushko et al. 2013).

Valid when sigma_phi^2 < 0.1 rad^2 (weak scatter regime).

Born backscatter cross-section per unit volume:
  sigma_v(Delta_k) = r_e^2 * (2*pi)^3 * Phi_Ne(Delta_k)   [m^-1]

where:
  r_e = 2.818e-15 m   (classical electron radius)
  Delta_k = 2*k0      (Bragg condition for backscatter, |Delta_k| = 2*k0 = 4*pi*f/c)
  Phi_Ne(kappa) = C_s_abs * (kappa^2 + kappa_0^2)^(-(p+3)/2) * exp(-kappa^2/kappa_m^2)
                             [von Karman 3-D spectrum, [m^-3]]

Dimensional analysis:
  [sigma_v] = [r_e^2] * [Phi_Ne] = m^2 * m^-3 * (2pi)^3  --> need m^-1
  Factor: (2*pi)^3 * r_e^2 * [Phi_Ne = m^-3] = m^-1 (using FT convention
  where Phi = integral R(r) exp(ik.r) d^3r without 1/(2pi)^3 factor)

Absolute spectral coefficient C_s_abs is derived from the relative
fluctuation amplitude Cs_rel (dimensionless, same as config.SPREAD_F.Cs):
  C_s_abs = Cs_rel^2 * Ne0^2 * p * (2*pi)^2 * kappa_0^p   [m^(-p-6)]

Received scatter power (one-way radar equation, monostatic):
  Pr_scatter = Pt * Gt * Gr * lambda^2 / (4*pi)^3 * sigma_v * V_eff / R^4

where V_eff is the resolution volume and R is the one-way range.

References:
  Booker & Gordon (1950) Proc. IRE 38:401-412
  Galushko et al. (2013) Radio Sci. 48:577-586
  Ishimaru (1978) "Wave Propagation and Scattering in Random Media" Vol. 2
"""
import numpy as np
from scipy.integrate import quad

R_E_M  = 2.818e-15   # classical electron radius [m]
C_MS   = 2.998e8     # speed of light [m/s]
TWO_PI = 2.0 * np.pi


def _c_s_abs(Cs_rel: float,
             Ne0_m3: float,
             p_spec: float,
             L_outer_km: float) -> float:
    """
    Compute absolute spectral coefficient C_s_abs [m^(-p-6)] from
    relative fluctuation amplitude Cs_rel (dimensionless):

      C_s_abs = Cs_rel^2 * Ne0^2 * p * (2*pi)^2 * kappa_0^p

    Derivation: normalise von Karman variance integral
      integral Phi_Ne d^3k / (2*pi)^3 = delta_Ne_rms^2 = (Cs_rel * Ne0)^2
    for kappa >> kappa_0 (outer-scale dominated variance).
    """
    kappa_0 = TWO_PI / (L_outer_km * 1e3)              # [rad/m]
    return (Cs_rel * Ne0_m3)**2 * p_spec * TWO_PI**2 * kappa_0**p_spec


def phi_von_karman(kappa:       float,
                   C_s_abs:     float,
                   p_spec:      float,
                   kappa_0:     float,
                   kappa_m:     float = 0.0) -> float:
    """
    3-D von Karman power spectrum of Ne fluctuations [m^-3].

    Phi_Ne(kappa) = C_s_abs * (kappa^2 + kappa_0^2)^(-(p+3)/2)
                            * exp(-kappa^2 / kappa_m^2)

    kappa_m = 2*pi / L_inner (inner-scale cutoff; 0 = no inner-scale cutoff)
    """
    Phi = C_s_abs * (kappa**2 + kappa_0**2) ** (-(p_spec + 3.0) / 2.0)
    if kappa_m > 0.0:
        Phi *= np.exp(-(kappa / kappa_m)**2)
    return max(Phi, 0.0)


def born_sigma_v(freq_MHz:    float,
                 Cs_rel:      float,
                 p_spec:      float,
                 L_outer_km:  float,
                 Ne0_m3:      float,
                 L_inner_m:   float = 0.1) -> float:
    """
    Born backscatter cross-section per unit volume sigma_v [m^-1].

    For monostatic HF radar backscatter (Bragg condition: Delta_k = 2*k0):
      sigma_v = r_e^2 * (2*pi)^3 * Phi_Ne(2*k0)

    Parameters
    ----------
    freq_MHz   : radar frequency [MHz]
    Cs_rel     : relative rms fluctuation delta_Ne/Ne0 (dimensionless)
    p_spec     : 3-D spectral index p (von Karman exponent; typ. 2.5-4.0)
    L_outer_km : outer irregularity scale [km]
    Ne0_m3     : local background electron density [m^-3]
    L_inner_m  : inner-scale cutoff [m] (0.1 m default = no effective cutoff)

    Returns
    -------
    sigma_v : [m^-1]
    """
    k0      = TWO_PI * freq_MHz * 1e6 / C_MS           # free-space wavenumber [rad/m]
    kappa_B = 2.0 * k0                                  # Bragg wavenumber [rad/m]
    kappa_0 = TWO_PI / (L_outer_km * 1e3)              # outer-scale cutoff [rad/m]
    kappa_m = TWO_PI / max(L_inner_m, 1e-6)            # inner-scale cutoff [rad/m]
    C_abs   = _c_s_abs(Cs_rel, Ne0_m3, p_spec, L_outer_km)
    Phi     = phi_von_karman(kappa_B, C_abs, p_spec, kappa_0, kappa_m)
    return R_E_M**2 * TWO_PI**3 * Phi


def born_scatter(ray:         dict,
                 Pt_W:        float,
                 Gt:          float,
                 Gr:          float,
                 freq_MHz:    float,
                 Cs_rel:      float,
                 p_spec:      float,
                 L_outer_km:  float,
                 Ne0_m3:      float,
                 beam_width_km: float = 50.0) -> dict:
    """
    Compute Born scatter output fields for one ray.

    Scatter power integrated along the ray path:
      Pr_scatter = Pt*Gt*Gr*lambda^2 / (4*pi)^3
                   * sigma_v * (beam_width_km*1e3)^2 * ds * N_steps / R^4

    where R = group_path_km (one-way), ds = step size [km], and
    beam_width_km^2 * ds is the resolution volume element.

    New output fields:
      Pr_scatter_W   : Born scatter received power [W]
      sigma_v_m1     : volume backscatter coefficient [m^-1]
      tau_scatter_ms : group delay at mid-scatter (same as ray tau_ms for backscatter)
      az_scatter_deg : estimated azimuth spread [deg] from Bragg condition geometry
      scatter_method : 'Born'
    """
    lam_m   = C_MS / (freq_MHz * 1e6)
    R_m     = ray['group_path_km'] * 1e3
    sv      = born_sigma_v(freq_MHz, Cs_rel, p_spec, L_outer_km, Ne0_m3)

    # Resolution volume: cross-section = beam_width^2, depth = range resolution
    # Use full-path depth as upper bound (conservative)
    V_eff_m3 = (beam_width_km * 1e3)**2 * (R_m / 3.0)   # ~1/3 path in ionosphere

    Pr_W  = (Pt_W * Gt * Gr * lam_m**2 / (4.0 * np.pi)**3
             * sv * V_eff_m3 / max(R_m**4, 1.0))

    # Azimuth spread from Bragg geometry: delta_az ~ lambda / (2 * L_outer)
    az_spread = float(np.degrees(lam_m / (2.0 * L_outer_km * 1e3)))

    return {
        'Pr_scatter_W'   : float(Pr_W),
        'Pr_scatter_dBW' : float(10.0 * np.log10(max(Pr_W, 1e-30))),
        'sigma_v_m1'     : float(sv),
        'tau_scatter_ms' : float(ray['tau_ms']),
        'az_scatter_deg' : float(az_spread),
        'scatter_method' : 'Born',
    }
