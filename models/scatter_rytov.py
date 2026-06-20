"""
Module C2: Rytov scintillation theory for ionospheric HF scatter.

Computes scintillation index S4, phase jitter sigma_phi, and time-delay
jitter sigma_tau using the Rytov approximation (weak-to-moderate scatter).

Valid range: S4 < 0.6 (i.e., sigma_phi^2 < ~10 rad^2).
For S4 > 0.6 use the Rino (1979 Part 2) modified Rytov or MPS (Phase 7).

Intensity scintillation index S4:
  S4^2 = sigma_I^2 / <I>^2 (normalised intensity variance)

Rytov formula (Yeh & Liu 1982, Eq. 10.42; Rino 1979 Part 1):
  S4^2 = 2 * k0^4 * r_e^2 * Delta_L * integral_0^inf
           Phi_Ne(k_perp, kz=0) * sin^2(k_perp^2 * z_eff / (2*k0)) * k_perp dk_perp

where:
  k0      = 2*pi*f/c   (free-space wavenumber [rad/m])
  r_e     = classical electron radius [m]
  Delta_L = scattering layer thickness (path length through layer) [m]
  z_eff   = propagation distance from scatter layer to observer [m]
            (= h_scatter / cos(zenith_angle) for oblique paths)
  Phi_Ne(k_perp, kz=0) = C_s_abs * (k_perp^2 + kappa_0^2)^(-(p+3)/2)
                          (von Karman 3-D spectrum at kz=0, kappa=k_perp)
  Fresnel filter: sin^2(k_perp^2 * z_eff / (2*k0))

Phase jitter:
  sigma_phi^2 = (r_e * lambda)^2 * delta_Ne_rms^2 * Delta_L * L_outer  [rad^2]
  (identical to scatter_selector; included here for completeness)

Time-delay jitter:
  sigma_tau = sigma_phi / (2*pi*f)   [s]  (dispersion-free approximation)

Fresnel radius:
  r_F = sqrt(lambda * z_eff / (2*pi))   [m]
  Dominant scatter from scales near r_F.

References:
  Rino (1979) Radio Sci. 14(6), 1135-1145 (Part 1: weak scatter).
  Yeh & Liu (1982) Proc. IEEE 70(4), 324-360.
  Tatarski (1961) "Wave Propagation in a Turbulent Medium".
"""
import numpy as np
from scipy.integrate import quad

R_E_M  = 2.818e-15   # classical electron radius [m]
C_MS   = 2.998e8     # speed of light [m/s]
TWO_PI = 2.0 * np.pi


def rytov_s4(freq_MHz:    float,
             Cs_rel:      float,
             p_spec:      float,
             L_outer_km:  float,
             Ne0_m3:      float,
             z_eff_km:    float,
             DeltaL_km:   float) -> float:
    """
    Rytov scintillation index S4 (intensity).

    S4^2 = 2*k0^4*r_e^2 * C_s_abs * Delta_L
           * integral_kappa0^k_max
             (k_perp^2+kappa_0^2)^(-(p+3)/2) * sin^2(k_perp^2*z_eff/2k0)
             * k_perp dk_perp

    Integration from outer-scale wavenumber kappa_0 to k0 (Nyquist of HF).

    Parameters
    ----------
    freq_MHz   : radar frequency [MHz]
    Cs_rel     : relative rms fluctuation delta_Ne/Ne0
    p_spec     : 3-D spectral index
    L_outer_km : outer irregularity scale [km]
    Ne0_m3     : background Ne at scatter layer [m^-3]
    z_eff_km   : effective propagation distance from layer to observer [km]
    DeltaL_km  : scatter layer thickness [km]

    Returns
    -------
    S4 : scintillation index (0 to ~1; capped at 1 for strong scatter)
    """
    from .scatter_born import _c_s_abs

    k0       = TWO_PI * freq_MHz * 1e6 / C_MS       # [rad/m]
    kappa_0  = TWO_PI / (L_outer_km * 1e3)          # outer scale [rad/m]
    z_eff_m  = z_eff_km * 1e3                        # [m]
    DeltaL_m = DeltaL_km * 1e3                       # [m]
    C_abs    = _c_s_abs(Cs_rel, Ne0_m3, p_spec, L_outer_km)

    def integrand(k_perp):
        Phi    = C_abs * (k_perp**2 + kappa_0**2) ** (-(p_spec + 3.0) / 2.0)
        filt   = np.sin(k_perp**2 * z_eff_m / (2.0 * k0)) ** 2
        return Phi * filt * k_perp

    # Integrate from kappa_0 to k0 (HF can only sense scales down to ~lambda/2)
    k_max = k0
    try:
        I, _ = quad(integrand, kappa_0, k_max,
                    limit=200, epsabs=1e-40, epsrel=1e-6)
    except Exception:
        I = 0.0

    S4_sq = 2.0 * k0**4 * R_E_M**2 * DeltaL_m * I
    return float(min(np.sqrt(max(S4_sq, 0.0)), 1.0))


def rytov_sigma_phi(freq_MHz:    float,
                    Cs_rel:      float,
                    Ne0_m3:      float,
                    DeltaL_km:   float,
                    L_outer_km:  float) -> float:
    """
    RMS phase jitter sigma_phi [rad] from Rytov approximation.

    sigma_phi^2 = (r_e * lambda)^2 * (Cs_rel*Ne0)^2 * Delta_L * L_outer
    Ref: Tatarski (1961) Ch. 6; consistent with scatter_selector formula.
    """
    lam_m    = C_MS / (freq_MHz * 1e6)
    delta_Ne = Cs_rel * max(Ne0_m3, 1.0)
    DeltaL_m = DeltaL_km * 1e3
    L_m      = L_outer_km * 1e3
    sphi_sq  = (R_E_M * lam_m)**2 * delta_Ne**2 * DeltaL_m * L_m
    return float(np.sqrt(max(sphi_sq, 0.0)))


def rytov_sigma_tau(freq_MHz:    float,
                    Cs_rel:      float,
                    Ne0_m3:      float,
                    DeltaL_km:   float,
                    L_outer_km:  float) -> float:
    """
    RMS time-delay jitter sigma_tau [ms] from phase jitter.

    sigma_tau = sigma_phi / omega_0   [s]
    omega_0 = 2*pi*f
    Dispersion-free approximation (valid when pulse bandwidth << f).
    """
    omega0   = TWO_PI * freq_MHz * 1e6
    sig_phi  = rytov_sigma_phi(freq_MHz, Cs_rel, Ne0_m3, DeltaL_km, L_outer_km)
    return float(sig_phi / omega0 * 1e3)   # [ms]


def fresnel_radius_km(freq_MHz: float, z_eff_km: float) -> float:
    """
    Fresnel scale r_F = sqrt(lambda * z_eff / (2*pi)) [km].
    Dominant scatter contribution from irregularities near this scale.
    """
    lam_m  = C_MS / (freq_MHz * 1e6)
    r_F_m  = np.sqrt(lam_m * z_eff_km * 1e3 / TWO_PI)
    return r_F_m / 1e3


def rytov_full(freq_MHz:    float,
               Cs_rel:      float,
               p_spec:      float,
               L_outer_km:  float,
               Ne0_m3:      float,
               z_eff_km:    float,
               DeltaL_km:   float) -> dict:
    """
    Compute all Rytov scintillation outputs.

    Returns dict with:
      S4              : intensity scintillation index (0-1)
      sigma_phi_rad   : rms phase jitter [rad]
      sigma_tau_ms    : rms time-delay jitter [ms]
      fresnel_r_km    : Fresnel radius [km]
      scatter_method  : 'Rytov'
    """
    S4      = rytov_s4(freq_MHz, Cs_rel, p_spec, L_outer_km, Ne0_m3,
                       z_eff_km, DeltaL_km)
    sphi    = rytov_sigma_phi(freq_MHz, Cs_rel, Ne0_m3, DeltaL_km, L_outer_km)
    stau_ms = rytov_sigma_tau(freq_MHz, Cs_rel, Ne0_m3, DeltaL_km, L_outer_km)
    r_F_km  = fresnel_radius_km(freq_MHz, z_eff_km)

    return {
        'S4'             : S4,
        'sigma_phi_rad'  : sphi,
        'sigma_tau_ms'   : stau_ms,
        'fresnel_r_km'   : r_F_km,
        'scatter_method' : 'Rytov',
    }
