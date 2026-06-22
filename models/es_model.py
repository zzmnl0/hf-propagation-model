"""
Module M3 – Es-layer reflection / scattering  (Hao et al. 2017).
Three-segment model: pure reflection | mixed | pure scattering,
determined by the ratio foEs/f.

compute_power() returns:
    {
      'Pr_W'           : total received power [W],
      'Pr_reflect_W'   : reflection contribution [W],
      'Pr_scatter_W'   : scattering contribution [W],
      'mode'           : 'reflect' | 'scatter' | 'mixed',
      'alpha'          : reflection weight (0-1),
      'tau_ms'         : group delay (one-way path D_km) [ms],
      'delta_tau_ms'   : delay spread from scattering [ms],
    }
"""
import numpy as np
from math import factorial
from config import C_MS, C_KMS, ES, GEOMAG


class EsLayerModel:
    """
    Es-layer three-segment propagation model (Hao et al. 2017).

    Parameters (all from config.ES if not overridden):
        foEs_MHz  : Es plasma frequency  [MHz]
        h_Es_km   : Es center height     [km]
        delta_h_m : half-thickness       [m]   (Hao 2017 best-fit: 115 m)
        n_exp     : density exponent     (Hao 2017: 5)
        L1_m, L2_m: horizontal irregularity scales [m]
        L3_m      : vertical  irregularity scale   [m]
        delta_N_N : relative density fluctuation dN/N
        fr        : reflection threshold foEs/f  (0.25)
        fs        : scattering threshold foEs/f  (0.10)
    """

    def __init__(self,
                 foEs_MHz:  float = ES['foEs_MHz'],
                 h_Es_km:   float = ES['h_Es_km'],
                 delta_h_m: float = ES['delta_h_m'],
                 n_exp:     int   = ES['n_exp'],
                 L1_m:      float = ES['L1_m'],
                 L2_m:      float = ES['L2_m'],
                 L3_m:      float = ES['L3_m'],
                 delta_N_N: float = ES['delta_N_N'],
                 fr:        float = ES['fr'],
                 fs:        float = ES['fs']):
        self.foEs    = foEs_MHz
        self.h_Es    = h_Es_km
        self.delta_h = delta_h_m
        self.n_exp   = n_exp
        self.L1      = L1_m
        self.L2      = L2_m
        self.L3      = L3_m
        self.dNN     = delta_N_N
        self.fr      = fr
        self.fs      = fs

    # ── Mode classification ───────────────────────────────────────────────────

    def classify(self, f_MHz: float) -> tuple[str, float]:
        """
        Decide propagation regime from foEs/f.

        Returns
        -------
        mode  : 'reflect' | 'scatter' | 'mixed'
        alpha : reflection weight (1.0 for reflect, 0.0 for scatter)
        """
        ratio = self.foEs / f_MHz
        if ratio > self.fr:
            return 'reflect', 1.0
        elif ratio < self.fs:
            return 'scatter', 0.0
        else:
            alpha = (ratio - self.fs) / (self.fr - self.fs)
            return 'mixed', alpha

    # ── Reflection coefficient ────────────────────────────────────────────────

    def reflection_coeff_sq(self, theta_rad: float, f_MHz: float) -> float:
        """
        |rho|^2 – squared reflection coefficient from the Es thin layer.
        Hao (2017) Eq. 3, n_exp = 5.

        theta_rad : grazing angle (complement of incidence angle) [rad]
        f_MHz     : operating frequency [MHz]
        """
        lam = C_MS / (f_MHz * 1e6)          # free-space wavelength [m]
        n   = self.n_exp                     # = 5
        fN  = self.foEs                      # [MHz], same units as f_MHz
        dh  = self.delta_h                   # [m]

        L = 4.0 * np.pi * np.sin(theta_rad) * dh / lam

        # Σ_{k=0}^{n-1} L^(2k)/(2k)! * (-1)^k * (cos L + L*sin L/(2k+1))
        sumval = 0.0
        cosL   = np.cos(L)
        sinL   = np.sin(L)
        for k in range(n):
            coeff   = (L ** (2 * k)) / factorial(2 * k) * ((-1) ** k)
            sumval += coeff * (cosL + L * sinL / (2 * k + 1))

        bracket = ((-1) ** n) * sumval + 1.0
        prefac  = n * factorial(2 * n - 2) * (fN / f_MHz) ** 2
        denom   = (theta_rad ** 2) * (L ** (2 * n)) + 1e-60
        rho_sq  = (prefac / denom * bracket) ** 2
        return rho_sq

    # ── Scatter cross-section ─────────────────────────────────────────────────

    def scatter_cross_section(self, theta_rad: float, f_MHz: float) -> float:
        """
        sigma_5 [m^-1] – anisotropic Es irregularity scatter cross section.
        Hao (2017) Eq. 7, m=5.

        theta_rad : grazing angle [rad]
        f_MHz     : operating frequency [MHz]
        """
        lam = C_MS / (f_MHz * 1e6)          # [m]
        fN  = self.foEs
        L1, L2, L3 = self.L1, self.L2, self.L3
        dNN = self.dNN
        st  = np.sin(theta_rad)

        factor1 = 1.03e2 * dNN ** 2 * (fN / f_MHz) ** 4 * L1 * L2
        bracket = 1.0 + (3.5 * lam / (4.0 * np.pi * L3 * st)) ** 2
        sigma5  = (factor1
                   * bracket ** (-6.5)
                   * (2.0 * np.pi * L3 / lam) ** 9
                   / (L3 ** 3 * st ** 13))
        return sigma5

    # ── Combined three-segment power ──────────────────────────────────────────

    def compute_power(self,
                      Pt_W:      float,
                      Gt:        float,
                      Gr:        float,
                      f_MHz:     float,
                      D_km:      float,
                      theta_rad: float) -> dict:
        """
        Three-segment Es received power (Hao 2017).

        Parameters
        ----------
        Pt_W      : transmit power [W]
        Gt, Gr    : TX / RX antenna gains (linear; isotropic = 1)
        f_MHz     : operating frequency [MHz]
        D_km      : total path length TX->Es->RX [km]
        theta_rad : Es-layer grazing angle [rad]

        Returns a result dict (see module docstring).
        """
        lam  = C_MS / (f_MHz * 1e6)     # [m]
        D_m  = D_km * 1e3               # [m]

        mode, alpha = self.classify(f_MHz)

        # Group delay for the total path
        tau_ms = D_km / C_KMS * 1e3

        # Reflection power
        rho_sq     = self.reflection_coeff_sq(theta_rad, f_MHz)
        Pr_reflect = Pt_W * Gt * Gr * lam ** 2 / ((4.0 * np.pi) ** 2 * D_m ** 2) * rho_sq

        # Scatter power
        sigma5     = self.scatter_cross_section(theta_rad, f_MHz)
        Pr_scatter = (Pt_W * self.delta_h * Gt * Gr * lam ** 2
                      / (2.0 * np.pi ** 2 * D_m ** 2 * np.sin(theta_rad))
                      * sigma5)

        # Weighted combination
        Pr = alpha * Pr_reflect + (1.0 - alpha) * Pr_scatter

        # Delay spread: scatter contributes L3/c_ms RMS delay
        delta_tau_ms = (1.0 - alpha) * (self.L3 / C_MS) * 1e3

        return {
            'Pr_W'         : Pr,
            'Pr_reflect_W' : Pr_reflect,
            'Pr_scatter_W' : Pr_scatter,
            'mode'         : mode,
            'alpha'        : alpha,
            'tau_ms'       : tau_ms,
            'delta_tau_ms' : delta_tau_ms,
        }

    # ── PE interface helper ───────────────────────────────────────────────────

    def transmission_amplitude(self, theta_rad: float, f_MHz: float) -> float:
        """
        Field amplitude transmission coefficient T = sqrt(1 - |rho|^2).
        Clamped to [0, 1] for use as a PE thin-layer boundary condition.
        """
        rho_sq = self.reflection_coeff_sq(theta_rad, f_MHz)
        return float(np.sqrt(max(0.0, min(1.0, 1.0 - rho_sq))))


# ── Booker-Gordon (1956) anisotropic irregularity scatter ────────────────────
#
# Reference: Booker, H.G. (1956). A theory of scattering by nonisotropic
#   irregularities with application to radar reflections from the aurora.
#   J. Atm. Terr. Phys. 8, 204-221.
# Reference: Booker, H.G. & Gordon, W.E. (1950). Proc. IRE, 38, 401-412.
# Reference: Fejer & Kan (1969); Hao (2017) [Hao model kept for comparison].


def booker_gordon_sigma_v(freq_MHz:   float,
                           Cs_rel:     float,
                           Ne_m3:      float,
                           L_perp_m:   float,
                           L_para_m:   float,
                           p_es:       float = 1.5,
                           aspect_deg: float = 0.0) -> float:
    """
    Booker-Gordon (1956) volume scatter cross section [m^-1] for anisotropic
    Es irregularities.

    The anisotropic PSD separates into:
      (a) an isotropic von Karman spectrum in the plane perpendicular to B
      (b) a Gaussian field-aligned aspect filter parallel to B

    Combined formula (Booker 1956; upgrade_plan_3d_scatter Section C4):
      sigma_v = r_e^2 * (2*pi)^3 * C_s_abs
                * (Dk_perp^2 + kappa_0^2)^(-(p+3)/2)
                * exp(-Dk_para^2 * L_para^2 / 4)

    where:
      Dk      = 2*k0            (backscatter Bragg wavenumber [rad/m])
      Dk_perp = Dk * cos(asp)   (component perpendicular to B)
      Dk_para = Dk * sin(asp)   (component along B)
      kappa_0 = 2*pi / L_perp   (outer-scale wavenumber for perpendicular turbulence)
      C_s_abs = Cs_rel^2 * Ne^2 * p * (2*pi)^2 * kappa_0^p  (Born convention)

    The Gaussian aspect filter exp(-Dk_para^2*L_para^2/4) enforces the
    Es "aspect sensitivity": scatter is suppressed unless Dk_para*L_para << 1,
    i.e., unless the Bragg vector is nearly perpendicular to B.

    Parameters
    ----------
    freq_MHz   : operating frequency [MHz]
    Cs_rel     : relative rms density fluctuation delta_Ne / Ne  (0-1)
    Ne_m3      : background electron density [m^-3]
    L_perp_m   : outer irregularity scale perpendicular to B [m]  (kappa_0 = 2pi/L_perp)
    L_para_m   : field-aligned scale for aspect filter [m]  (large -> narrow aspect)
    p_es       : power-law spectral index for Es turbulence (default 1.5)
    aspect_deg : angle between Dk and B-perp plane [deg]
                 (0 = optimal perpendicular scatter; 90 = fully field-aligned)

    Returns
    -------
    sigma_v [m^-1]
    """
    r_e     = 2.818e-15   # classical electron radius [m]
    TWO_PI  = 2.0 * np.pi
    k0      = TWO_PI * freq_MHz * 1e6 / C_MS
    Dk      = 2.0 * k0    # backscatter Bragg wavenumber [rad/m]

    asp_r   = np.radians(float(aspect_deg))
    Dk_perp = Dk * np.cos(asp_r)
    Dk_para = Dk * np.sin(asp_r)

    # Von Karman outer-scale wavenumber and absolute spectral strength
    kappa_0 = TWO_PI / float(L_perp_m)
    C_s_abs = (float(Cs_rel) * float(Ne_m3))**2 * float(p_es) * TWO_PI**2 * kappa_0**float(p_es)

    # Isotropic von Karman spectrum at perpendicular Bragg wavenumber
    phi_iso = C_s_abs * (Dk_perp**2 + kappa_0**2) ** (-0.5 * (float(p_es) + 3.0))

    # Field-aligned (aspect) suppression filter
    gauss   = np.exp(-Dk_para**2 * float(L_para_m)**2 / 4.0)

    return float(r_e**2 * TWO_PI**3 * phi_iso * gauss)


def es_aspect_angle_deg(theta_rad: float,
                         dip_deg:   float = GEOMAG['dip_deg']) -> float:
    """
    Booker-Gordon aspect angle for Es backscatter at grazing angle theta [deg].

    For oblique incidence at grazing angle theta, the Bragg wavevector is
    directed upward at ~90-theta from horizontal.  The aspect angle with
    respect to the horizontal geomagnetic field direction is approximately:
      aspect = |90 - dip - theta_deg|

    For vertical incidence (theta=90 deg) and dip=48.7 deg: aspect = 41.3 deg.
    The exponential suppression exp(-Dk_para^2*L_para^2/4) makes scatter
    essentially zero unless aspect << arcsin(sqrt(4*ln(100))/(Dk*L_para)).

    Parameters
    ----------
    theta_rad : grazing angle from horizontal [rad]
    dip_deg   : magnetic dip angle [deg] (positive downward)

    Returns
    -------
    aspect_deg [deg]  (0 = perpendicular to B, optimal scatter)
    """
    theta_deg  = float(np.degrees(theta_rad))
    aspect_deg = abs(90.0 - float(dip_deg) - theta_deg)
    return float(np.clip(aspect_deg, 0.0, 90.0))


def es_scatter_valid(aspect_deg: float,
                     L_para_m:   float,
                     freq_MHz:   float,
                     threshold:  float = 0.01) -> bool:
    """
    Check whether the Booker-Gordon aspect condition allows significant scatter.

    Returns True when exp(-Dk_para^2*L_para^2/4) > threshold, meaning the
    field-aligned scale L_para is large enough relative to the Bragg half-wavelength
    that scatter is not suppressed below threshold level.

    Parameters
    ----------
    aspect_deg : aspect angle [deg]
    L_para_m   : field-aligned irregularity scale [m]
    freq_MHz   : operating frequency [MHz]
    threshold  : minimum Gaussian suppression factor (default 0.01)
    """
    k0      = 2.0 * np.pi * freq_MHz * 1e6 / C_MS
    Dk      = 2.0 * k0
    Dk_para = Dk * np.sin(np.radians(float(aspect_deg)))
    gauss   = np.exp(-Dk_para**2 * float(L_para_m)**2 / 4.0)
    return bool(gauss > float(threshold))
