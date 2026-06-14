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
from config import C_MS, C_KMS, ES


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
