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
      'tau_ms'         : group delay (single path) [ms],
      'delta_tau_ms'   : delay spread from scattering [ms],
    }

Implemented in Part 4.
"""
import numpy as np
from math import factorial
from config import C_MS, C_KMS, ES


class EsLayerModel:
    """
    Es-layer three-segment propagation model.

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

    # ── Public interface ──────────────────────────────────────────────────────

    def classify(self, f_MHz: float) -> tuple[str, float]:
        """
        Decide propagation regime from foEs/f.

        Returns
        -------
        mode  : 'reflect' | 'scatter' | 'mixed'
        alpha : reflection weight  (1.0 for reflect, 0.0 for scatter)
        """
        raise NotImplementedError("Implemented in Part 4")

    def reflection_coeff_sq(self, theta_rad: float, f_MHz: float) -> float:
        """
        |rho|^2 – squared reflection coefficient from the Es thin layer.
        Hao (2017) Eq. 3 with n = n_exp = 5.

        Parameters
        ----------
        theta_rad : grazing angle (complement of incidence angle) [rad]
        f_MHz     : operating frequency [MHz]
        """
        raise NotImplementedError("Implemented in Part 4")

    def scatter_cross_section(self, theta_rad: float, f_MHz: float) -> float:
        """
        sigma_5 – Es irregularity scatter cross section [m^-1].
        Hao (2017) Eq. 7, anisotropic irregularity model.

        Parameters
        ----------
        theta_rad : grazing angle [rad]
        f_MHz     : operating frequency [MHz]
        """
        raise NotImplementedError("Implemented in Part 4")

    def compute_power(self,
                      Pt_W: float,
                      Gt: float,
                      Gr: float,
                      f_MHz: float,
                      D_km: float,
                      theta_rad: float) -> dict:
        """
        Combined three-segment power calculation.

        Parameters
        ----------
        Pt_W      : transmit power [W]
        Gt, Gr    : TX / RX antenna gains (linear, isotropic = 1)
        f_MHz     : operating frequency [MHz]
        D_km      : one-way path length TX->Es->RX [km]
        theta_rad : Es grazing angle [rad]

        Returns a result dict (see module docstring).
        """
        raise NotImplementedError("Implemented in Part 4")

    # ── Boundary-condition helper for PE interface ────────────────────────────

    def transmission_amplitude(self, theta_rad: float, f_MHz: float) -> float:
        """
        Amplitude transmission coefficient T = sqrt(1 - |rho|^2).
        Used when PE propagates through Es height as a thin-layer BC.
        """
        raise NotImplementedError("Implemented in Part 4")
