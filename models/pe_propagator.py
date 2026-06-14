"""
Module M4 – Parabolic-equation (PE) propagator for plasma-bubble scattering.
Algorithm: Split-Step Fourier (SSF), wide-angle PE.

Key functions / class:
    PEPropagator.propagate()    – run the SSF loop
    PEPropagator.analyze()      – extract AOA spectrum + scatter modes
    construct_incident_field()  – RT->PE Gaussian-beam interface

Per-step SSF (Carrano 2020):
    Step A  (refraction, spatial):
        u_half = u * exp[i k0 (n - 1) dx]
    Step B  (diffraction, spectral):
        U_half = FFT(u_half)
        kx_eff = sqrt(k0^2 - kz^2) - k0        (wide-angle propagator)
        U_next = U_half * exp[i kx_eff dx]
        u_next = IFFT(U_next)

Implemented in Part 5.
"""
import numpy as np
from config import C_MS, C_KMS, FREQ_MHZ, PE
from utils import freq_to_k0


class PEPropagator:
    """
    2-D wide-angle PE / SSF propagator.

    Parameters
    ----------
    freq_MHz  : operating frequency [MHz]
    pe_params : PE config dict (from config.PE)
    """

    def __init__(self,
                 freq_MHz: float = FREQ_MHZ,
                 pe_params: dict = PE):
        self.freq    = freq_MHz
        self.k0      = freq_to_k0(freq_MHz)     # [km^-1]
        self.params  = pe_params

    # ── Domain helpers ────────────────────────────────────────────────────────

    def extract_domain(self,
                       Ne_2d: np.ndarray,
                       x_array: np.ndarray,
                       z_array: np.ndarray,
                       x_range: tuple,
                       z_range: tuple
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Extract the plasma-bubble sub-domain from the background grid
        and interpolate to the fine PE grid (dz = pe_params['dz_m'] / 1000 km).

        Parameters
        ----------
        Ne_2d   : (Nx, Nz) background density [m^-3]
        x_array : (Nx,) [km]
        z_array : (Nz,) [km]
        x_range : (x_min, x_max) [km]
        z_range : (z_min, z_max) [km]

        Returns
        -------
        n_pe   : (Nx_pe, Nz_pe) refractive-index on fine PE grid
        x_pe   : (Nx_pe,) [km]
        z_pe   : (Nz_pe,) [km]
        """
        raise NotImplementedError("Implemented in Part 5")

    # ── Core SSF step ─────────────────────────────────────────────────────────

    @staticmethod
    def ssf_step(u: np.ndarray,
                 n_half: np.ndarray,
                 k0: float,
                 dz_km: float,
                 dx_km: float) -> np.ndarray:
        """
        One SSF step: u(x) -> u(x + dx).

        Parameters
        ----------
        u      : (Nz,) complex field at current x slice
        n_half : (Nz,) refractive index at x + dx/2
        k0     : free-space wavenumber [km^-1]
        dz_km  : vertical sampling [km]
        dx_km  : step size [km]

        Returns u_next : (Nz,) complex field at x + dx.
        """
        raise NotImplementedError("Implemented in Part 5")

    @staticmethod
    def apply_pml(u: np.ndarray, n_pml: int, sigma: float) -> np.ndarray:
        """
        Apply exponential-taper PML at both ends of the z axis.
        Prevents artificial reflections from domain edges.
        """
        raise NotImplementedError("Implemented in Part 5")

    # ── Main propagation loop ─────────────────────────────────────────────────

    def propagate(self,
                  u_init: np.ndarray,
                  n_field: np.ndarray,
                  dx_km: float,
                  dz_km: float) -> tuple[np.ndarray, np.ndarray | None]:
        """
        Propagate u_init through n_field with SSF.

        Parameters
        ----------
        u_init  : (Nz,) complex initial field at x = x_entry
        n_field : (Nx, Nz) refractive-index field of PE domain
        dx_km   : propagation step [km]
        dz_km   : vertical sampling [km]

        Returns
        -------
        u_out     : (Nz,) field at x = x_exit
        u_history : (Nx, Nz) or None  (only when store_history = True)
        """
        raise NotImplementedError("Implemented in Part 5")

    # ── Output analysis ───────────────────────────────────────────────────────

    def analyze(self,
                u_out: np.ndarray,
                z_array_km: np.ndarray,
                dx_total_km: float) -> dict:
        """
        Compute AOA spectrum and scatter parameters from exit field.

        Returns
        -------
        {
          'aoa_deg'          : (Nz,) arrival-angle array [deg],
          'power_aoa'        : (Nz,) power angular spectrum [|U|^2],
          'mean_aoa_deg'     : dominant arrival angle [deg],
          'delta_tau_ms'     : RMS delay spread [ms],
          'tau_extra_mean_ms': mean extra delay vs free-space [ms],
        }
        """
        raise NotImplementedError("Implemented in Part 5")

    def extract_scatter_modes(self,
                               aoa_deg: np.ndarray,
                               power_aoa: np.ndarray,
                               aoa_inc_deg: float) -> list[dict]:
        """
        Identify distinct scatter peaks in the AOA spectrum.
        Each peak -> one scatter mode dict:
            {'aoa_deg', 'delta_aoa_deg', 'power'}
        Sorted by descending power.
        """
        raise NotImplementedError("Implemented in Part 5")


# ── RT -> PE interface ─────────────────────────────────────────────────────────

def construct_incident_field(A_inc: float,
                              beta_inc_deg: float,
                              z_inc_km: float,
                              z_array_km: np.ndarray,
                              k0_per_km: float,
                              w0_km: float = PE['w0_km']
                              ) -> np.ndarray:
    """
    Build the PE initial field at the bubble entry plane (x = x_entry).
    Uses a Gaussian beam centred at z_inc_km with tilt angle beta_inc_deg.

    u(x0, z) = A_inc * exp[-(z - z_inc)^2/(2w0^2)] * exp[-i k0 sin(beta) * z]

    Parameters
    ----------
    A_inc        : field amplitude (derived from RT power)
    beta_inc_deg : ray elevation at bubble entry [deg]
    z_inc_km     : ray centre height at entry [km]
    z_array_km   : PE z-axis grid [km]
    k0_per_km    : free-space wavenumber [km^-1]
    w0_km        : Gaussian beam waist [km]

    Returns u_init : (Nz,) complex initial field.
    """
    raise NotImplementedError("Implemented in Part 5")
