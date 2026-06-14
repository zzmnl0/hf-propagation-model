"""
Module M2 (basic) – 2-D ray tracer.
Integrates the Haselgrove equations (isotropic, Cartesian 2-D)
with an RK4 stepper; supports a fan of rays from a single TX point.

Each ray is returned as a dict:
    {
      'beta_deg'        : launch elevation [deg],
      'trajectory'      : list of (x, z, kx, kz) state vectors,
      'group_path_km'   : accumulated group path [km],
      'tau_ms'          : group delay [ms],
      'h_reflect_km'    : maximum height reached [km],
      'at_Es'           : dict with params at Es height, or None,
      'at_bubble'       : dict with params at bubble bottom, or None,
    }

Implemented in Part 2.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from config import C_KMS, FREQ_MHZ, RT, ES, BUBBLE
from utils import freq_to_k0


class RefractiveIndex:
    """
    Interpolable refractive-index model backed by a 2-D Ne grid.

    Parameters
    ----------
    Ne_2d   : (Nx, Nz) electron density [m^-3]
    x_array : (Nx,)   [km]
    z_array : (Nz,)   [km]
    freq_MHz: float
    """

    def __init__(self,
                 Ne_2d: np.ndarray,
                 x_array: np.ndarray,
                 z_array: np.ndarray,
                 freq_MHz: float = FREQ_MHZ):
        raise NotImplementedError("Implemented in Part 2")

    def n2(self, x: float, z: float) -> float:
        """n^2(x, z)  – interpolated from grid."""
        raise NotImplementedError("Implemented in Part 2")

    def n(self, x: float, z: float) -> float:
        """n(x, z)."""
        raise NotImplementedError("Implemented in Part 2")

    def grad_n2(self, x: float, z: float,
                dx: float = 0.5, dz: float = 0.5
                ) -> tuple[float, float]:
        """
        Numerical gradient (dn^2/dx, dn^2/dz) via central difference.
        Step sizes dx, dz in [km].
        """
        raise NotImplementedError("Implemented in Part 2")


def haselgrove_deriv(state: np.ndarray,
                     n_model: RefractiveIndex,
                     k0: float) -> np.ndarray:
    """
    Right-hand side of 2-D Haselgrove equations.

    State vector: [x, z, kx, kz]  (positions in km, wavenumbers in km^-1).
    Returns d(state)/dP where P is group path.

    Equations (isotropic, 2-D Cartesian):
        dx/dP  =  kx / |k|
        dz/dP  =  kz / |k|
        dkx/dP = (k0^2/2|k|) * dn^2/dx
        dkz/dP = (k0^2/2|k|) * dn^2/dz
    """
    raise NotImplementedError("Implemented in Part 2")


def rk4_step(state: np.ndarray,
             n_model: RefractiveIndex,
             k0: float,
             ds: float) -> np.ndarray:
    """Single RK4 step of length ds [km]. Returns new state."""
    raise NotImplementedError("Implemented in Part 2")


def trace_single_ray(tx_pos: tuple,
                     beta_deg: float,
                     n_model: RefractiveIndex,
                     freq_MHz: float = FREQ_MHZ,
                     rt_params: dict = RT,
                     h_Es_km: float | None = None,
                     h_bubble_bot_km: float | None = None
                     ) -> dict:
    """
    Trace one ray from tx_pos at elevation beta_deg.

    Parameters
    ----------
    tx_pos          : (x0, z0) transmitter position [km, km]
    beta_deg        : launch elevation above horizon [deg]
    n_model         : RefractiveIndex instance
    freq_MHz        : operating frequency
    rt_params       : RT config dict
    h_Es_km         : record ray state when first crossing this height (Es)
    h_bubble_bot_km : record ray state when first crossing this height (bubble)

    Returns a ray dict (see module docstring).
    """
    raise NotImplementedError("Implemented in Part 2")


def shoot_rays_fan(tx_pos: tuple,
                   n_model: RefractiveIndex,
                   freq_MHz: float = FREQ_MHZ,
                   rt_params: dict = RT,
                   h_Es_km: float | None = None,
                   h_bubble_bot_km: float | None = None
                   ) -> list[dict]:
    """
    Launch a fan of rays spanning beta_min..beta_max with n_fan rays.

    Returns list of ray dicts (see module docstring).
    """
    raise NotImplementedError("Implemented in Part 2")
