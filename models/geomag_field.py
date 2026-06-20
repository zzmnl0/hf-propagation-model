"""
Module: 2-D geomagnetic field grid via ppigrf IGRF-14.

GeomagnField2D pre-computes fH_MHz(x, z) and dip_deg(x, z) on the background
grid by converting x-coordinates to geographic (lat, lon) along the link
great-circle and querying ppigrf at each (lat, lon, alt) point.

Phase 5 note: altitude variation of B is included (full 2-D grid). The
'fixed-angle approximation' (alpha = pi/2 - dip) in RefractiveIndexAH is
preserved; only the spatial variation of fH and dip is upgraded.

Reference: IGRF-14 via ppigrf (Finlay et al. 2010, updated 2024).
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator


class GeomagnField2D:
    """
    Position-dependent geomagnetic parameters on a 2-D (x, z) grid.

    Parameters
    ----------
    x_km        : (Nx,) horizontal distances from TX along link [km]
    z_km        : (Nz,) heights above Earth surface [km]
    tx_lat      : TX latitude [deg N]
    tx_lon      : TX longitude [deg E]
    bearing_deg : great-circle bearing from TX to RX [deg, CW from North]
    dt          : datetime for IGRF epoch
    """

    def __init__(self,
                 x_km:        np.ndarray,
                 z_km:        np.ndarray,
                 tx_lat:      float,
                 tx_lon:      float,
                 bearing_deg: float,
                 dt):
        import ppigrf
        from utils import destination_point

        Nx, Nz = len(x_km), len(z_km)

        # Geographic coordinates for each x position along link
        lats = np.array([destination_point(tx_lat, tx_lon, bearing_deg, float(x))[0]
                         for x in x_km])  # (Nx,)
        lons = np.array([destination_point(tx_lat, tx_lon, bearing_deg, float(x))[1]
                         for x in x_km])  # (Nx,)

        # Build flat (Nx*Nz,) arrays for ppigrf vectorized call
        lat_flat = np.repeat(lats, Nz)                      # each lat repeated Nz times
        lon_flat = np.repeat(lons, Nz)
        alt_flat = np.tile(z_km, Nx)                         # z_km repeated Nx times

        Be, Bn, Bu = ppigrf.igrf(lon_flat, lat_flat, alt_flat, dt)
        Be = np.asarray(Be).ravel()
        Bn = np.asarray(Bn).ravel()
        Bu = np.asarray(Bu).ravel()

        F       = np.sqrt(Be**2 + Bn**2 + Bu**2)            # total field [nT]
        fH_flat = 2.7994e-5 * F                              # gyrofrequency [MHz]
        # Dip angle: positive downward (magnetic inclination, NED convention)
        dip_flat = np.degrees(np.arctan2(-Bu, np.sqrt(Be**2 + Bn**2)))

        fH_2d  = fH_flat.reshape(Nx, Nz)
        dip_2d = dip_flat.reshape(Nx, Nz)

        kw = dict(method='linear', bounds_error=False, fill_value=None)
        self._fH_interp  = RegularGridInterpolator((x_km, z_km), fH_2d,  **kw)
        self._dip_interp = RegularGridInterpolator((x_km, z_km), dip_2d, **kw)

    def fH_MHz_batch(self, pts: np.ndarray) -> np.ndarray:
        """(M, 2) array of (x_km, z_km) -> (M,) fH_MHz."""
        return self._fH_interp(pts)

    def dip_deg_batch(self, pts: np.ndarray) -> np.ndarray:
        """(M, 2) array of (x_km, z_km) -> (M,) dip angle [deg, positive down]."""
        return self._dip_interp(pts)
