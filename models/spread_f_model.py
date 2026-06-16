"""
SpreadFModel - Thin random phase screen for spread-F Ne irregularities.

Physical basis: Rino (1979) power-law phase screen model.
Reference: Rino, C.L. (1979). A power law phase screen model for ionospheric
           scintillation. Radio Science, 14(6), 1135-1145.

The model adds spatially-structured Ne perturbations at h_screen_km using a
1-D power-law random field, modulated by a Gaussian vertical envelope.
"""
import numpy as np


class SpreadFModel:
    """
    Rino (1979) power-law thin phase screen for spread-F irregularities.

    Parameters
    ----------
    Cs           : phase spectral strength (Rino 1979 definition)
    p            : power-law spectral index (typical 2.5 to 4.0)
    h_screen_km  : screen centre height [km]
    L0_km        : outer irregularity scale [km]
    seed         : random seed (int or None) for reproducibility
    """

    def __init__(self, Cs: float = 1e-3, p: float = 3.0,
                 h_screen_km: float = 300.0, L0_km: float = 50.0,
                 seed=None):
        self.Cs          = float(Cs)
        self.p           = float(p)
        self.h_screen_km = float(h_screen_km)
        self.L0_km       = float(L0_km)
        self.seed        = seed

    def apply(self,
              Ne_2d:  np.ndarray,
              x_km:   np.ndarray,
              z_km:   np.ndarray) -> np.ndarray:
        """
        Add spread-F Ne perturbation to Ne_2d.

        dNe(x, z) = Cs * phi(x) * Ne_bg(z) * G(z - h_screen)

        where phi(x) is a unit-RMS power-law random field and G is a
        Gaussian vertical envelope (half-width 20 km).

        Returns updated Ne_2d with negative values clamped to zero.
        """
        if self.Cs <= 0.0:
            return Ne_2d

        phi = self._gen_screen(x_km)                          # (Nx,) unit-RMS

        h_w = 20.0                                            # Gaussian half-width [km]
        G_z = np.exp(-0.5 * ((z_km - self.h_screen_km) / h_w) ** 2)  # (Nz,)

        Ne_bg = Ne_2d.mean(axis=0)                            # (Nz,) horizontal mean
        dNe   = (self.Cs
                 * phi[:, np.newaxis]
                 * Ne_bg[np.newaxis, :]
                 * G_z[np.newaxis, :])

        return np.maximum(Ne_2d + dNe, 0.0)

    def _gen_screen(self, x_km: np.ndarray) -> np.ndarray:
        """
        Generate 1-D power-law random field along x (unit RMS).

        Power spectrum: S(k) = (k^2 + k0^2)^(-(p+1)/2)
        where k0 = 2pi/L0 is the outer-scale wavenumber.
        """
        N  = len(x_km)
        dx = float(x_km[1] - x_km[0])

        k_cyc = np.fft.rfftfreq(N, d=dx)          # cycles/km
        k     = 2.0 * np.pi * k_cyc               # rad/km
        k0    = 2.0 * np.pi / self.L0_km          # outer-scale wavenumber

        S_k = np.where(k > 0,
                       (k ** 2 + k0 ** 2) ** (-(self.p + 1) / 2.0),
                       0.0)

        rng   = np.random.default_rng(self.seed)
        n_pos = len(k)
        amp   = np.sqrt(S_k / (N * dx))
        xi    = rng.standard_normal(n_pos) + 1j * rng.standard_normal(n_pos)
        phi_k = amp * xi
        phi_x = np.fft.irfft(phi_k, n=N)

        rms = float(np.std(phi_x))
        if rms > 1e-30:
            phi_x = phi_x / rms
        return phi_x
