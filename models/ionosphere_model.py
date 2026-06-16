"""
Module M1 – Ionospheric electron-density field.
Layers (added in order):
    1. IRI-2016 background  (iri2016)
    2. TID perturbation      (Hooke 1968 / Koval 2018)
    3. Es thin layer         (Hao 2017)
    4. Plasma bubble         (analytic Gaussian depletion)
Output: Ne_2d [m^-3] on the 2-D (x, z) Cartesian grid.
"""
import numpy as np
import iri2016
from config import (K_FP, IRI_DT, IRI_LAT, IRI_LON,
                    TID, ES, BUBBLE, SPREAD_F, FREQ_MHZ)
from utils import ne_to_n2, n2_to_n


class IonosphereModel:
    """
    Build the 2-D electron-density field Ne(x, z) by stacking background
    IRI + optional TID / Es / plasma-bubble perturbations.

    Parameters
    ----------
    iri_params    : dict   keys: 'dt', 'lat', 'lon'
    tid_params    : dict   from config.TID
    es_params     : dict   from config.ES
    bubble_params : dict   from config.BUBBLE
    freq_MHz      : float  operating frequency (needed for n^2 computation)
    """

    def __init__(self,
                 iri_params:      dict | None = None,
                 tid_params:      dict | None = None,
                 es_params:       dict | None = None,
                 bubble_params:   dict | None = None,
                 spread_f_params: dict | None = None,
                 freq_MHz:        float = FREQ_MHZ):
        self.iri      = iri_params      or {'dt': IRI_DT, 'lat': IRI_LAT, 'lon': IRI_LON}
        self.tid      = tid_params      or TID
        self.es       = es_params       or ES
        self.bubble   = bubble_params   or BUBBLE
        self.spread_f = spread_f_params or SPREAD_F
        self.freq     = freq_MHz

    # ── Public interface ──────────────────────────────────────────────────────

    def build_Ne_field(self,
                       x_array: np.ndarray,
                       z_array: np.ndarray,
                       t: float = 0.0
                       ) -> tuple[np.ndarray, np.ndarray]:
        """
        Build the electron-density field on a 2-D grid.

        Parameters
        ----------
        x_array : (Nx,) horizontal distance from TX [km]
        z_array : (Nz,) height above Earth surface  [km]
        t       : snapshot time [s]  (used by TID only)

        Returns
        -------
        Ne_2d : (Nx, Nz) electron density [m^-3]
        n_2d  : (Nx, Nz) refractive index n (isotropic, real)
        """
        # 1. IRI 1D background
        Ne_1d = self._iri_background(z_array)

        # 2. Broadcast to (Nx, Nz)
        Ne_2d = np.tile(Ne_1d, (len(x_array), 1))

        # 3. TID – uses pure IRI as background (x-dependent)
        if self.tid.get('enable', False):
            Ne_2d = self._add_tid(Ne_2d, x_array, z_array, t)

        # 4. Es – uniform in x; compute contribution from IRI-only 1D
        if self.es.get('enable', False):
            Ne_1d_es = self._add_es_layer(Ne_1d.copy(), z_array)
            Ne_2d   += (Ne_1d_es - Ne_1d)[np.newaxis, :]

        # 5. Plasma bubble (Gaussian depletion)
        if self.bubble.get('enable', False):
            Ne_2d = self._add_plasma_bubble(Ne_2d, x_array, z_array)

        # 6. Spread-F phase screen (applied last)
        if self.spread_f.get('enable', False):
            Ne_2d = self._add_spread_f(Ne_2d, x_array, z_array)

        Ne_2d = np.maximum(Ne_2d, 0.0)
        n2_2d = ne_to_n2(Ne_2d, self.freq)
        n_2d  = n2_to_n(n2_2d)
        return Ne_2d, n_2d

    # ── Private helpers ───────────────────────────────────────────────────────

    def _iri_background(self, z_array: np.ndarray) -> np.ndarray:
        """
        Call iri2016 and return 1-D Ne profile Ne_1d(z) [m^-3].
        The same profile is broadcast horizontally (homogeneous background).
        """
        dt  = self.iri.get('dt',  IRI_DT)
        lat = self.iri.get('lat', IRI_LAT)
        lon = self.iri.get('lon', IRI_LON)
        z_min = float(z_array[0])
        z_max = float(z_array[-1])
        dz    = float(round(z_array[1] - z_array[0], 6))
        res   = iri2016.IRI(dt, (z_min, z_max, dz), lat, lon)
        Ne_1d = res['ne'].values.astype(float)   # [m^-3]
        assert len(Ne_1d) == len(z_array), (
            f"iri2016 returned {len(Ne_1d)} pts but z_array has {len(z_array)} pts "
            f"(z_min={z_min}, z_max={z_max}, dz={dz})"
        )
        return np.maximum(Ne_1d, 0.0)

    def _add_tid(self,
                 Ne_2d:   np.ndarray,
                 x_array: np.ndarray,
                 z_array: np.ndarray,
                 t:       float) -> np.ndarray:
        """
        Add TID perturbation dNe(x, z, t) to Ne_2d.
        Uses Hooke (1968) formula; see Koval et al. (2018) Eq. 1.
        Returns updated Ne_2d.
        """
        p       = self.tid
        lam_h_m = p['lambda_h_km'] * 1e3           # horizontal wavelength [m]
        T_s     = p['T_s']                          # period [s]
        amp     = p['amplitude']                    # desired max |dNe/Ne0|
        I_rad   = np.deg2rad(p['I_dip_deg'])
        H_m     = p['H_km'] * 1e3                  # scale height [m]
        omega_b = p.get('omega_b_rad_s', 2.0 * np.pi / 1200.0)

        kx    = 2.0 * np.pi / lam_h_m              # [rad/m]
        omega = 2.0 * np.pi / T_s                  # [rad/s]

        kz2 = kx**2 * (omega_b**2 / omega**2 - 1.0) - 1.0 / (4.0 * H_m**2)
        if kz2 <= 0.0:
            return Ne_2d                            # evanescent – skip TID

        kz     = np.sqrt(kz2)                       # [rad/m]
        k_para = kx * np.cos(I_rad) + kz * np.sin(I_rad)
        sinI   = np.sin(I_rad)

        # Background 1D (all columns identical at this stage)
        Ne_bg  = Ne_2d[0, :].copy()                # (Nz,)
        dz_m   = (z_array[1] - z_array[0]) * 1e3  # [m]
        dNe_dz = np.gradient(Ne_bg, dz_m)          # (Nz,) [m^-3/m]

        # Phase field (Nx, Nz)
        X_m   = x_array * 1e3                      # (Nx,) [m]
        Z_m   = z_array * 1e3                      # (Nz,) [m]
        phase = (kx * X_m[:, np.newaxis]
                 + kz * Z_m[np.newaxis, :]
                 - omega * t)                       # (Nx, Nz)

        # Hooke (1968) un-normalised dNe
        dNe_raw = (-dNe_dz[np.newaxis, :] * sinI * np.sin(phase)
                   - k_para * Ne_bg[np.newaxis, :] * np.cos(phase))  # (Nx, Nz)

        # Normalise so max|dNe/Ne0| = amp, evaluated only over significant-
        # density altitudes (> 1% of F2 peak). Without this restriction the
        # near-zero D-layer values dominate the ratio and suppress F-region
        # perturbation to effectively zero.
        Ne_peak = Ne_bg.max()
        f_mask  = Ne_bg > Ne_peak * 0.01          # (Nz,) significant-Ne mask
        if not f_mask.any():
            return Ne_2d
        max_ratio = np.max(
            np.abs(dNe_raw[:, f_mask]) / Ne_bg[np.newaxis, f_mask]
        )
        if max_ratio < 1e-30:
            return Ne_2d
        dNe = dNe_raw * (amp / max_ratio)

        return np.maximum(Ne_2d + dNe, 0.0)

    def _add_es_layer(self,
                      Ne_1d:   np.ndarray,
                      z_array: np.ndarray) -> np.ndarray:
        """
        Add Es thin-layer peak to a 1-D Ne profile.
        Profile shape: Ne_Es(z') = Nmax*[1-(z'/dh)^(2n)], |z'| <= dh
        Returns updated Ne_1d.
        """
        p         = self.es
        foEs_MHz  = p['foEs_MHz']
        h_Es_km   = p['h_Es_km']
        delta_h_m = p['delta_h_m']                 # half-thickness [m]
        n_exp     = int(p['n_exp'])

        Nmax      = (foEs_MHz * 1e6 / K_FP) ** 2   # [m^-3]

        z_prime_m = (z_array - h_Es_km) * 1e3      # distance from Es centre [m]
        mask      = np.abs(z_prime_m) <= delta_h_m
        Ne_es     = np.zeros_like(Ne_1d)
        zp        = z_prime_m[mask]
        Ne_es[mask] = Nmax * (1.0 - (zp / delta_h_m) ** (2 * n_exp))

        return Ne_1d + Ne_es

    def _add_plasma_bubble(self,
                           Ne_2d:   np.ndarray,
                           x_array: np.ndarray,
                           z_array: np.ndarray) -> np.ndarray:
        """
        Subtract Gaussian density depletion (plasma bubble).
        dNe = Ne_bg * delta_max * exp(-Dx^2/Lx^2 - Dz^2/Lz^2)
        Returns updated Ne_2d (negative values clamped to 0).
        """
        p         = self.bubble
        delta_max = p['delta_max']
        x0, z0    = p['x0_km'],  p['z0_km']
        Lx, Lz    = p['Lx_km'],  p['Lz_km']

        dx = (x_array - x0)[:, np.newaxis]         # (Nx, 1)
        dz = (z_array - z0)[np.newaxis, :]         # (1, Nz)

        depletion = delta_max * np.exp(-(dx / Lx) ** 2 - (dz / Lz) ** 2)
        return np.maximum(Ne_2d * (1.0 - depletion), 0.0)

    def _add_spread_f(self,
                      Ne_2d:   np.ndarray,
                      x_array: np.ndarray,
                      z_array: np.ndarray) -> np.ndarray:
        """
        Add spread-F Ne irregularities via Rino (1979) power-law phase screen.
        Returns updated Ne_2d.
        """
        from .spread_f_model import SpreadFModel
        p   = self.spread_f
        sfm = SpreadFModel(
            Cs         = float(p.get('Cs',          1e-3)),
            p          = float(p.get('p',           3.0)),
            h_screen_km= float(p.get('h_screen_km', 300.0)),
            L0_km      = float(p.get('L0_km',       50.0)),
            seed       = p.get('seed', None),
        )
        return sfm.apply(Ne_2d, x_array, z_array)
