"""
Module M1 - Ionospheric electron-density field.
Layers (added in order):
    1. IRI-2016 background  (iri2016)         [Phase 5: optional lateral sampling]
    2. TID perturbation      (Hooke 1968)      [Phase 5: multi-direction superposition]
    3. Es thin layer         (Hao 2017)
    4. Plasma bubble         (analytic Gaussian depletion)
Output: Ne_2d [m^-3] on the 2-D (x, z) Cartesian grid.
"""
import numpy as np
import iri2016
from scipy.interpolate import RegularGridInterpolator
from config import (K_FP, IRI_DT, IRI_LAT, IRI_LON,
                    TID, ES, BUBBLE, SPREAD_F, FREQ_MHZ,
                    TX_LAT, TX_LON, LINK_BEARING_DEG, IRI_LATERAL)
from utils import ne_to_n2, n2_to_n


class IonosphereModel:
    """
    Build the 2-D electron-density field Ne(x, z) by stacking background
    IRI + optional TID / Es / plasma-bubble perturbations.

    Parameters
    ----------
    iri_params         : dict   keys: 'dt', 'lat', 'lon'
    tid_params         : dict   from config.TID
    es_params          : dict   from config.ES
    bubble_params      : dict   from config.BUBBLE
    spread_f_params    : dict   from config.SPREAD_F
    lateral_iri_params : dict   Phase 5 lateral sampling; keys: 'enable', 'spacing_km',
                                'tx_lat', 'tx_lon', 'bearing_deg'
    freq_MHz           : float  operating frequency
    """

    def __init__(self,
                 iri_params:         dict | None = None,
                 tid_params:         dict | None = None,
                 es_params:          dict | None = None,
                 bubble_params:      dict | None = None,
                 spread_f_params:    dict | None = None,
                 lateral_iri_params: dict | None = None,
                 freq_MHz:           float = FREQ_MHZ):
        self.iri         = iri_params      or {'dt': IRI_DT, 'lat': IRI_LAT, 'lon': IRI_LON}
        self.tid         = tid_params      or TID
        self.es          = es_params       or ES
        self.bubble      = bubble_params   or BUBBLE
        self.spread_f    = spread_f_params or SPREAD_F
        self.lateral_iri = lateral_iri_params or IRI_LATERAL
        self.freq        = freq_MHz

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
        # 1. IRI background (single profile broadcast, or lateral 2-D sampling)
        if self.lateral_iri.get('enable', False):
            Ne_2d = self._iri_lateral_background(x_array, z_array)
        else:
            Ne_1d = self._iri_background(z_array)
            Ne_2d = np.tile(Ne_1d, (len(x_array), 1))

        # 2. TID (single component or multi-direction superposition)
        if self.tid.get('enable', False):
            Ne_2d = self._add_tid(Ne_2d, x_array, z_array, t)

        # 3. Es – uniform in x; delta relative to homogeneous IRI background
        if self.es.get('enable', False):
            Ne_1d_ref = Ne_2d[len(x_array) // 2, :]   # mid-path column as ref
            Ne_1d_es  = self._add_es_layer(Ne_1d_ref.copy(), z_array)
            Ne_2d    += (Ne_1d_es - Ne_1d_ref)[np.newaxis, :]

        # 4. Plasma bubble (Gaussian depletion)
        if self.bubble.get('enable', False):
            Ne_2d = self._add_plasma_bubble(Ne_2d, x_array, z_array)

        # 5. Spread-F phase screen (applied last)
        if self.spread_f.get('enable', False):
            Ne_2d = self._add_spread_f(Ne_2d, x_array, z_array)

        Ne_2d = np.maximum(Ne_2d, 0.0)
        n2_2d = ne_to_n2(Ne_2d, self.freq)
        n_2d  = n2_to_n(n2_2d)
        return Ne_2d, n_2d

    # ── Private helpers ───────────────────────────────────────────────────────

    def _iri_at_point(self,
                      lat: float, lon: float,
                      z_array: np.ndarray, dt) -> np.ndarray:
        """Single iri2016 query at (lat, lon); return Ne_1d(z) [m^-3]."""
        z_min = float(z_array[0])
        z_max = float(z_array[-1])
        dz    = float(round(z_array[1] - z_array[0], 6))
        res   = iri2016.IRI(dt, (z_min, z_max, dz), lat, lon)
        Ne_1d = res['ne'].values.astype(float)
        assert len(Ne_1d) == len(z_array), (
            f"iri2016 returned {len(Ne_1d)} pts but z_array has {len(z_array)} pts"
        )
        return np.maximum(Ne_1d, 0.0)

    def _iri_background(self, z_array: np.ndarray) -> np.ndarray:
        """Single mid-path IRI profile (backward-compatible path)."""
        dt  = self.iri.get('dt',  IRI_DT)
        lat = self.iri.get('lat', IRI_LAT)
        lon = self.iri.get('lon', IRI_LON)
        return self._iri_at_point(lat, lon, z_array, dt)

    def _iri_lateral_background(self,
                                 x_array: np.ndarray,
                                 z_array: np.ndarray) -> np.ndarray:
        """
        Sample IRI along the link great-circle at ~spacing_km intervals.

        Geographic coordinates at each sample x are derived from the TX
        position and link bearing via the spherical-Earth destination formula.
        Profiles are linearly interpolated onto the full x_array.

        Ref: Cervera & Harris (2014) use position-dependent IRI for 3-D tracing.
        """
        from utils import destination_point

        spacing   = float(self.lateral_iri.get('spacing_km', 50.0))
        tx_lat    = float(self.lateral_iri.get('tx_lat',      TX_LAT))
        tx_lon    = float(self.lateral_iri.get('tx_lon',      TX_LON))
        bearing   = float(self.lateral_iri.get('bearing_deg', LINK_BEARING_DEG))
        dt        = self.iri.get('dt', IRI_DT)

        x_min, x_max = float(x_array[0]), float(x_array[-1])
        # Extend sample range slightly beyond grid to avoid edge extrapolation
        x_samp = np.arange(x_min, x_max + spacing * 0.5, spacing)
        # Ensure endpoints are included
        if x_samp[-1] < x_max:
            x_samp = np.append(x_samp, x_max)
        x_samp = np.clip(x_samp, x_min, x_max)
        x_samp = np.unique(x_samp)

        Ne_rows = []
        for xs in x_samp:
            lat_s, lon_s = destination_point(tx_lat, tx_lon, bearing, float(xs))
            Ne_rows.append(self._iri_at_point(lat_s, lon_s, z_array, dt))

        Ne_samp = np.array(Ne_rows)  # (N_samp, Nz)

        if len(x_samp) == 1:
            return np.tile(Ne_samp[0], (len(x_array), 1))

        # Bilinear interpolation onto full x_array
        interp = RegularGridInterpolator(
            (x_samp, z_array), Ne_samp,
            method='linear', bounds_error=False, fill_value=None
        )
        pts   = np.array([[x, z] for x in x_array for z in z_array])
        Ne_2d = np.maximum(
            interp(pts).reshape(len(x_array), len(z_array)), 0.0
        )
        return Ne_2d

    def _tid_one_component(self,
                            Ne_2d:       np.ndarray,
                            x_array:     np.ndarray,
                            z_array:     np.ndarray,
                            t:           float,
                            amplitude:   float,
                            lambda_h_km: float,
                            T_s:         float,
                            I_dip_deg:   float,
                            H_km:        float,
                            omega_b:     float,
                            kx:          float | None = None
                            ) -> np.ndarray:
        """
        Apply one TID component to Ne_2d using the Hooke (1968) formula.

        kx [rad/m]: horizontal wave-number projected onto the 2-D plane (x-axis).
            If None, defaults to 2pi/lambda_h (TID aligned with link).
            If given, must be the projected value; lambda_h_km is still used
            for the dispersion relation (to get kz from the full |k_h|).

        Hooke (1968) Eq. for AGW-driven Ne perturbation (Koval 2018 Eq. 1):
            dNe = A * (-dNe0/dz * sinI * sin(phase)
                       - k_para * Ne0 * cos(phase))
            k_para = kx * cos(I) + kz * sin(I)
            kz^2   = kh^2 * (wb^2/w^2 - 1) - 1/(4H^2)

        The normalisation ensures max|dNe/Ne0| = amplitude over the
        significant-density region (Ne > 1% of NmF2), following the
        fix from Koval et al. (2018).
        """
        lam_h_m = lambda_h_km * 1e3              # [m]
        kh_full = 2.0 * np.pi / lam_h_m          # full horizontal |k| for dispersion
        if kx is None:
            kx = kh_full                          # TID aligned with link

        I_rad   = np.deg2rad(I_dip_deg)
        H_m     = H_km * 1e3
        omega   = 2.0 * np.pi / T_s

        kz2 = kh_full**2 * (omega_b**2 / omega**2 - 1.0) - 1.0 / (4.0 * H_m**2)
        if kz2 <= 0.0:
            return Ne_2d                          # evanescent mode – skip

        kz     = np.sqrt(kz2)
        sinI   = np.sin(I_rad)
        k_para = kx * np.cos(I_rad) + kz * np.sin(I_rad)

        X_m = x_array * 1e3                       # (Nx,) [m]
        Z_m = z_array * 1e3                        # (Nz,) [m]
        phase = (kx * X_m[:, np.newaxis]
                 + kz * Z_m[np.newaxis, :]
                 - omega * t)                      # (Nx, Nz)

        # Vertical gradient of background Ne at each column (handles lateral IRI)
        dz_m      = (z_array[1] - z_array[0]) * 1e3
        dNe_dz_2d = np.gradient(Ne_2d, dz_m, axis=1)   # (Nx, Nz)

        dNe_raw = (-dNe_dz_2d * sinI * np.sin(phase)
                   - k_para * Ne_2d * np.cos(phase))    # (Nx, Nz)

        # Normalise over significant-density region (Ne > 1% of peak)
        Ne_peak = Ne_2d.max()
        f_mask  = Ne_2d > Ne_peak * 0.01
        if not f_mask.any():
            return Ne_2d
        max_ratio = np.max(np.abs(dNe_raw[f_mask]) / Ne_2d[f_mask])
        if max_ratio < 1e-30:
            return Ne_2d
        dNe = dNe_raw * (amplitude / max_ratio)

        return np.maximum(Ne_2d + dNe, 0.0)

    def _add_tid(self,
                 Ne_2d:   np.ndarray,
                 x_array: np.ndarray,
                 z_array: np.ndarray,
                 t:       float) -> np.ndarray:
        """
        Dispatcher: single-component (backward compat) or
        multi-direction superposition (Phase 5).
        """
        p      = self.tid
        n_comp = int(p.get('n_components', 1))
        omega_b = p.get('omega_b_rad_s', 2.0 * np.pi / 1200.0)

        if n_comp == 1:
            return self._tid_one_component(
                Ne_2d, x_array, z_array, t,
                amplitude   = p['amplitude'],
                lambda_h_km = p['lambda_h_km'],
                T_s         = p['T_s'],
                I_dip_deg   = p['I_dip_deg'],
                H_km        = p['H_km'],
                omega_b     = omega_b,
            )

        # Multi-component: each component projected onto link direction (2-D plane)
        az_list  = p.get('az_deg_list',      [0.0]           * n_comp)
        amp_list = p.get('amplitude_list',   [p['amplitude']] * n_comp)
        T_list   = p.get('period_s_list',    [p['T_s']]       * n_comp)
        lam_list = p.get('lambda_h_km_list', [p['lambda_h_km']] * n_comp)
        link_az  = float(p.get('link_bearing_deg', LINK_BEARING_DEG))

        Ne_out = Ne_2d.copy()
        for i in range(n_comp):
            # Project TID wave vector onto the link (x-axis) direction
            # kx_eff = k_h * cos(az_TID - az_link)
            # Using kh_full for dispersion; projected kx for phase and k_para
            az_rel   = np.radians(float(az_list[i]) - link_az)
            cos_proj = np.cos(az_rel)
            if abs(cos_proj) < 1e-6:
                # TID nearly perpendicular to link: no x-variation in 2-D cross-section
                continue
            kx_i = (2.0 * np.pi / (float(lam_list[i]) * 1e3)) * cos_proj

            Ne_out = self._tid_one_component(
                Ne_out, x_array, z_array, t,
                amplitude   = float(amp_list[i]),
                lambda_h_km = float(lam_list[i]),
                T_s         = float(T_list[i]),
                I_dip_deg   = p['I_dip_deg'],
                H_km        = p['H_km'],
                omega_b     = omega_b,
                kx          = kx_i,
            )
        return Ne_out

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
