"""
Module M4 – Parabolic-equation (PE) propagator for plasma-bubble scattering.
Algorithm: Split-Step Fourier (SSF), wide-angle PE (Carrano 2020).

Per-step SSF:
    Step A  (refraction, spatial):
        u_half = u * exp[i k0 (n - 1) dx]
    Step B  (diffraction, spectral):
        U_half = FFT(u_half)
        kx_eff = sqrt(k0^2 - kz^2) - k0        (wide-angle propagator)
        U_next = U_half * exp[i kx_eff dx]
        u_next = IFFT(U_next)
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import find_peaks
from config import C_MS, C_KMS, FREQ_MHZ, PE, RE_KM
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
        self.freq   = freq_MHz
        self.k0     = freq_to_k0(freq_MHz)    # [km^-1]
        self.params = pe_params

    # ── Domain extraction ─────────────────────────────────────────────────────

    def extract_domain(self,
                       Ne_2d:    np.ndarray,
                       x_array:  np.ndarray,
                       z_array:  np.ndarray,
                       x_range:  tuple,
                       z_range:  tuple
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Extract the plasma-bubble sub-domain from the background Ne grid
        and interpolate to the fine PE grid (dz = pe_params['dz_m']/1000 km).

        Parameters
        ----------
        Ne_2d   : (Nx, Nz) background electron density [m^-3]
        x_array : (Nx,) [km]
        z_array : (Nz,) [km]
        x_range : (x_min, x_max) [km]
        z_range : (z_min, z_max) [km]

        Returns
        -------
        n_pe : (Nx_pe, Nz_pe) refractive-index field on fine PE grid
        x_pe : (Nx_pe,) [km]
        z_pe : (Nz_pe,) [km]
        """
        dx_km   = self.params['dx_km']
        dz_km   = self.params['dz_m'] / 1000.0    # m -> km
        freq_Hz = self.freq * 1e6

        # Refractive index on background grid
        fp2  = (8.98 ** 2) * Ne_2d
        n2   = np.maximum(1.0 - fp2 / freq_Hz ** 2, 1e-6)
        n_bg = np.sqrt(n2)

        interp = RegularGridInterpolator(
            (x_array, z_array), n_bg,
            method='linear', bounds_error=False, fill_value=1.0
        )

        x_pe = np.arange(x_range[0], x_range[1] + dx_km * 0.5, dx_km)
        z_pe = np.arange(z_range[0], z_range[1] + dz_km * 0.5, dz_km)

        XX, ZZ = np.meshgrid(x_pe, z_pe, indexing='ij')
        pts    = np.column_stack([XX.ravel(), ZZ.ravel()])
        n_pe   = interp(pts).reshape(len(x_pe), len(z_pe))

        if self.params.get('earth_flat', True):
            # Earth-flattening: n_eff(x,z) = n(x,z)*(1 + z/R_E)
            # Accounts for coordinate curvature in Cartesian PE frame (Part 7, Ch.8 A2).
            n_pe = np.maximum(n_pe * (1.0 + ZZ / RE_KM), 1e-6)

        return n_pe, x_pe, z_pe

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
        # Step A: refraction (spatial domain)
        u_half = u * np.exp(1j * k0 * (n_half - 1.0) * dx_km)

        # Step B: diffraction (spectral domain, wide-angle)
        Nz     = len(u_half)
        U_half = np.fft.fft(u_half)
        kz     = np.fft.fftfreq(Nz, d=dz_km) * 2.0 * np.pi   # [km^-1]
        kx_eff = np.sqrt(np.maximum(k0 ** 2 - kz ** 2, 0.0)) - k0
        U_next = U_half * np.exp(1j * kx_eff * dx_km)
        return np.fft.ifft(U_next)

    @staticmethod
    def apply_pml(u: np.ndarray,
                  n_pml: int,
                  sigma: float) -> np.ndarray:
        """
        Apply exponential-taper PML at both ends of the z axis.
        Prevents artificial reflections from domain edges.

        n_pml : number of grid points in the absorbing layer
        sigma : maximum damping coefficient (dimensionless, 0.3~1.0)
        """
        N      = len(u)
        n_act  = min(n_pml, N // 2)
        window = np.ones(N, dtype=float)
        for i in range(n_act):
            frac          = (n_act - i) / n_act
            att           = np.exp(-sigma * frac ** 2)
            window[i]    *= att
            window[N-1-i] *= att
        return u * window

    # ── Main propagation loop ─────────────────────────────────────────────────

    def propagate(self,
                  u_init:  np.ndarray,
                  n_field: np.ndarray,
                  dx_km:   float,
                  dz_km:   float
                  ) -> tuple[np.ndarray, np.ndarray | None]:
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
        u_history : (Nx-1, Nz) or None  (only when store_history = True)
        """
        n_pml   = self.params['n_pml']
        sigma   = self.params['sigma_pml']
        store   = self.params['store_history']

        Nx, _Nz = n_field.shape
        u        = u_init.astype(complex).copy()
        history  = []

        for ix in range(Nx - 1):
            n_half = 0.5 * (n_field[ix] + n_field[ix + 1])
            u      = self.ssf_step(u, n_half, self.k0, dz_km, dx_km)
            u      = self.apply_pml(u, n_pml, sigma)
            if store:
                history.append(u.copy())

        return u, (np.array(history) if store else None)

    # ── Output analysis ───────────────────────────────────────────────────────

    def analyze(self,
                u_out:       np.ndarray,
                z_array_km:  np.ndarray,
                dx_total_km: float) -> dict:
        """
        Compute AOA spectrum and scatter parameters from exit field.

        Returns
        -------
        {
          'aoa_deg'          : (Nz,) arrival-angle array [deg]  (NaN for evanescent),
          'power_aoa'        : (Nz,) power angular spectrum [|U|^2],
          'mean_aoa_deg'     : dominant arrival angle [deg],
          'delta_tau_ms'     : RMS delay spread [ms],
          'tau_extra_mean_ms': mean extra delay vs free-space [ms],
        }
        """
        Nz    = len(u_out)
        dz_km = float(z_array_km[1] - z_array_km[0])

        U_out  = np.fft.fft(u_out)
        kz     = np.fft.fftfreq(Nz, d=dz_km) * 2.0 * np.pi   # [km^-1]

        valid   = np.abs(kz) <= self.k0
        kz_safe = np.where(valid, kz, 0.0)
        aoa_deg = np.where(valid,
                           np.degrees(np.arcsin(kz_safe / self.k0)),
                           np.nan)

        power_aoa         = np.abs(U_out) ** 2
        power_aoa[~valid] = 0.0

        # Extra path for each angular component: dx / cos(aoa) - dx
        aoa_rad_safe = np.radians(np.where(valid, aoa_deg, 0.0))
        cos_aoa      = np.where(valid, np.cos(aoa_rad_safe), 1.0)
        cos_aoa      = np.clip(cos_aoa, 1e-6, None)
        delta_path   = np.where(valid,
                                dx_total_km / cos_aoa - dx_total_km,
                                0.0)
        tau_extra_ms = delta_path / C_KMS * 1e3

        P_total      = float(np.sum(power_aoa)) + 1e-30
        tau_mean_ms  = float(np.nansum(power_aoa * tau_extra_ms) / P_total)
        tau_sq_mean  = float(np.nansum(power_aoa * tau_extra_ms ** 2) / P_total)
        delta_tau_ms = float(np.sqrt(max(tau_sq_mean - tau_mean_ms ** 2, 0.0)))

        peak_idx = int(np.argmax(power_aoa))
        mean_aoa = float(aoa_deg[peak_idx]) if valid[peak_idx] else 0.0

        return {
            'aoa_deg'           : aoa_deg,
            'power_aoa'         : power_aoa,
            'mean_aoa_deg'      : mean_aoa,
            'delta_tau_ms'      : delta_tau_ms,
            'tau_extra_mean_ms' : tau_mean_ms,
        }

    def extract_scatter_modes(self,
                               aoa_deg:    np.ndarray,
                               power_aoa:  np.ndarray,
                               aoa_inc_deg: float) -> list[dict]:
        """
        Identify distinct scatter peaks in the AOA spectrum.
        Each peak -> one mode dict: {'aoa_deg', 'delta_aoa_deg', 'power'}.
        Sorted by descending power.
        """
        min_frac   = self.params['min_power_frac']
        peak_power = float(np.nanmax(power_aoa)) if np.any(power_aoa > 0) else 0.0
        threshold  = peak_power * min_frac

        valid_power = np.where(np.isnan(aoa_deg), 0.0, power_aoa)
        peaks, _    = find_peaks(valid_power, height=threshold, distance=5)

        modes = []
        for pk in peaks:
            if np.isnan(aoa_deg[pk]):
                continue
            modes.append({
                'aoa_deg'       : float(aoa_deg[pk]),
                'delta_aoa_deg' : float(aoa_deg[pk]) - aoa_inc_deg,
                'power'         : float(power_aoa[pk]),
            })
        return sorted(modes, key=lambda m: -m['power'])


# ── RT -> PE interface ─────────────────────────────────────────────────────────

def construct_incident_field(A_inc:        float,
                              beta_inc_deg: float,
                              z_inc_km:     float,
                              z_array_km:   np.ndarray,
                              k0_per_km:    float,
                              w0_km:        float = PE['w0_km']
                              ) -> np.ndarray:
    """
    Build the PE initial field at the bubble entry plane (x = x_entry).
    Uses a Gaussian beam centred at z_inc_km tilted at beta_inc_deg.

        u(x0, z) = A_inc * exp[-(z-z_inc)^2/(2w0^2)] * exp[-i k0 sin(beta) * (z-z_inc)]

    Parameters
    ----------
    A_inc        : field amplitude
    beta_inc_deg : ray elevation at bubble entry [deg]
    z_inc_km     : ray centre height at entry [km]
    z_array_km   : PE z-axis grid [km]
    k0_per_km    : free-space wavenumber [km^-1]
    w0_km        : Gaussian beam waist [km]

    Returns u_init : (Nz,) complex initial field.
    """
    beta_rad = np.radians(beta_inc_deg)
    dz       = z_array_km - z_inc_km
    envelope = np.exp(-dz ** 2 / (2.0 * w0_km ** 2))
    # Positive sign: upward ray (beta>0) has kz = k0*sin(beta) > 0 in PE envelope
    # Full field E = u*exp(ik0*x); for kz>0 plane wave: u ~ exp(+ikz*z)
    phase    = k0_per_km * np.sin(beta_rad) * dz
    return (A_inc * envelope * np.exp(1j * phase)).astype(complex)
