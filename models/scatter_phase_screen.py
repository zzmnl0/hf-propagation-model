"""
models/scatter_phase_screen.py
Rino (1979) 2-D power-law phase screen and Knepp (1983) MPS propagator.

7d (C3): RinoPhaseScreen2D  -- 2-D anisotropic phase screen drawn from Rino spectrum
7e (C5): MultiPhaseScreen   -- Knepp (1983) multi-plane PE propagator

Physics summary
---------------
Rino (1979) 2-D phase PSD (height-integrated from 3-D Born):
  Phi_phi(kx,ky) = (r_e*lambda)^2 * C_s_abs * DeltaL * F(p)
                   * (kx^2 + ky^2 + kappa_0^2)^(-(p+2)/2)

  F(p) = sqrt(pi) * Gamma((p+2)/2) / Gamma((p+3)/2)   [Rino 1979 Eq. 14]
  C_s_abs = Cs_rel^2 * Ne^2 * p * (2pi)^2 * kappa_0^p  [scatter_born convention]

Phase variance (closed form):
  sigma_phi^2 = (r_e*lambda)^2 * C_s_abs * DeltaL * F(p) * kappa_0^(-p) / (2*pi*p)

Spectral synthesis (discrete):
  F[m,n] = (N/dx) * sqrt(Phi_phi(kx_m, ky_n)) * xi[m,n]    xi ~ CN(0,1)
  phi(x,y) = IFFT2(F)   [numpy convention: (1/N^2)*sum F*exp]

Knepp (1983) MPS parabolic equation (Eq. 2-4):
  Phase screen:       u -> u * exp(i*phi_i)
  Free-space step:    U -> U * exp(-i*(kx^2+ky^2)*dz/(2*k0))

References
----------
Rino, C.L. (1979). A power law phase screen model for ionospheric scintillation.
  Part 1: Weak scatter. Radio Sci. 14(6), 1135-1145.
Rino, C.L. (1979). Part 2: Strong scatter. Radio Sci. 14(6), 1147-1155.
Knepp, D.L. (1983). Multiple phase-screen calculation of the temporal behavior
  of stochastic waves. Proc. IEEE, 71(6), 722-737.
"""
import numpy as np
from scipy.special import gamma as gamma_fn
from config import C_MS, FREQ_MHZ, SCATTER, PHASE_SCREEN


# ── Module constants ──────────────────────────────────────────────────────────
R_E_M  = 2.818e-15   # classical electron radius [m]
TWO_PI = 2.0 * np.pi


# ── 7d: Rino (1979) 2-D phase screen ─────────────────────────────────────────

class RinoPhaseScreen2D:
    """
    2-D random phase screen with Rino (1979) power-law spectrum.

    The 2-D phase PSD is derived by integrating the 3-D Born volume spectrum
    along the propagation direction (height integration):

      Phi_phi(kx, ky) = (r_e*lambda)^2 * C_s_abs * DeltaL * F(p)
                        * (kx^2 + ky^2 + kappa_0^2)^(-(p+2)/2)

    Screen realisations are generated via spectral synthesis (IFFT of
    sqrt(PSD) * complex Gaussian noise), correctly normalised so that the
    measured phase variance equals the theoretical sigma_phi^2.

    Parameters
    ----------
    freq_MHz   : operating frequency [MHz]
    Cs_rel     : relative rms Ne fluctuation delta_Ne/Ne (0-1)
    p          : power-law spectral index (F-layer typical: 3.0)
    L_outer_km : outer irregularity scale [km]
    Ne_m3      : background electron density [m^-3]
    DeltaL_km  : scatter layer thickness [km]
    n_pts      : transverse grid size (square: n_pts x n_pts)
    dx_m       : transverse grid spacing [m]
    seed       : random seed (int or None)
    """

    def __init__(self,
                 freq_MHz:   float = FREQ_MHZ,
                 Cs_rel:     float = SCATTER['Cs_rel'],
                 p:          float = SCATTER['p'],
                 L_outer_km: float = SCATTER['L_outer_km'],
                 Ne_m3:      float = 1e11,
                 DeltaL_km:  float = SCATTER['DeltaL_km'],
                 n_pts:      int   = PHASE_SCREEN['n_pts'],
                 dx_m:       float = PHASE_SCREEN['dx_m'],
                 seed              = PHASE_SCREEN['seed']):
        self.freq_MHz   = float(freq_MHz)
        self.Cs_rel     = float(Cs_rel)
        self.p          = float(p)
        self.L_outer_km = float(L_outer_km)
        self.Ne_m3      = float(Ne_m3)
        self.DeltaL_km  = float(DeltaL_km)
        self.n_pts      = int(n_pts)
        self.dx_m       = float(dx_m)
        self.seed       = seed

        lam_m       = C_MS / (freq_MHz * 1e6)
        self.lam_m  = lam_m
        self.k0     = TWO_PI / lam_m          # free-space wavenumber [rad/m]

    # ── Spectral helpers ──────────────────────────────────────────────────────

    def _c_s_abs(self) -> float:
        """
        Absolute spectral coefficient C_s_abs [m^(-p-6)].
        Convention from scatter_born.py (outer-scale dominated):
          C_s_abs = Cs_rel^2 * Ne^2 * p * (2pi)^2 * kappa_0^p
        """
        kappa_0 = TWO_PI / (self.L_outer_km * 1e3)
        return (self.Cs_rel * self.Ne_m3)**2 * self.p * TWO_PI**2 * kappa_0**self.p

    def _F_factor(self) -> float:
        """F(p) = sqrt(pi)*Gamma((p+2)/2)/Gamma((p+3)/2). Rino (1979) Eq. 14."""
        return (np.sqrt(np.pi)
                * gamma_fn(0.5 * (self.p + 2.0))
                / gamma_fn(0.5 * (self.p + 3.0)))

    def phi_2d_spectrum(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        """
        2-D phase PSD Phi_phi(kx, ky) [rad^2 * m^2].
        kx, ky in [rad/m]; arrays must be broadcastable.
        """
        kappa_0 = TWO_PI / (self.L_outer_km * 1e3)
        DeltaL  = self.DeltaL_km * 1e3
        A       = (R_E_M * self.lam_m)**2 * self._c_s_abs() * DeltaL * self._F_factor()
        denom   = (kx**2 + ky**2 + kappa_0**2) ** (0.5 * (self.p + 2.0))
        return A / denom

    def sigma_phi_rms(self) -> float:
        """
        Theoretical RMS phase variance [rad] from the closed-form integral.

          sigma_phi^2 = A * kappa_0^(-p) / (2*pi*p)

          where A = (r_e*lambda)^2 * C_s_abs * DeltaL * F(p)
        """
        kappa_0 = TWO_PI / (self.L_outer_km * 1e3)
        DeltaL  = self.DeltaL_km * 1e3
        A       = (R_E_M * self.lam_m)**2 * self._c_s_abs() * DeltaL * self._F_factor()
        sigma_sq = A * kappa_0**(-self.p) / (TWO_PI * self.p)
        return float(np.sqrt(max(sigma_sq, 0.0)))

    # ── Screen generation ─────────────────────────────────────────────────────

    def generate(self) -> np.ndarray:
        """
        Draw one phase screen realisation phi(x, y) [rad].

        Spectral synthesis with IFFT:
          F[m,n] = (N / dx) * sqrt(Phi(kx_m, ky_n)) * xi[m,n]
          phi(x,y) = Re{ IFFT2(F) }

        Variance: E[|phi|^2] = sum Phi * dk^2 / (N*dx)^2 = sigma_phi^2.
        Returns (n_pts, n_pts) float array [rad].
        """
        N  = self.n_pts
        dx = self.dx_m

        kx1d   = np.fft.fftfreq(N, d=dx) * TWO_PI   # [rad/m]
        ky1d   = np.fft.fftfreq(N, d=dx) * TWO_PI
        KX, KY = np.meshgrid(kx1d, ky1d, indexing='ij')

        S     = self.phi_2d_spectrum(KX, KY)           # (N,N) [rad^2 * m^2]
        amp   = (N / dx) * np.sqrt(np.maximum(S, 0.0)) # normalised amplitude

        rng   = np.random.default_rng(self.seed)
        xi    = ((rng.standard_normal((N, N))
                  + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0))
        phi_k  = amp * xi
        phi_k[0, 0] = 0.0                              # enforce zero mean

        return np.real(np.fft.ifft2(phi_k))            # (N,N) [rad]


# ── 7e: Knepp (1983) Multi-Phase-Screen propagator ───────────────────────────

class MultiPhaseScreen:
    """
    Knepp (1983) multi-phase-screen (MPS) parabolic-equation propagator.

    Splits the propagation path through a scattering volume into N_screen
    thin phase screens with free-space diffraction between them.

    Each screen step (Knepp 1983, Eq. 2-4):
      1. Apply phase screen:   u -> u * exp(i * phi_i(x,y))
      2. Angular spectrum FFT: U = FFT2(u)
      3. Free-space step:      U -> U * exp(-i*(kx^2+ky^2)*dz/(2*k0))
      4. Back to spatial:      u = IFFT2(U)

    The S4 scintillation index is computed from the output intensity statistics:
      S4 = std(I) / mean(I)

    Parameters
    ----------
    freq_MHz   : operating frequency [MHz]
    Cs_rel     : relative rms Ne fluctuation
    p          : spectral index
    L_outer_km : outer scale [km]
    Ne_m3      : background electron density [m^-3]
    DeltaL_km  : total scatter layer thickness [km]
    N_screen   : number of phase screens
    n_pts      : transverse grid size
    dx_m       : transverse spacing [m]
    seed       : random seed base (None = random; screen i uses seed+i)
    """

    def __init__(self,
                 freq_MHz:   float = FREQ_MHZ,
                 Cs_rel:     float = SCATTER['Cs_rel'],
                 p:          float = SCATTER['p'],
                 L_outer_km: float = SCATTER['L_outer_km'],
                 Ne_m3:      float = 1e11,
                 DeltaL_km:  float = SCATTER['DeltaL_km'],
                 N_screen:   int   = PHASE_SCREEN['N_screen'],
                 n_pts:      int   = PHASE_SCREEN['n_pts'],
                 dx_m:       float = PHASE_SCREEN['dx_m'],
                 seed              = PHASE_SCREEN['seed']):
        self.N_screen  = int(N_screen)
        self.n_pts     = int(n_pts)
        self.dx_m      = float(dx_m)
        lam_m          = C_MS / (float(freq_MHz) * 1e6)
        self.k0        = TWO_PI / lam_m
        self.dz_m      = float(DeltaL_km) * 1e3 / max(N_screen, 1)

        self.screens = [
            RinoPhaseScreen2D(
                freq_MHz=freq_MHz, Cs_rel=Cs_rel, p=p,
                L_outer_km=L_outer_km, Ne_m3=Ne_m3,
                DeltaL_km=float(DeltaL_km) / max(N_screen, 1),
                n_pts=n_pts, dx_m=dx_m,
                seed=None if seed is None else int(seed) + i)
            for i in range(N_screen)
        ]

    def propagate(self, u0: np.ndarray | None = None) -> dict:
        """
        Propagate a field through N_screen phase screens.

        Parameters
        ----------
        u0 : (n_pts, n_pts) complex array   Input field (plane wave if None).

        Returns
        -------
        dict:
          u_out   : output field (n_pts, n_pts) complex
          I_mean  : mean output intensity (normalised to input mean)
          S4      : sqrt(var(I)) / mean(I)  scintillation index
          Pr_norm : mean(I_out) / mean(I_in)  (energy conservation check, ~1)
        """
        N  = self.n_pts
        dx = self.dx_m

        if u0 is None:
            u = np.ones((N, N), dtype=complex)
        else:
            u = u0.astype(complex).copy()

        I0 = float(np.mean(np.abs(u)**2))

        kx1d   = np.fft.fftfreq(N, d=dx) * TWO_PI
        ky1d   = np.fft.fftfreq(N, d=dx) * TWO_PI
        KX, KY = np.meshgrid(kx1d, ky1d, indexing='ij')
        H_free = np.exp(-1j * (KX**2 + KY**2) * self.dz_m / (2.0 * self.k0))

        for scr in self.screens:
            phi = scr.generate()
            u   = u * np.exp(1j * phi)
            U   = np.fft.fft2(u)
            U   = U * H_free
            u   = np.fft.ifft2(U)

        I_out  = np.abs(u)**2
        I_mean = float(np.mean(I_out))
        S4     = float(np.std(I_out) / max(I_mean, 1e-30))
        Pr_norm = I_mean / max(I0, 1e-30)

        return {
            'u_out'   : u,
            'I_mean'  : I_mean,
            'S4'      : S4,
            'Pr_norm' : Pr_norm,
        }


def _delta_kernel(tau_axis: np.ndarray, tau0: float) -> np.ndarray:
    out = np.zeros_like(tau_axis, dtype=float)
    if out.size:
        out[int(np.argmin(np.abs(tau_axis - float(tau0))))] = 1.0
    return out


def _gaussian_kernel(tau_axis: np.ndarray, tau0: float, sigma_ms: float) -> np.ndarray:
    sigma = max(float(sigma_ms), 1e-9)
    k = np.exp(-0.5 * ((tau_axis - float(tau0)) / sigma) ** 2)
    s = float(np.sum(k))
    return _delta_kernel(tau_axis, tau0) if s <= 0.0 else k / s


def make_path_mps_delay_kernel(mode: dict,
                               n3d,
                               tau_axis: np.ndarray,
                               scatter_params: dict | None = None,
                               phase_params: dict | None = None,
                               freq_MHz: float = FREQ_MHZ,
                               t: float = 0.0) -> dict:
    """Build a path-related MPS delay kernel from 3-D ray samples.

    This is a stage-P3 MPS statistical kernel, not a complete 3-D full-wave
    solution. It uses Rino phase-screen statistics and Knepp MPS propagation
    to obtain S4/energy diagnostics, then maps phase variance to delay spread.
    """
    sp = {**SCATTER, **(scatter_params or {})}
    pp = {**PHASE_SCREEN, **(phase_params or {})}
    tau0 = float(mode.get("tau_ms", 0.0))
    sigma_phi = max(float(sp.get("sigma_phi_rad", sp.get("sigma_phi", 0.0))), 0.0)
    if sigma_phi == 0.0 and float(sp.get("Cs_rel", 0.0)) <= 0.0:
        kernel = _delta_kernel(tau_axis, tau0)
        return {
            "scatter_method": "MPS",
            "mps_status": "sigma_phi_zero_delta",
            "mps_delay_kernel": kernel,
            "kernel_normalization": float(np.sum(kernel)),
            "sigma_phi_rad": 0.0,
            "sigma_phi_sq": 0.0,
            "sigma_tau_ms": 0.0,
            "S4": 0.0,
            "Pr_norm": 1.0,
            "energy_conservation_error": 0.0,
            "scatter_solver_dimensionality": "path_mps_statistical",
            "scatter_validity": "phase_screen_not_fullwave",
        }

    ecef_path = mode.get("trajectory_ecef_km")
    pts = np.asarray(ecef_path if ecef_path is not None else [], dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 2:
        kernel = _delta_kernel(tau_axis, tau0)
        return {
            "scatter_method": "MPS_required",
            "mps_status": "empty_path",
            "mps_delay_kernel": kernel,
            "kernel_normalization": float(np.sum(kernel)),
            "scatter_validity": "no_path_no_mps",
        }

    sample_step = max(int(pts.shape[0] // 24), 1)
    samples = pts[::sample_step, :3]
    if not hasattr(n3d, "get_Ne_ecef"):
        raise TypeError("formal MPS path sampling requires an ECEF NeProvider")
    radius = float(mode["earth_radius_km"])
    ne_vals = [max(float(n3d.get_Ne_ecef(p, radius, t=t)), 0.0)
               for p in samples]
    Ne_eff = float(np.median(ne_vals)) if ne_vals else 1.0e11
    DeltaL_km = float(sp.get("DeltaL_km", mode.get("scatter_DeltaL_km", 100.0)))
    mps = MultiPhaseScreen(
        freq_MHz=freq_MHz,
        Cs_rel=float(sp.get("Cs_rel", 0.01)),
        p=float(sp.get("p", 3.0)),
        L_outer_km=float(sp.get("L_outer_km", 10.0)),
        Ne_m3=Ne_eff,
        DeltaL_km=DeltaL_km,
        N_screen=int(pp.get("N_screen", 10)),
        n_pts=int(pp.get("n_pts", 128)),
        dx_m=float(pp.get("dx_m", 200.0)),
        seed=pp.get("seed", None),
    )
    res = mps.propagate()
    screen_sigma = mps.screens[0].sigma_phi_rms() if mps.screens else 0.0
    sigma_phi = max(sigma_phi, screen_sigma * np.sqrt(max(int(pp.get("N_screen", 10)), 1)))
    sigma_tau = float(mode.get("sigma_tau_ms", 0.0))
    if sigma_tau <= 0.0:
        sigma_tau = max(sigma_phi, float(res.get("S4", 0.0))) * DeltaL_km / 2.998e5 * 1e3
    kernel = _gaussian_kernel(tau_axis, tau0, max(sigma_tau, float(np.median(np.diff(tau_axis))) if len(tau_axis) > 1 else 0.02))
    pr_norm = float(res.get("Pr_norm", 1.0))
    return {
        "scatter_method": "MPS",
        "mps_status": "ok",
        "mps_delay_kernel": kernel,
        "kernel_normalization": float(np.sum(kernel)),
        "sigma_phi_rad": float(sigma_phi),
        "sigma_phi_sq": float(sigma_phi * sigma_phi),
        "sigma_tau_ms": float(sigma_tau),
        "S4": float(res.get("S4", 0.0)),
        "Pr_norm": pr_norm,
        "energy_conservation_error": abs(pr_norm - 1.0),
        "scatter_solver_dimensionality": "path_mps_statistical",
        "scatter_validity": "phase_screen_not_fullwave",
        "scatter_path_coordinate_model": "spherical_ECEF",
        "theory_reference": "Rino 1979 phase screen; Knepp 1983 multiple phase screen",
    }
