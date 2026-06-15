"""
Module M5-M7 – Hybrid propagation model  (top-level integration).
Orchestrates the full pipeline:

  M1 (ionosphere) -> M2 (RT / P2P) -> M3 (Es) -> M4 (PE)
  -> M5 (synthesis) -> M6 (dual-path) -> M7 (main-mode ID)

Primary entry point:
    HybridPropagationModel.compute(tx, rx) -> (mode_results, tau, pd, main_mode)

mode_results element keys:
    label, tau_ms, delta_tau_ms, Pr_W, Pr_dBW,
    h_reflect_km, group_path_km, beta_deg, phi_deg
"""
import numpy as np
from scipy.signal import find_peaks
from config import (FREQ_MHZ, PT_W, GT, GR, ES, BUBBLE,
                    BG_X, BG_Z, PE, P2P, MODE, TX_POS, RX_POS, C_MS, EPS0)
from .ionosphere_model import IonosphereModel
from .ray_tracer import RefractiveIndex
from .point_to_point import (find_all_rays_p2p, classify_mode,
                              extract_es_params, extract_bubble_entry)
from .es_model import EsLayerModel
from .pe_propagator import PEPropagator, construct_incident_field
from utils import free_space_loss_dB, to_dBW, freq_to_k0


# ── Field-amplitude helper ────────────────────────────────────────────────────

def _amplitude_from_power(Pr_W: float) -> float:
    """Approximate field amplitude [V/m] from received power [W]."""
    return float(np.sqrt(2.0 * max(Pr_W, 0.0) / (EPS0 * C_MS)))


# ── P-D spectrum ──────────────────────────────────────────────────────────────

def build_pd_spectrum(mode_results: list,
                      tau_axis: np.ndarray | None = None,
                      tau_res_ms: float = 0.02
                      ) -> tuple[np.ndarray, np.ndarray]:
    """
    Construct the power-delay (P-D) spectrum by summing Gaussian-spread
    contributions from each propagation mode.

    Each mode i contributes:
        Pr_i * exp(-(tau - tau_i)^2 / (2 * sigma_i^2))

    Parameters
    ----------
    mode_results : list of mode dicts (must contain 'tau_ms', 'Pr_W', 'delta_tau_ms')
    tau_axis     : optional pre-defined delay axis [ms]; auto-generated if None
    tau_res_ms   : delay-axis resolution [ms]

    Returns (tau_axis, pd_W) both as 1-D arrays.
    """
    if not mode_results:
        return np.array([0.0]), np.array([0.0])

    taus = [m['tau_ms'] for m in mode_results]
    if tau_axis is None:
        tau_min  = max(0.0, min(taus) - 2.0)
        tau_max  = max(taus) + 2.0
        tau_axis = np.arange(tau_min, tau_max, tau_res_ms)

    pd = np.zeros(len(tau_axis))
    for m in mode_results:
        tau_c = m['tau_ms']
        Pr    = float(m.get('Pr_W', 0.0))
        sigma = max(float(m.get('delta_tau_ms', tau_res_ms)), tau_res_ms)
        pd   += Pr * np.exp(-(tau_axis - tau_c) ** 2 / (2.0 * sigma ** 2))

    return tau_axis, pd


# ── Main-mode identification ──────────────────────────────────────────────────

def identify_main_mode(tau_axis: np.ndarray,
                        pd_W: np.ndarray,
                        mode_results: list
                        ) -> tuple[dict | None, list]:
    """
    Detect peaks in the P-D spectrum and rank them by power.

    Returns (main_mode_dict, all_modes_ranked) where main_mode_dict is the
    mode whose tau is closest to the highest P-D peak.
    Returns (None, []) if no peaks found.
    """
    if len(pd_W) == 0 or not mode_results:
        return None, []

    peak_max  = float(np.max(pd_W))
    threshold = max(peak_max * 0.01, 1e-30)
    peaks, _  = find_peaks(pd_W, height=threshold, distance=5)

    if len(peaks) == 0:
        # fall back: use global maximum
        peaks = np.array([int(np.argmax(pd_W))])

    peak_powers = pd_W[peaks]
    sorted_idx  = np.argsort(-peak_powers)

    def _nearest(pk_idx):
        tau_pk = tau_axis[pk_idx]
        return min(mode_results, key=lambda m: abs(m['tau_ms'] - tau_pk))

    ranked = [_nearest(peaks[i]) for i in sorted_idx]
    return ranked[0], ranked


# ── Main class ────────────────────────────────────────────────────────────────

class HybridPropagationModel:
    """
    End-to-end hybrid HF propagation model.

    Parameters
    ----------
    iono_params   : dict with optional keys:
                      'iri_params'    -> dict {dt, lat, lon}
                      'tid_params'    -> config.TID-style dict
                      'es_params'     -> config.ES-style dict  (+ 'enable' key)
                      'bubble_params' -> config.BUBBLE-style dict (+ 'enable' key)
    radar_params  : dict with keys: Pt_W, Gt, Gr, freq_MHz
    x_array, z_array : background grid [km]
    """

    def __init__(self,
                 iono_params:  dict,
                 radar_params: dict,
                 x_array: np.ndarray = BG_X,
                 z_array: np.ndarray = BG_Z):
        self.iono_params  = iono_params
        self.radar_params = radar_params
        self.x_array      = x_array
        self.z_array      = z_array
        self.freq         = float(radar_params.get('freq_MHz', FREQ_MHZ))
        self.k0           = freq_to_k0(self.freq)

        # Ionosphere model (passes iri/tid/es/bubble sub-dicts)
        self.iono = IonosphereModel(freq_MHz=self.freq, **_split_iono(iono_params))

        # Es model (optional)
        es_cfg = iono_params.get('es_params') or {}
        if es_cfg.get('enable', False):
            kw       = {k: v for k, v in es_cfg.items() if k != 'enable'}
            self.es  = EsLayerModel(**kw)
        else:
            self.es  = None

        # PE propagator
        self.pe = PEPropagator(freq_MHz=self.freq)

    def compute(self,
                tx_km: tuple = TX_POS,
                rx_km: tuple = RX_POS,
                t: float = 0.0,
                p2p_params: dict | None = None
                ) -> tuple[list, np.ndarray, np.ndarray, dict | None]:
        """
        Run the full hybrid pipeline for a single TX-RX link.

        Pipeline
        --------
        1.  Build electron-density field  Ne(x, z)
        2.  Solve P2P ray paths  (variational, Nosikov 2020)
        3.  For each ray:
              3a.  Apply Es model       if ray crosses h_Es_km
              3b.  Apply PE/SSF         if ray enters plasma bubble
        4.  Build P-D spectrum
        5.  Identify dominant mode

        Parameters
        ----------
        p2p_params : override P2P config (e.g. reduce n_init for tests)

        Returns
        -------
        mode_results : list of mode dicts
        tau_axis     : delay axis [ms]
        pd_W         : P-D spectrum [W]
        main_mode    : dict of dominant mode (or None)
        """
        if p2p_params is None:
            p2p_params = P2P

        rp = self.radar_params
        Pt = float(rp.get('Pt_W', PT_W))
        Gt = float(rp.get('Gt',   GT))
        Gr = float(rp.get('Gr',   GR))

        # ── Step 1: build density field ───────────────────────────────────────
        Ne_2d, _  = self.iono.build_Ne_field(self.x_array, self.z_array, t)
        n_model   = RefractiveIndex(Ne_2d, self.x_array, self.z_array, self.freq)

        # ── Step 2: P2P variational ray paths ─────────────────────────────────
        rays = find_all_rays_p2p(tx_km, rx_km, n_model,
                                  self.freq, p2p_params=p2p_params)

        mode_results = []
        for ray in rays:
            gp           = float(ray['group_path_km'])
            tau_ms       = float(ray['tau_ms'])
            h_r          = float(ray['h_reflect_km'])
            label        = str(ray.get('label', classify_mode(ray)))
            pts          = ray['points']
            delta_tau_ms = 0.0

            # Free-space path loss
            L_free_dB = free_space_loss_dB(gp, self.freq)
            Pr_W      = Pt * Gt * Gr * 10.0 ** (-L_free_dB / 10.0)

            # ── Step 3a: Es model ─────────────────────────────────────────────
            if self.es is not None:
                es_cfg = self.iono_params.get('es_params') or {}
                h_Es   = float(es_cfg.get('h_Es_km', ES['h_Es_km']))
                es_par = extract_es_params(pts, h_Es)
                if es_par is not None:
                    theta = np.radians(max(float(es_par['theta_Es_deg']), 0.5))
                    es_res        = self.es.compute_power(
                        Pr_W, Gt, Gr, self.freq, gp / 2.0, theta)
                    Pr_W          = float(es_res['Pr_W'])
                    delta_tau_ms += float(es_res['delta_tau_ms'])
                    label         = 'Es_' + es_res['mode']

            # ── Step 3b: Plasma bubble (PE/SSF) ──────────────────────────────
            bub_cfg = self.iono_params.get('bubble_params') or {}
            if bub_cfg.get('enable', False):
                z0    = float(bub_cfg.get('z0_km', BUBBLE['z0_km']))
                Lz    = float(bub_cfg.get('Lz_km', BUBBLE['Lz_km']))
                Lx    = float(bub_cfg.get('Lx_km', BUBBLE['Lx_km']))
                h_bot = z0 - Lz

                bub_entry = extract_bubble_entry(pts, h_bot)
                if bub_entry is not None:
                    x_ent  = float(bub_entry['x'])
                    z_ent  = float(bub_entry['z'])
                    b_deg  = float(bub_entry['beta_inc_deg'])
                    A_inc  = _amplitude_from_power(Pr_W)

                    dx_km = float(PE['dx_km'])
                    dz_km = float(PE['dz_m']) / 1000.0

                    x_rng = (max(x_ent, float(self.x_array[0])),
                             min(x_ent + 2.0 * Lx, float(self.x_array[-1])))
                    z_rng = (max(h_bot, float(self.z_array[0])),
                             min(z0 + Lz, float(self.z_array[-1])))

                    n_pe, x_pe, z_pe = self.pe.extract_domain(
                        Ne_2d, self.x_array, self.z_array, x_rng, z_rng)

                    u_init    = construct_incident_field(
                        A_inc, b_deg, z_ent, z_pe, self.k0)
                    u_out, _  = self.pe.propagate(u_init, n_pe, dx_km, dz_km)

                    dx_tot    = (len(x_pe) - 1) * dx_km
                    pe_res    = self.pe.analyze(u_out, z_pe, dx_tot)
                    delta_tau_ms += float(pe_res['delta_tau_ms'])

                    P_in  = float(np.sum(np.abs(u_init) ** 2))
                    P_out = float(np.sum(np.abs(u_out)  ** 2))
                    if P_in > 0.0:
                        Pr_W *= P_out / P_in

            mode_results.append({
                'label'         : label,
                'tau_ms'        : tau_ms,
                'delta_tau_ms'  : delta_tau_ms,
                'Pr_W'          : Pr_W,
                'Pr_dBW'        : to_dBW(Pr_W),
                'h_reflect_km'  : h_r,
                'group_path_km' : gp,
                'beta_deg'      : float(ray.get('beta_deg', 0.0)),
                'phi_deg'       : 0.0,
                'points'        : pts,
            })

        # ── Steps 4-5: P-D spectrum + main mode ───────────────────────────────
        tau_axis, pd_W = build_pd_spectrum(mode_results)
        main_mode, _   = identify_main_mode(tau_axis, pd_W, mode_results)

        return mode_results, tau_axis, pd_W, main_mode


# ── Internal helper ───────────────────────────────────────────────────────────

def _split_iono(iono_params: dict) -> dict:
    """Extract IonosphereModel keyword arguments from iono_params."""
    keys = ('iri_params', 'tid_params', 'es_params', 'bubble_params')
    return {k: iono_params[k] for k in keys if k in iono_params}
