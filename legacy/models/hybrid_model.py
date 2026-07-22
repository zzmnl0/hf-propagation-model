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
from config import (FREQ_MHZ, PT_W, GT, GR, ES, BUBBLE,
                    BG_X, BG_Z, PE, P2P, RT, MODE, TX_POS, RX_POS, C_MS, EPS0,
                    RADAR, LINK_BEARING_DEG, GEOMAG, ABSORPTION, TUBE_TRACER,
                    UPDATE2)
from models.ionosphere_model import IonosphereModel
from .ray_tracer import RefractiveIndex
from .point_to_point import (find_all_rays_p2p, classify_mode,
                              extract_es_params, extract_bubble_entry)
from models.es_model import EsLayerModel
from .pe_propagator import PEPropagator, construct_incident_field
from .hybrid_ray_pe import apply_ray_pe_coupling
from .tube_tracer import TubeRayTracer
from .radar_clutter import apply_clutter_to_mode
from utils import (free_space_loss_dB, to_dBW, freq_to_k0,
                   radar_equation_W, d_layer_absorption_dB)
from models.pd_spectrum import build_pd_spectrum, identify_main_mode


# ── Field-amplitude helper ────────────────────────────────────────────────────

def _amplitude_from_power(Pr_W: float) -> float:
    """Approximate field amplitude [V/m] from received power [W]."""
    return float(np.sqrt(2.0 * max(Pr_W, 0.0) / (EPS0 * C_MS)))


# ── P-D spectrum ──────────────────────────────────────────────────────────────

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
                 iono_params:        dict,
                 radar_params:       dict,
                 x_array:            np.ndarray = BG_X,
                 z_array:            np.ndarray = BG_Z,
                 radar_mode:         bool = False,
                 geomag_params:      dict | None = None,
                 absorption_params:  dict | None = None,
                 tube_mode:          bool = False):
        self.iono_params  = iono_params
        self.radar_params = radar_params
        self.x_array      = x_array
        self.z_array      = z_array
        self.radar_mode   = radar_mode
        self.tube_mode    = tube_mode
        self.geomag       = geomag_params    # None -> isotropic; dict -> AH O/X
        self.absorption   = absorption_params if absorption_params is not None else {}
        self.freq         = float(radar_params.get('freq_MHz', FREQ_MHZ))
        self.k0           = freq_to_k0(self.freq)
        self.mode_summary = {}               # populated by compute()

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

        # ── Step 2: ray solver — tube tracer or P2P variational ──────────────
        if self.tube_mode and self.radar_mode:
            mode_results = self._compute_tube_modes(n_model)
        else:
            wm = ('both' if (self.geomag is not None and
                              self.geomag.get('enable_OX', False)) else None)
            rays = find_all_rays_p2p(tx_km, rx_km, n_model,
                                      self.freq, p2p_params=p2p_params,
                                      wave_mode=wm, geomag=self.geomag)

            mode_results = []
            for ray in rays:
                gp           = float(ray['group_path_km'])
                tau_ms       = float(ray['tau_ms'])
                h_r          = float(ray['h_reflect_km'])
                label        = str(ray.get('label', classify_mode(ray)))
                pts          = ray['points']
                delta_tau_ms = 0.0

                # Power: radar equation (monostatic) or Friis (communication)
                if self.radar_mode:
                    sigma = float(self.radar_params.get('sigma_rcs_m2',
                                                         RADAR['sigma_rcs_m2']))
                    Pr_W = radar_equation_W(Pt, Gt, Gr, self.freq, gp, sigma)
                else:
                    L_free_dB = free_space_loss_dB(gp, self.freq)
                    Pr_W = Pt * Gt * Gr * 10.0 ** (-L_free_dB / 10.0)

                # -- D-layer / complex collision absorption --------------------
                abs_oneway_dB = 0.0
                abs_total_dB = 0.0
                abs_model = 'none'
                abs_diag = {'abs_status': 'disabled'}
                complex_cfg = UPDATE2.get('ABSORPTION_COMPLEX', {})
                if complex_cfg.get('enable', False):
                    from .absorption_collision import complex_absorption_path_dB
                    from scipy.interpolate import RegularGridInterpolator
                    ne_interp = RegularGridInterpolator(
                        (self.x_array, self.z_array), Ne_2d,
                        bounds_error=False, fill_value=0.0)
                    def _sample_ne(x, y, z):
                        return float(ne_interp([[x, z]])[0])
                    abs_diag = complex_absorption_path_dB(
                        pts, self.freq, _sample_ne, complex_cfg,
                        two_way=bool(self.radar_mode),
                        wave_mode=ray.get('wave_mode', 'O'),
                        geomag=self.geomag or GEOMAG)
                    abs_oneway_dB = float(abs_diag.get('abs_complex_dB', 0.0))
                    abs_total_dB = float(abs_diag.get('abs_total_dB', abs_oneway_dB))
                    abs_model = 'appleton_hartree_collision_path'
                    Pr_W *= 10.0 ** (-abs_total_dB / 10.0)
                elif self.absorption.get('enable', False):
                    geomag  = self.geomag or GEOMAG
                    fH_L    = (float(geomag.get('fH_MHz', 1.197))
                               * float(np.cos(np.radians(geomag.get('dip_deg', 48.7)))))
                    ab_kw   = {k: v for k, v in self.absorption.items()
                               if k != 'enable'}
                    ab_kw['fH_L_MHz'] = fH_L
                    A_dB    = d_layer_absorption_dB(
                        self.freq, float(ray.get('beta_deg', 30.0)), **ab_kw)
                    n_cross = 2 if self.radar_mode else 1
                    abs_oneway_dB = float(A_dB)
                    abs_total_dB = float(A_dB * n_cross)
                    abs_model = 'legacy_d_layer'
                    abs_diag = {'abs_status': 'ok'}
                    Pr_W   *= 10.0 ** (-abs_total_dB / 10.0)
                # ── Step 3a: Es model ─────────────────────────────────────────
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

                # ── Step 3b: Plasma bubble (PE/SSF) ──────────────────────────
                pe_diag = {}
                bub_cfg = self.iono_params.get('bubble_params') or {}
                if bub_cfg.get('enable', False):
                    z0 = float(bub_cfg.get('z0_km', BUBBLE['z0_km']))
                    Lz = float(bub_cfg.get('Lz_km', BUBBLE['Lz_km']))
                    Lx = float(bub_cfg.get('Lx_km', BUBBLE['Lx_km']))
                    h_bot = z0 - Lz

                    bub_entry = extract_bubble_entry(pts, h_bot)
                    if bub_entry is not None:
                        x_ent = float(bub_entry['x'])
                        z_ent = float(bub_entry['z'])
                        b_deg = float(bub_entry['beta_inc_deg'])
                        pe_cfg = UPDATE2.get('PE_COUPLING', {})
                        if pe_cfg.get('enable', False):
                            Pr_W, pe_diag = apply_ray_pe_coupling(
                                self.pe, Ne_2d, self.x_array, self.z_array,
                                Pr_W, self.k0, x_ent, z_ent, b_deg, bub_cfg,
                                pe_cfg)
                            delta_tau_ms += float(pe_diag.get('pe_delta_tau_ms', 0.0))


                mode_results.append({
                    'label'         : label,
                    'tau_ms'        : tau_ms,
                    'tau_2way_ms'   : 2.0 * tau_ms if self.radar_mode else None,
                    'delta_tau_ms'  : delta_tau_ms,
                    'Pr_W'          : Pr_W,
                    'Pr_dBW'        : to_dBW(Pr_W),
                    'h_reflect_km'  : h_r,
                    'group_path_km' : gp,
                    'beta_deg'      : float(ray.get('beta_deg', 0.0)),
                    'phi_deg'       : float(self.radar_params.get('bearing_deg',
                                                                   LINK_BEARING_DEG)),
                    'wave_mode'     : ray.get('wave_mode', 'iso'),
                    'abs_dB'        : abs_oneway_dB,
                    'abs_total_dB'  : abs_total_dB,
                    'abs_model'     : abs_model,
                    'abs_complex_dB': float(abs_diag.get('abs_complex_dB', 0.0)),
                    'abs_status'    : abs_diag.get('abs_status', ''),
                    'collision_model': abs_diag.get('collision_model', ''),
                    'abs_index_model': abs_diag.get('abs_index_model', ''),
                    'points'        : pts,
                })

                if pe_diag:
                    mode_results[-1].update(pe_diag)
        # ── Steps 4-5: P-D spectrum + main mode ───────────────────────────────
        _pd_method = 'tube' if self.tube_mode else 'gaussian'
        tau_axis, pd_W = build_pd_spectrum(mode_results, method=_pd_method)
        main_mode, _, mode_summary = identify_main_mode(tau_axis, pd_W, mode_results)
        self.mode_summary = mode_summary

        return mode_results, tau_axis, pd_W, main_mode


    def _compute_tube_modes(self, n_model) -> list:
        """
        Flux tube radar pipeline (tube_mode=True, radar_mode=True).

        1. Extract sigma0 and x_tgt from RADAR config.
        2. Run TubeRayTracer.compute() to get physical backscatter modes.
        3. Apply D-layer absorption (2x traversal) and Es correction
           using the same logic as the P2P path.
        """
        rp      = self.radar_params
        Pt      = float(rp.get('Pt_W', PT_W))
        Gt      = float(rp.get('Gt',   GT))
        Gr      = float(rp.get('Gr',   GR))
        x_tgt   = float(RADAR['target_range_km'])
        s0_dB   = float(RADAR.get('sigma0_ground_dB', -20.0))
        sigma0  = 10.0 ** (s0_dB / 10.0)

        tracer     = TubeRayTracer(n_model, self.freq,
                                   tube_params=TUBE_TRACER, rt_params=RT)
        tube_modes = tracer.compute((0.0, 0.0), x_tgt, Pt, Gt, Gr, sigma0)

        for m in tube_modes:
            Pr_W = m['Pr_W']

            # D-layer absorption (2 traversals for radar round trip)
            if self.absorption.get('enable', False):
                geomag  = self.geomag or GEOMAG
                fH_L    = (float(geomag.get('fH_MHz', 1.197))
                           * float(np.cos(np.radians(geomag.get('dip_deg', 48.7)))))
                ab_kw   = {k: v for k, v in self.absorption.items()
                           if k != 'enable'}
                ab_kw['fH_L_MHz'] = fH_L
                A_dB    = d_layer_absorption_dB(
                    self.freq, float(m.get('beta_deg', 30.0)), **ab_kw)
                Pr_W   *= 10.0 ** (-(A_dB * 2) / 10.0)

            # Es correction (if ray path crosses Es layer)
            if self.es is not None:
                es_cfg = self.iono_params.get('es_params') or {}
                h_Es   = float(es_cfg.get('h_Es_km', ES['h_Es_km']))
                pts    = m.get('points')
                if pts is not None and len(pts) > 0:
                    es_par = extract_es_params(pts, h_Es)
                    if es_par is not None:
                        theta  = np.radians(max(float(es_par['theta_Es_deg']), 0.5))
                        es_res = self.es.compute_power(
                            Pr_W, Gt, Gr, self.freq,
                            m['group_path_km'] / 2.0, theta)
                        Pr_W              = float(es_res['Pr_W'])
                        m['delta_tau_ms'] = (float(m.get('delta_tau_ms', 0.0))
                                             + float(es_res['delta_tau_ms']))
                        m['label']        = 'Es_' + es_res['mode']

            m['Pr_W']   = Pr_W
            m['Pr_dBW'] = to_dBW(Pr_W)

            clutter_cfg = (UPDATE2.get('RADAR_CLUTTER', {})
                           if isinstance(UPDATE2, dict) else {})
            if clutter_cfg.get('enable', False):
                params = {k: v for k, v in clutter_cfg.items()
                          if k not in ('enable', 'apply_to_power')}
                updated = apply_clutter_to_mode(
                    m,
                    freq_MHz=self.freq,
                    params=params,
                    Pt_W=Pt,
                    Gt=Gt,
                    Gr=Gr,
                    apply_to_power=bool(clutter_cfg.get('apply_to_power', False)),
                )
                m.clear()
                m.update(updated)

        return tube_modes


# ── Internal helper ───────────────────────────────────────────────────────────

def _split_iono(iono_params: dict) -> dict:
    """Extract IonosphereModel keyword arguments from iono_params."""
    keys = ('iri_params', 'tid_params', 'es_params', 'bubble_params', 'spread_f_params')
    return {k: iono_params[k] for k in keys if k in iono_params}




