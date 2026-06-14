"""
Module M5-M7 – Hybrid propagation model  (top-level integration).
Orchestrates the full pipeline:

  M1 (ionosphere) -> M2 (RT / P2P) -> M3 (Es) -> M4 (PE)
  -> M5 (synthesis) -> M6 (dual-path) -> M7 (main-mode ID)

Primary entry point:
    HybridPropagationModel.compute(tx, rx) -> (mode_results, tau, pd, main_mode)

Implemented in Part 6.
"""
import numpy as np
from scipy.signal import find_peaks
from config import (FREQ_MHZ, PT_W, GT, GR, ES, BUBBLE,
                    BG_X, BG_Z, PE, P2P, MODE)
from .ionosphere_model import IonosphereModel
from .ray_tracer import RefractiveIndex
from .point_to_point import find_all_rays_p2p, classify_mode, extract_es_params, extract_bubble_entry
from .es_model import EsLayerModel
from .pe_propagator import PEPropagator, construct_incident_field
from utils import free_space_loss_dB, to_dBW, group_delay_ms, freq_to_k0


# ── Loss helpers ──────────────────────────────────────────────────────────────

def _amplitude_from_power(Pr_W: float) -> float:
    """Approximate field amplitude from power [W] (isotropic far-field)."""
    from config import C_MS, EPS0
    return np.sqrt(2.0 * Pr_W / (EPS0 * C_MS * 1.0))   # unit area assumed


# ── P-D spectrum ──────────────────────────────────────────────────────────────

def build_pd_spectrum(mode_results: list,
                      tau_axis: np.ndarray | None = None,
                      tau_res_ms: float = 0.02
                      ) -> tuple[np.ndarray, np.ndarray]:
    """
    Construct the power-delay (P-D) spectrum by summing Gaussian-spread
    contributions from each propagation mode.

    Each mode contributes:
        Pr_i * G(tau - tau_i, dtau_i)   where G is a normalised Gaussian.

    Parameters
    ----------
    mode_results : list of mode dicts (must contain 'tau_ms', 'Pr_W', 'delta_tau_ms')
    tau_axis     : optional pre-defined delay axis [ms]; auto-generated if None
    tau_res_ms   : delay-axis resolution [ms]

    Returns (tau_axis, pd_W) both as 1-D arrays.
    """
    raise NotImplementedError("Implemented in Part 6")


def identify_main_mode(tau_axis: np.ndarray,
                        pd_W: np.ndarray,
                        mode_results: list
                        ) -> tuple[dict | None, list]:
    """
    Detect peaks in the P-D spectrum and map them back to mode dicts.

    Returns (main_mode_dict, all_modes_ranked_by_power).
    main_mode_dict is None if no peaks found.
    """
    raise NotImplementedError("Implemented in Part 6")


# ── Main class ────────────────────────────────────────────────────────────────

class HybridPropagationModel:
    """
    End-to-end hybrid HF propagation model.

    Parameters
    ----------
    iono_params   : dict – passed to IonosphereModel (TID/Es/bubble flags)
    radar_params  : dict – Pt_W, Gt, Gr, freq_MHz
    x_array, z_array : background grid arrays [km]
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
        self.freq         = radar_params.get('freq_MHz', FREQ_MHZ)
        self.k0           = freq_to_k0(self.freq)

        self.iono  = IonosphereModel(freq_MHz=self.freq, **_split_iono(iono_params))
        self.es    = (EsLayerModel(**iono_params['es'])
                      if iono_params.get('es', {}).get('enable') else None)
        self.pe    = PEPropagator(freq_MHz=self.freq)

    def compute(self,
                tx_km: tuple = (0.0, 0.0),
                rx_km: tuple = (1169.0, 0.0),
                t: float = 0.0
                ) -> tuple[list, np.ndarray, np.ndarray, dict | None]:
        """
        Run the full hybrid pipeline for a single TX-RX link.

        Pipeline
        --------
        1.  Build density field  Ne(x, z)
        2.  Solve P2P ray paths  (variational method)
        3.  For each ray:
              3a.  Apply Es model      if ray crosses Es height
              3b.  Apply PE/SSF        if ray enters plasma bubble
        4.  Synthesise P-D spectrum
        5.  Identify main mode

        Returns
        -------
        mode_results : list of mode dicts
        tau_axis     : delay axis [ms]
        pd_W         : P-D spectrum [W]
        main_mode    : dict of dominant mode  (or None)
        """
        raise NotImplementedError("Implemented in Part 6")


# ── Internal helper ───────────────────────────────────────────────────────────

def _split_iono(iono_params: dict) -> dict:
    """Extract IonosphereModel keyword arguments from iono_params."""
    keys = ('iri_params', 'tid_params', 'es_params', 'bubble_params')
    return {k: iono_params.get(k) for k in keys if k in iono_params}
