"""
Module M2 (P2P) – Point-to-point ray solver.
Two methods:
  A. Newton shooting   – fast, works well without strong irregularities.
  B. Variational (Nosikov 2020) – finds ALL modes (high + low ray)
     even when the medium is heavily disturbed (TID / plasma bubble).

Each returned ray dict has the same schema as ray_tracer.py, plus:
    'label'  : mode classification string (e.g. '1F_high', 'Es_reflect')
    'points' : (n_ctrl+2, 2) control-point array for variational rays

Implemented in Part 3.
"""
import numpy as np
from config import P2P, ES, BUBBLE, MODE
from .ray_tracer import RefractiveIndex, trace_single_ray


def find_ray_newton(tx: tuple,
                    rx: tuple,
                    n_model: RefractiveIndex,
                    freq_MHz: float,
                    beta_init_deg: float = 30.0,
                    max_iter: int = 50,
                    tol_km: float = 1.0
                    ) -> dict | None:
    """
    Newton shooting method.
    Iteratively adjust launch elevation until the ray lands on rx.

    Returns a ray dict on convergence, None otherwise.
    """
    raise NotImplementedError("Implemented in Part 3")


def optical_path(points: np.ndarray, n_model: RefractiveIndex) -> float:
    """
    Discrete optical-path functional S = sum n(midpoint) * |segment|.
    points : (N+1, 2) array of (x, z) control points  [km].
    """
    raise NotImplementedError("Implemented in Part 3")


def optical_path_gradient(points: np.ndarray,
                           n_model: RefractiveIndex,
                           fd_step: float = 0.2
                           ) -> np.ndarray:
    """
    Gradient dS/dr_i for interior control points i = 1..N-1.
    Returns array of shape (N-1, 2).
    (See Nosikov 2020, Eq. for Fermat's principle gradient.)
    """
    raise NotImplementedError("Implemented in Part 3")


def remove_tangential(grad: np.ndarray,
                      points: np.ndarray) -> np.ndarray:
    """
    Project gradient onto the plane perpendicular to the local path tangent.
    Keeps arc-length nearly constant during iteration.
    Returns grad_perp of shape (N-1, 2).
    """
    raise NotImplementedError("Implemented in Part 3")


def variational_find_ray(tx: tuple,
                          rx: tuple,
                          n_model: RefractiveIndex,
                          beta_init_deg: float,
                          is_high_ray: bool = True,
                          p2p_params: dict = P2P
                          ) -> tuple[np.ndarray, float]:
    """
    Variational (Nosikov 2020) solver starting from a single initial guess.

    high ray -> gradient descent  (S is minimum along path)
    low ray  -> spring-force flip (S is saddle point; flip sign to attract)

    Returns (points, group_path_km).
        points : (n_ctrl+2, 2) converged control-point array
        group_path_km : optical path length [km]
    """
    raise NotImplementedError("Implemented in Part 3")


def find_all_rays_p2p(tx: tuple,
                       rx: tuple,
                       n_model: RefractiveIndex,
                       freq_MHz: float,
                       p2p_params: dict = P2P
                       ) -> list[dict]:
    """
    Systematic search for ALL propagation modes between tx and rx.

    Steps:
        1. Scan n_init elevation angles -> initial straight-line paths.
        2. Run variational solver (both high and low variants).
        3. Cluster & deduplicate by (h_reflect, tau).
        4. Classify each unique mode.

    Returns list of ray dicts sorted by group delay.
    """
    raise NotImplementedError("Implemented in Part 3")


def classify_mode(ray_dict: dict) -> str:
    """
    Assign a text label to a propagation mode.

    Rules (from config.MODE):
        h_r < MODE['h_Es_km']               -> 'Es'
        MODE['h_Es_km'] <= h_r < h_E_km    -> 'E'
        h_r >= h_E_km, short delay          -> '1F_low'
        h_r >= h_E_km, long  delay          -> '1F_high'  (or '2F')
    """
    raise NotImplementedError("Implemented in Part 3")


def extract_es_params(points: np.ndarray,
                      h_Es_km: float) -> dict | None:
    """
    Find where the control-point path crosses h_Es_km (upward crossing).
    Returns dict with keys: x [km], z [km], theta_Es_deg [deg].
    Returns None if path never reaches h_Es_km.
    """
    raise NotImplementedError("Implemented in Part 3")


def extract_bubble_entry(points: np.ndarray,
                          z_bubble_bot_km: float) -> dict | None:
    """
    Find where the control-point path first crosses z_bubble_bot_km upward.
    Returns dict: x [km], z [km], beta_inc_deg [deg] (incidence elevation).
    Returns None if path never reaches that height.
    """
    raise NotImplementedError("Implemented in Part 3")
