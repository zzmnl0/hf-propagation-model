"""
Module: 3-D point-to-point ray solver (two-parameter Newton shooting).

Finds rays that travel from TX (x=0, y=0, z=0) to RX (x=x_rx, y=0, z=0)
by iterating over launch (elevation beta, azimuth az) using Newton's method
with a 2x2 Jacobian computed via central differences.

Algorithm:
  1. Seed elevations from a coarse scan (same betas as 2-D P2P).
  2. For each seed, start with az_0 = 0 (in-plane ray).
  3. Newton iteration: update (beta, az) to minimise landing error
     F = (x_land - x_rx, y_land - y_rx).
  4. Deduplicate by (h_reflect, tau) tolerances.

New ray dict fields vs 2-D P2P:
  y_land_km, az_deg_tx, az_deg_rx, azimuth_deflect_deg,
  reflect_x_km, reflect_y_km, h_reflect_km (physical),
  elevation_tx_deg (= beta_deg_tx), elevation_rx_deg (= beta_deg_rx).

Reference: Cervera & Harris (2014) for 3-D P2P Newton shooting approach.
"""
import numpy as np
from config import C_KMS, P2P, RT, MODE, FREQ_MHZ
from .ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from .point_to_point import classify_mode as classify_mode_2d


# ── Newton shooter ─────────────────────────────────────────────────────────────

def _shoot_3d(tx: tuple,
              rx: tuple,
              n3d: RefractiveIndex3D,
              freq_MHz: float,
              beta_init: float,
              az_init:   float = 0.0,
              max_iter:  int   = 20,
              tol_km:    float = 5.0,
              d_beta:    float = 0.05,
              d_az:      float = 0.10,
              t:         float = 0.0
              ) -> dict | None:
    """
    Newton 2-parameter shooting: iterate (beta, az) to hit (x_rx, y_rx=0).

    Jacobian (2x2) via central differences:
      J11 = dx_land/dbeta,  J12 = dx_land/daz
      J21 = dy_land/dbeta,  J22 = dy_land/daz

    Step sizes: d_beta=0.05 deg, d_az=0.10 deg (Cervera & Harris 2014).
    Convergence: ||(x_err, y_err)|| < tol_km.
    """
    beta  = float(beta_init)
    az    = float(az_init)
    x_rx  = float(rx[0])
    y_rx  = 0.0

    for _ in range(max_iter):
        ray_c  = trace_single_ray_3d(tx, beta,        az,        n3d, freq_MHz, t=t)
        ray_bh = trace_single_ray_3d(tx, beta+d_beta, az,        n3d, freq_MHz, t=t)
        ray_bl = trace_single_ray_3d(tx, beta-d_beta, az,        n3d, freq_MHz, t=t)
        ray_ah = trace_single_ray_3d(tx, beta,        az+d_az,   n3d, freq_MHz, t=t)
        ray_al = trace_single_ray_3d(tx, beta,        az-d_az,   n3d, freq_MHz, t=t)

        x_c = ray_c['x_land_km'];  y_c = ray_c['y_land_km']
        dx_db = (ray_bh['x_land_km'] - ray_bl['x_land_km']) / (2.0 * d_beta)
        dy_db = (ray_bh['y_land_km'] - ray_bl['y_land_km']) / (2.0 * d_beta)
        dx_da = (ray_ah['x_land_km'] - ray_al['x_land_km']) / (2.0 * d_az)
        dy_da = (ray_ah['y_land_km'] - ray_al['y_land_km']) / (2.0 * d_az)

        J  = np.array([[dx_db, dx_da],
                        [dy_db, dy_da]])
        F  = np.array([x_rx - x_c, y_rx - y_c])
        err = np.linalg.norm(F)

        if err < tol_km:
            ray_c['beta_deg_tx']      = beta
            ray_c['az_deg_tx']        = az
            ray_c['elevation_tx_deg'] = beta
            ray_c['elevation_rx_deg'] = ray_c['beta_deg_rx']
            return ray_c

        det = J[0,0]*J[1,1] - J[0,1]*J[1,0]
        if abs(det) < 1e-12:
            break
        J_inv   = np.array([[ J[1,1], -J[0,1]],
                              [-J[1,0],  J[0,0]]]) / det
        delta   = J_inv @ F
        beta   += float(delta[0])
        az     += float(delta[1])
        beta    = float(np.clip(beta, 2.0, 88.0))

    return None


# ── Deduplication ─────────────────────────────────────────────────────────────

def _deduplicate_3d(candidates: list,
                    clust_h_km:   float,
                    clust_tau_ms: float) -> list:
    unique = []
    for c in candidates:
        is_dup = any(
            abs(c['h_reflect_km'] - u['h_reflect_km']) < clust_h_km and
            abs(c['tau_ms']       - u['tau_ms'])       < clust_tau_ms
            for u in unique
        )
        if not is_dup:
            unique.append(c)
    return unique


# ── Mode classification for 3-D rays ─────────────────────────────────────────

def classify_mode_3d(ray: dict) -> str:
    """
    Classify 3-D ray using same height/tau thresholds as 2-D,
    with optional O/X suffix from wave_mode field.
    """
    h   = ray.get('h_reflect_km', 0.0)
    tau = ray.get('tau_ms', 0.0)
    wm  = ray.get('wave_mode', 'iso')

    if h < MODE['h_Es_km']:
        base = 'Es'
    elif h < MODE['h_E_km']:
        base = 'E'
    elif h < 300.0:
        base = '1F_low' if tau < 5.0 else '1F_high'
    else:
        base = '2F'

    if wm in ('O', 'X'):
        return '{}_{}'.format(base, wm)
    return base


# ── Main P2P search ───────────────────────────────────────────────────────────

def find_all_rays_3d(tx:        tuple,
                     rx:        tuple,
                     n3d:       RefractiveIndex3D,
                     freq_MHz:  float = FREQ_MHZ,
                     p2p_params: dict = P2P,
                     t:         float = 0.0
                     ) -> list[dict]:
    """
    3-D point-to-point ray search using Newton two-parameter shooting.

    Seeds: n_init elevation angles in [5, 80] deg, azimuth seed = 0 deg.
    Each seed runs Newton iteration to converge to (x_rx, y_rx=0).
    Results are deduplicated by (h_reflect, tau) and sorted by tau.

    Parameters
    ----------
    tx          : (x0, z0) TX position [km, km]  (y=0 implicit)
    rx          : (x_rx, z_rx) RX position  (y=0 implicit)
    n3d         : RefractiveIndex3D instance
    freq_MHz    : float
    p2p_params  : dict  (uses n_init, clust_h_km, clust_tau_ms from P2P)
    t           : snapshot time [s]

    Returns list of ray dicts with 3-D fields including:
      tau_ms, h_reflect_km, group_path_km, label,
      elevation_tx_deg, elevation_rx_deg, azimuth_deflect_deg,
      reflect_x_km, reflect_y_km, y_land_km
    """
    n_init    = p2p_params.get('n_init', 18)
    clust_h   = p2p_params.get('clust_h_km', 10.0)
    clust_tau = p2p_params.get('clust_tau_ms', 0.05)
    tol_km    = p2p_params.get('tol_km', 5.0)

    betas = np.linspace(5.0, 80.0, n_init)
    candidates = []

    for beta_i in betas:
        for az_i in (0.0,):   # start in-plane; Newton corrects azimuth
            ray = _shoot_3d(tx, rx, n3d, freq_MHz,
                            beta_init=float(beta_i),
                            az_init=az_i,
                            tol_km=tol_km,
                            t=t)
            if ray is None:
                continue
            if ray['h_reflect_km'] > RT['z_stop_km']:
                continue
            if ray['tau_ms'] < 0.1:
                continue
            ray['wave_mode'] = 'iso'
            ray['label']     = classify_mode_3d(ray)
            candidates.append(ray)

    unique = _deduplicate_3d(candidates, clust_h, clust_tau)
    unique.sort(key=lambda r: r['tau_ms'])
    return unique
