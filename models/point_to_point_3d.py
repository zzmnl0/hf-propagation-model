"""
Module: 3-D point-to-point ray solver (two-parameter Newton shooting).

Finds rays that travel from TX (x=0, y=0, z=0) to RX (x=x_rx, y=0, z=0)
by iterating over launch (elevation beta, azimuth az) using Newton's method
with a 2x2 Jacobian computed via forward or central finite differences.

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
from config import C_KMS, P2P_3D, RT, MODE, FREQ_MHZ
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
              jacobian:  str   = 'forward',
              t:         float = 0.0
              ) -> dict | None:
    """
    Newton 2-parameter shooting: iterate (beta, az) to hit (x_rx, y_rx=0).

    Jacobian (2x2) via finite differences:
      J11 = dx_land/dbeta,  J12 = dx_land/daz
      J21 = dy_land/dbeta,  J22 = dy_land/daz

    Step sizes: d_beta=0.05 deg, d_az=0.10 deg (Cervera & Harris 2014).
    Convergence: ||(x_err, y_err)|| < tol_km.
    """
    beta  = float(beta_init)
    az    = float(az_init)
    x_rx  = float(rx[0])
    y_rx  = 0.0

    jacobian = str(jacobian).lower()
    if jacobian not in ('forward', 'central'):
        jacobian = 'forward'

    # Forward Jacobian: 3 ray traces per iteration; central uses 5.
    # Early termination: abandon seed if error is not decreasing after iter 3.
    prev_err = float('inf')
    for k in range(max_iter):
        ray_c  = trace_single_ray_3d(tx, beta,        az,        n3d, freq_MHz, t=t)

        x_c = ray_c['x_land_km'];  y_c = ray_c['y_land_km']
        if jacobian == 'central':
            ray_bp = trace_single_ray_3d(tx, beta+d_beta, az,        n3d, freq_MHz, t=t)
            ray_bm = trace_single_ray_3d(tx, beta-d_beta, az,        n3d, freq_MHz, t=t)
            ray_ap = trace_single_ray_3d(tx, beta,        az+d_az,   n3d, freq_MHz, t=t)
            ray_am = trace_single_ray_3d(tx, beta,        az-d_az,   n3d, freq_MHz, t=t)
            dx_db = (ray_bp['x_land_km'] - ray_bm['x_land_km']) / (2.0*d_beta)
            dy_db = (ray_bp['y_land_km'] - ray_bm['y_land_km']) / (2.0*d_beta)
            dx_da = (ray_ap['x_land_km'] - ray_am['x_land_km']) / (2.0*d_az)
            dy_da = (ray_ap['y_land_km'] - ray_am['y_land_km']) / (2.0*d_az)
        else:
            ray_bh = trace_single_ray_3d(tx, beta+d_beta, az,        n3d, freq_MHz, t=t)
            ray_ah = trace_single_ray_3d(tx, beta,        az+d_az,   n3d, freq_MHz, t=t)
            dx_db = (ray_bh['x_land_km'] - x_c) / d_beta
            dy_db = (ray_bh['y_land_km'] - y_c) / d_beta
            dx_da = (ray_ah['x_land_km'] - x_c) / d_az
            dy_da = (ray_ah['y_land_km'] - y_c) / d_az

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

        # Abandon diverging seeds early (saves ~60% of wasted iterations)
        if k >= 3 and err >= prev_err:
            break
        prev_err = err

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


def _seed_betas_3d(rx: tuple, p2p_params: dict, rays_fan: list | None = None) -> tuple[list[float], str]:
    """Return P2P seed elevations and the effective seed_source."""
    n_init = int(p2p_params.get('n_init', 18))
    seed_source = str(p2p_params.get('seed_source', 'internal_scan')).lower()
    if seed_source == 'fan_warm_start':
        if rays_fan:
            x_rx_f = float(rx[0])
            fan_sorted = sorted(rays_fan,
                                key=lambda r: abs(float(r.get('x_land_km', 0.0)) - x_rx_f))
            return [float(r.get('beta_deg_tx', r.get('beta_deg', 0.0)))
                    for r in fan_sorted[:min(n_init, len(fan_sorted))]], seed_source
        return list(np.linspace(5.0, 80.0, n_init)), 'internal_scan'
    if seed_source == 'explicit_list':
        betas = [float(v) for v in p2p_params.get('seed_betas_deg', [])]
        if betas:
            return betas, seed_source
        return list(np.linspace(5.0, 80.0, n_init)), 'internal_scan'
    return list(np.linspace(5.0, 80.0, n_init)), 'internal_scan'


def _cause_layer_from_height(label: str, h_km: float) -> str:
    if label.startswith('Es') or h_km < float(MODE['h_Es_km']):
        return 'Es'
    if h_km < float(MODE['h_E_km']):
        return 'E'
    if h_km < 300.0:
        return 'F1'
    return 'F2'


def annotate_mode_causes_3d(modes: list,
                            p2p_params: dict | None = None,
                            seed_source: str | None = None,
                            environment: str | None = None) -> list:
    """Attach P1 cause_* fields for plotting/CSV interpretation."""
    params = p2p_params or P2P_3D
    low_max = float(params.get('launch_low_max_deg', 25.0))
    high_min = float(params.get('launch_high_min_deg', 55.0))
    seed = seed_source or str(params.get('seed_source', 'internal_scan'))
    env = str(environment or '').lower()
    for m in modes:
        beta = float(m.get('beta_deg_tx', m.get('elevation_tx_deg', m.get('beta_deg', 0.0))))
        if beta <= low_max:
            launch = 'low_elevation'
        elif beta >= high_min:
            launch = 'high_elevation'
        else:
            launch = 'medium_elevation'
        label = str(m.get('label', ''))
        layer = _cause_layer_from_height(label, float(m.get('h_reflect_km', 0.0)))
        wave_mode = str(m.get('wave_mode', 'iso'))
        if wave_mode in ('O', 'X'):
            mechanism = 'magnetoionic_split'
        elif bool(m.get('pe_coupled', False)):
            mechanism = 'ray_pe_coupled'
        elif layer == 'Es':
            mechanism = 'layered_reflection'
        elif env in ('bubble', 'spreadf', 'full') and m.get('scatter_method'):
            mechanism = 'irregularity_scatter'
        else:
            mechanism = 'elevation_branch'
        m['cause_layer'] = layer
        m['cause_launch'] = launch
        m['cause_mechanism'] = mechanism
        m['cause_summary'] = '{} {} branch'.format(layer, launch.replace('_', '-'))
        if mechanism in ('magnetoionic_split', 'irregularity_scatter', 'ray_pe_coupled'):
            m['cause_summary'] += ' ({})'.format(mechanism)
        m['p2p_seed_source'] = seed
    return modes


# ── Main P2P search ───────────────────────────────────────────────────────────

def find_all_rays_3d(tx:        tuple,
                     rx:        tuple,
                     n3d:       RefractiveIndex3D,
                     freq_MHz:  float = FREQ_MHZ,
                     p2p_params: dict = P2P_3D,
                     t:         float = 0.0,
                     rays_fan:  list  = None
                     ) -> list[dict]:
    """
    3-D point-to-point ray search using Newton two-parameter shooting.

    Seeds: controlled by p2p_params['seed_source']. Default internal_scan uses an independent elevation scan; fan_warm_start is explicit only.
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
    rays_fan    : optional list of fan rays for warm-start seeding; if provided
                  the n_init rays nearest to x_rx are used as Newton seeds

    Returns list of ray dicts with 3-D fields including:
      tau_ms, h_reflect_km, group_path_km, label,
      elevation_tx_deg, elevation_rx_deg, azimuth_deflect_deg,
      reflect_x_km, reflect_y_km, y_land_km
    """
    n_init    = p2p_params.get('n_init', 18)
    clust_h   = p2p_params.get('clust_h_km', 10.0)
    clust_tau = p2p_params.get('clust_tau_ms', 0.05)
    tol_km    = p2p_params.get('tol_km', 5.0)
    max_iter  = p2p_params.get('max_iter', 12)
    d_beta    = p2p_params.get('d_beta', 0.05)
    d_az      = p2p_params.get('d_az',   0.10)
    jacobian  = p2p_params.get('jacobian', 'forward')

    betas, seed_source = _seed_betas_3d(rx, p2p_params, rays_fan)
    candidates = []

    for beta_i in betas:
        for az_i in (0.0,):   # start in-plane; Newton corrects azimuth
            ray = _shoot_3d(tx, rx, n3d, freq_MHz,
                            beta_init=float(beta_i),
                            az_init=az_i,
                            max_iter=max_iter,
                            tol_km=tol_km,
                            d_beta=d_beta,
                            d_az=d_az,
                            jacobian=jacobian,
                            t=t)
            if ray is None:
                continue
            if ray['h_reflect_km'] > RT['z_stop_km']:
                continue
            if ray['tau_ms'] < 0.1:
                continue
            ray['wave_mode'] = 'iso'
            ray['label']     = classify_mode_3d(ray)
            ray['p2p_seed_source'] = seed_source
            candidates.append(ray)

    unique = _deduplicate_3d(candidates, clust_h, clust_tau)
    unique.sort(key=lambda r: r['tau_ms'])
    annotate_mode_causes_3d(unique, p2p_params, seed_source)
    return unique
