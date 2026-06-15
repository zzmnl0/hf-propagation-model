"""
Module M2 (P2P) – Point-to-point ray solver.
Two methods:
  A. Newton shooting   – fast, works well without strong irregularities.
  B. Variational (Nosikov 2020) – finds ALL modes (high + low ray)
     even when the medium is heavily disturbed (TID / plasma bubble).

Each returned ray dict has the same schema as ray_tracer.py, plus:
    'label'  : mode classification string (e.g. '1F_high', 'Es_reflect')
    'points' : (n_ctrl+2, 2) control-point array for variational rays
"""
import numpy as np
from config import C_KMS, P2P, RT, MODE
from .ray_tracer import RefractiveIndex, trace_single_ray


# ── Newton shooting ───────────────────────────────────────────────────────────

def find_ray_newton(tx: tuple,
                    rx: tuple,
                    n_model: RefractiveIndex,
                    freq_MHz: float,
                    beta_init_deg: float = 30.0,
                    max_iter: int = 50,
                    tol_km: float = 1.0
                    ) -> dict | None:
    """
    Newton shooting: iteratively adjust launch elevation until the ray
    lands within tol_km of rx.  Returns ray dict on convergence, else None.
    """
    beta   = float(beta_init_deg)
    d_beta = 1.0

    for _ in range(max_iter):
        ray  = trace_single_ray(tx, beta,          n_model, freq_MHz=freq_MHz)
        ray2 = trace_single_ray(tx, beta + d_beta, n_model, freq_MHz=freq_MHz)

        traj  = np.array(ray['trajectory'])
        traj2 = np.array(ray2['trajectory'])

        x_land  = traj[-1, 0]
        x_land2 = traj2[-1, 0]

        dx_dbeta = (x_land2 - x_land) / d_beta
        delta_x  = float(rx[0]) - x_land

        if abs(dx_dbeta) < 1e-6:
            break
        beta += delta_x / dx_dbeta

        if abs(delta_x) < tol_km:
            ray['label']  = classify_mode(ray)
            ray['points'] = None
            return ray
    return None


# ── Optical-path functional & gradient ───────────────────────────────────────

def optical_path(points: np.ndarray, n_model: RefractiveIndex) -> float:
    """
    Discrete Fermat functional  S = sum_i n(mid_i) * |seg_i|.
    points : (N, 2) array of (x, z) control points [km].
    """
    pts  = np.asarray(points, dtype=float)
    mids = 0.5 * (pts[:-1] + pts[1:])          # (N-1, 2)
    segs = np.linalg.norm(pts[1:] - pts[:-1], axis=1)
    return float(np.dot(n_model.n_batch(mids), segs))


def optical_path_gradient(points: np.ndarray,
                           n_model: RefractiveIndex,
                           fd_step: float = 0.2
                           ) -> np.ndarray:
    """
    Gradient dS/dr_i for interior points i = 1..N-2.
    Returns (N-2, 2). Uses a single vectorized n_batch call per invocation.

    dS/dr_i = n_L * dir_L - n_R * dir_R + 0.5*(|seg_L|+|seg_R|) * dn/dr_i
    """
    pts = np.asarray(points, dtype=float)
    N   = len(pts)
    M   = N - 2                              # number of interior points

    mids    = 0.5 * (pts[:-1] + pts[1:])    # (N-1, 2) segment midpoints
    pts_int = pts[1:-1]                      # (M, 2)   interior points

    # FD perturbations at interior points
    q_px = pts_int.copy(); q_px[:, 0] += fd_step
    q_mx = pts_int.copy(); q_mx[:, 0] -= fd_step
    q_pz = pts_int.copy(); q_pz[:, 1] += fd_step
    q_mz = pts_int.copy(); q_mz[:, 1] -= fd_step

    all_q = np.vstack([mids, q_px, q_mx, q_pz, q_mz])
    all_n  = n_model.n_batch(all_q)

    n_mid = all_n[:N-1]
    n_px  = all_n[N-1        : N-1+M]
    n_mx  = all_n[N-1+M      : N-1+2*M]
    n_pz  = all_n[N-1+2*M    : N-1+3*M]
    n_mz  = all_n[N-1+3*M    :]

    segs = pts[1:] - pts[:-1]
    lens = np.linalg.norm(segs, axis=1) + 1e-12
    dirs = segs / lens[:, None]

    n_L   = n_mid[:-1]    # (M,) left-segment midpoint n
    n_R   = n_mid[1:]     # (M,) right-segment midpoint n
    dir_L = dirs[:-1]     # (M, 2)
    dir_R = dirs[1:]      # (M, 2)
    len_L = lens[:-1]     # (M,)
    len_R = lens[1:]      # (M,)

    # Direction term
    grad = n_L[:, None] * dir_L - n_R[:, None] * dir_R

    # Density-gradient term (FD at interior point)
    dn_dx  = (n_px - n_mx) / (2.0 * fd_step)
    dn_dz  = (n_pz - n_mz) / (2.0 * fd_step)
    avg_l  = 0.5 * (len_L + len_R)
    grad  += avg_l[:, None] * np.column_stack([dn_dx, dn_dz])

    return grad   # (M, 2)


def remove_tangential(grad: np.ndarray,
                      points: np.ndarray) -> np.ndarray:
    """
    Remove the tangential component of grad so that updates don't
    redistribute control points along the same curve (Nosikov 2020).
    grad: (M, 2), points: (M+2, 2) including fixed endpoints.
    Returns grad_perp of shape (M, 2).
    """
    pts = np.asarray(points, dtype=float)
    tangents = pts[2:] - pts[:-2]          # central-diff tangent at each interior point
    t_norms  = tangents / (np.linalg.norm(tangents, axis=1, keepdims=True) + 1e-12)
    proj = np.sum(grad * t_norms, axis=1, keepdims=True)
    return grad - proj * t_norms


# ── Variational solver (Nosikov 2020) ─────────────────────────────────────────

def variational_find_ray(tx: tuple,
                          rx: tuple,
                          n_model: RefractiveIndex,
                          beta_init_deg: float,
                          is_high_ray: bool = True,
                          p2p_params: dict = P2P
                          ) -> tuple[np.ndarray, float]:
    """
    Variational solver from one initial elevation guess.
      is_high_ray=True  -> gradient descent  (S is a local minimum)
      is_high_ray=False -> sign-flipped grad (S saddle -> attractor)

    Returns (points, group_path_km):
        points : (n_ctrl+2, 2) converged control-point array  [km]
    """
    n_ctrl = p2p_params['n_ctrl']
    alpha  = p2p_params['alpha_km']
    k_spr  = p2p_params['k_spring']
    max_it = p2p_params['max_iter']
    tol    = p2p_params['tol_km']

    beta_rad = np.radians(float(beta_init_deg))
    t = np.linspace(0.0, 1.0, n_ctrl + 2)

    x_range = float(rx[0]) - float(tx[0])
    h_peak  = float(np.tan(beta_rad)) * x_range / 2.0
    h_peak  = float(np.clip(h_peak, 150.0, 500.0))

    # Parabolic arc: z = h_peak * 4t(1-t)
    pts = np.column_stack([
        float(tx[0]) + t * x_range,
        float(tx[1]) + h_peak * 4.0 * t * (1.0 - t),
    ])

    for _ in range(max_it):
        grad      = optical_path_gradient(pts, n_model)
        grad_perp = remove_tangential(grad, pts)
        spring    = pts[2:] - 2.0 * pts[1:-1] + pts[:-2]

        if is_high_ray:
            pts[1:-1] += -alpha * grad_perp + k_spr * spring
        else:
            pts[1:-1] +=  alpha * grad_perp + k_spr * spring

        pts[1:-1, 1] = np.maximum(pts[1:-1, 1], 0.0)

        if np.max(np.linalg.norm(grad_perp, axis=1)) < tol:
            break

    gp = optical_path(pts, n_model)
    return pts, gp


# ── P2P mode search ───────────────────────────────────────────────────────────

def _deduplicate(candidates: list,
                 clust_h_km: float,
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


def _ray_worker(args: tuple) -> dict | None:
    """
    Module-level worker for multiprocessing.Pool.
    Must stay at module level (not nested) so it is picklable on Windows.
    """
    tx, rx, n_model, beta, is_high, p2p_params = args
    pts, gp = variational_find_ray(tx, rx, n_model,
                                    float(beta),
                                    is_high_ray=is_high,
                                    p2p_params=p2p_params)
    tau_ms = gp / C_KMS * 1e3
    h_max  = float(np.max(pts[:, 1]))

    if h_max > RT['z_stop_km'] or tau_ms < 0.1:
        return None

    return {
        'beta_deg'      : float(beta),
        'points'        : pts,
        'trajectory'    : [np.array([p[0], p[1], 0.0, 0.0]) for p in pts],
        'group_path_km' : gp,
        'tau_ms'        : tau_ms,
        'h_reflect_km'  : h_max,
        'beta_recv_deg' : 0.0,
        'at_Es'         : None,
        'at_bubble'     : None,
        'L_bg_dB'       : 0.0,
        'is_high'       : is_high,
        'label'         : '',
    }


def find_all_rays_p2p(tx: tuple,
                       rx: tuple,
                       n_model: RefractiveIndex,
                       freq_MHz: float,
                       p2p_params: dict = P2P
                       ) -> list[dict]:
    """
    Systematic search: scan n_init elevations × {high, low} ray variants,
    deduplicate, classify, and sort by tau.

    Parallelism: controlled by p2p_params['n_workers'].
      1          -> sequential (safe in any calling context)
      0          -> auto-detect CPU count (requires if __name__=='__main__' guard)
      N > 1      -> use exactly N worker processes
    """
    n_init    = p2p_params['n_init']
    clust_h   = p2p_params['clust_h_km']
    clust_tau = p2p_params['clust_tau_ms']
    n_workers = int(p2p_params.get('n_workers', 1))

    betas = np.linspace(5.0, 80.0, n_init)
    tasks = [(tx, rx, n_model, float(b), is_high, p2p_params)
             for b in betas
             for is_high in (True, False)]

    if n_workers == 1:
        raw = [_ray_worker(t) for t in tasks]
    else:
        from multiprocessing import Pool, cpu_count
        n_proc = n_workers if n_workers > 0 else min(cpu_count(), len(tasks))
        with Pool(processes=n_proc) as pool:
            raw = pool.map(_ray_worker, tasks)

    candidates = [r for r in raw if r is not None]
    unique = _deduplicate(candidates, clust_h, clust_tau)
    for r in unique:
        r['label'] = classify_mode(r)
    unique.sort(key=lambda r: r['tau_ms'])
    return unique


# ── Mode classification ───────────────────────────────────────────────────────

def classify_mode(ray_dict: dict) -> str:
    """
    Assign a propagation mode label from reflection height and group delay.
    Thresholds follow config.MODE plus the 300 km boundary for 2F.
    """
    h    = ray_dict['h_reflect_km']
    tau  = ray_dict['tau_ms']
    h_Es = MODE['h_Es_km']
    h_E  = MODE['h_E_km']

    if h < h_Es:
        return 'Es'
    elif h < h_E:
        return 'E'
    elif h < 300.0:
        return '1F_low' if tau < 5.0 else '1F_high'
    else:
        return '2F'


# ── Es / bubble crossing extraction ──────────────────────────────────────────

def extract_es_params(points: np.ndarray,
                      h_Es_km: float) -> dict | None:
    """
    Find first upward crossing of h_Es_km in the control-point path.
    Returns {'x', 'z', 'theta_Es_deg'} or None.
    """
    pts = np.asarray(points, dtype=float)
    for i in range(1, len(pts)):
        if pts[i-1, 1] < h_Es_km <= pts[i, 1]:
            t    = (h_Es_km - pts[i-1, 1]) / (pts[i, 1] - pts[i-1, 1])
            x_es = pts[i-1, 0] + t * (pts[i, 0] - pts[i-1, 0])
            dx   = pts[i, 0] - pts[i-1, 0]
            dz   = pts[i, 1] - pts[i-1, 1]
            theta = float(np.degrees(np.arctan2(abs(dz), abs(dx) + 1e-12)))
            return {'x': float(x_es), 'z': float(h_Es_km),
                    'theta_Es_deg': theta}
    return None


def extract_bubble_entry(points: np.ndarray,
                          z_bubble_bot_km: float) -> dict | None:
    """
    Find first upward crossing of z_bubble_bot_km.
    Returns {'x', 'z', 'beta_inc_deg'} or None.
    """
    pts = np.asarray(points, dtype=float)
    for i in range(1, len(pts)):
        if pts[i-1, 1] < z_bubble_bot_km <= pts[i, 1]:
            t     = (z_bubble_bot_km - pts[i-1, 1]) / (pts[i, 1] - pts[i-1, 1])
            x_bub = pts[i-1, 0] + t * (pts[i, 0] - pts[i-1, 0])
            dx    = pts[i, 0] - pts[i-1, 0]
            dz    = pts[i, 1] - pts[i-1, 1]
            beta_inc = float(np.degrees(np.arctan2(abs(dz), abs(dx) + 1e-12)))
            return {'x': float(x_bub), 'z': float(z_bubble_bot_km),
                    'beta_inc_deg': beta_inc}
    return None
