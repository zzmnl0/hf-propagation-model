"""
models/tube_tracer.py
Flux tube ray tracer for OTH radar backscatter (Coleman 1997/1998).

Shoots a fan of rays from TX, identifies adjacent ray pairs (tubes) whose
ground landing points bracket the target range x_tgt, and computes physical
backscatter power including focusing/defocusing effects.

Key formulas (Coleman 1997, Radio Science 32(1)):
  F_focus    = R_eff * d_beta / d_x_land    [dimensionless, >1 = focusing]
  A_tube     = d_x_land * L_cross_km        [km^2, ground projection area]
  Pr         = Pt*Gt*Gr*lam^2*sigma0*A_tube / ((4*pi)^3 * P_one^4)
  f_pulse    = T_pulse / tau_spread  if tau_spread > T_pulse else 1.0

Each mode dict returned by compute() is compatible with
HybridPropagationModel mode_results; new fields vs. standard dict:
  tau_2way_ms, tau_spread_ms, F_focus, A_tube_km2
"""
import numpy as np
from .ray_tracer import trace_single_ray
from .point_to_point import classify_mode
from config import RT, FREQ_MHZ, C_MS, TUBE_TRACER
from utils import to_dBW


class TubeRayTracer:
    """
    Flux tube ray tracer for monostatic OTH radar backscatter.

    Parameters
    ----------
    n_model     : RefractiveIndex    electron density interpolator
    freq_MHz    : float              operating frequency [MHz]
    tube_params : dict               TUBE_TRACER config (None -> use config default)
    rt_params   : dict               RT config (None -> use config default)
    """

    def __init__(self, n_model, freq_MHz=FREQ_MHZ,
                 tube_params=None, rt_params=None):
        self.n_model     = n_model
        self.freq_MHz    = float(freq_MHz)
        self.tube_params = tube_params if tube_params is not None else TUBE_TRACER
        self.rt_params   = rt_params   if rt_params   is not None else RT
        self.lambda_m    = C_MS / (self.freq_MHz * 1e6)   # wavelength [m]

    # ── Private helpers ───────────────────────────────────────────────────────

    def _land_x(self, ray):
        """Return final x [km] of ray (ground landing point)."""
        traj = ray.get('trajectory', [])
        if not traj:
            return 0.0
        return float(traj[-1][0])

    def _trace(self, tx_pos, beta_deg):
        """Trace a single ray; thin wrapper around trace_single_ray."""
        return trace_single_ray(
            tx_pos, float(beta_deg), self.n_model,
            freq_MHz=self.freq_MHz,
            rt_params=self.rt_params,
        )

    # ── Step 2: Fan ray shooting ──────────────────────────────────────────────

    def shoot_fan(self, tx_pos):
        """
        Shoot fan of rays at delta_beta_deg spacing from beta_min to beta_max.

        Each returned ray dict is augmented with 'x_land' (final x [km]).
        Rays are returned in ascending beta order.
        """
        tp    = self.tube_params
        db    = float(tp.get('delta_beta_deg', 0.5))
        n_max = int(tp.get('n_tube_rays', 80))
        rt    = self.rt_params
        b_min = float(rt.get('beta_min', 5.0))
        b_max = float(rt.get('beta_max', 85.0))

        betas = np.arange(b_min, b_max + db * 0.5, db)
        betas = betas[:n_max]

        rays = []
        for b in betas:
            r = self._trace(tx_pos, float(b))
            r['x_land'] = self._land_x(r)
            rays.append(r)
        return rays

    # ── Step 3: Tube geometry ─────────────────────────────────────────────────

    def compute_tubes(self, fan_rays, x_tgt):
        """
        Form flux tubes from adjacent fan rays near target x_tgt.

        A tube is formed by pair (r_i, r_{i+1}) when at least one landing
        point satisfies |x_land - x_tgt| <= tol_km.

        Returns list of tube dicts:
          beta_deg, beta_deg_hi, x_land, tau_ms, group_path_km,
          h_reflect_km, delta_x_km, A_tube_km2, F_focus,
          tau_spread_ms, label, ray
        """
        tp      = self.tube_params
        tol_km  = float(tp.get('x_tgt_tol_km', 80.0))
        L_cross = float(tp.get('L_cross_km', 100.0))
        db_rad  = float(np.radians(tp.get('delta_beta_deg', 0.5)))

        tubes = []
        for i in range(len(fan_rays) - 1):
            r0 = fan_rays[i]
            r1 = fan_rays[i + 1]
            x0 = float(r0['x_land'])
            x1 = float(r1['x_land'])

            # At least one ray must land near x_tgt
            if abs(x0 - x_tgt) > tol_km and abs(x1 - x_tgt) > tol_km:
                continue

            delta_x = abs(x1 - x0)
            if delta_x < 1e-6:
                continue

            # F_focus = free-space divergence / actual divergence
            # R_eff = one-way group path so that R_eff*d_beta = free-space
            # lateral separation at target range (Coleman 1997, Eq. 15-17)
            R_eff_km      = float(r0['group_path_km'])
            A_tube_km2    = delta_x * L_cross
            F_focus       = (R_eff_km * db_rad) / delta_x
            tau_spread_ms = abs(float(r1['tau_ms']) - float(r0['tau_ms']))
            label         = classify_mode(r0)

            tubes.append({
                'beta_deg'      : float(r0['beta_deg']),
                'beta_deg_hi'   : float(r1['beta_deg']),
                'x_land'        : 0.5 * (x0 + x1),
                'tau_ms'        : float(r0['tau_ms']),
                'group_path_km' : float(r0['group_path_km']),
                'h_reflect_km'  : float(r0['h_reflect_km']),
                'delta_x_km'    : delta_x,
                'A_tube_km2'    : A_tube_km2,
                'F_focus'       : F_focus,
                'tau_spread_ms' : tau_spread_ms,
                'label'         : label,
                'ray'           : r0,
            })

        return tubes

    # ── Step 4: Backscatter power and pulse correction ────────────────────────

    def backscatter_power_W(self, Pt_W, Gt, Gr, P_one_km, sigma0, A_tube_km2):
        """
        Coleman (1997) backscatter power.

        Pr = Pt * Gt * Gr * lam^2 * sigma0 * A_tube
             / ((4*pi)^3 * P_one^4)

        Uses (4*pi)^3 to match radar_equation_W convention:
        sigma0 * A_tube serves as the effective target RCS [m^2].

        Parameters
        ----------
        P_one_km   : one-way group path [km]
        sigma0     : ground normalized RCS [linear]
        A_tube_km2 : ground projection area [km^2]
        """
        P_m  = float(P_one_km) * 1e3         # [m]
        A_m2 = float(A_tube_km2) * 1e6       # [m^2]

        if P_m <= 0.0 or A_m2 <= 0.0:
            return 0.0

        num   = Pt_W * Gt * Gr * (self.lambda_m ** 2) * sigma0 * A_m2
        denom = ((4.0 * np.pi) ** 3) * (P_m ** 4)
        return float(num / max(denom, 1e-300))

    def pulse_correction(self, tau_spread_ms):
        """
        Pulse broadening correction factor f_pulse in (0, 1].

        When intra-tube delay spread exceeds the radar pulse width,
        the effective echo amplitude is reduced by T_pulse / tau_spread.
        """
        T_pulse = float(self.tube_params.get('T_pulse_ms', 0.5))
        if tau_spread_ms > T_pulse and tau_spread_ms > 0.0:
            return T_pulse / tau_spread_ms
        return 1.0

    # ── Step 5: Newton landing-point refinement ───────────────────────────────

    def newton_refine(self, tx_pos, beta0, x_tgt,
                      max_iter=None, tol_km=None):
        """
        1-D Newton iteration: find beta such that x_land(beta) = x_tgt.

          F(beta) = x_land(beta) - x_tgt = 0
          beta_{n+1} = beta_n - F / (dF/dbeta)
          dF/dbeta estimated by central finite difference (d_beta = 0.01 deg)

        Returns (refined_ray, beta_final).
        Clips beta to [beta_min, beta_max]; returns best solution on timeout.
        """
        tp     = self.tube_params
        max_it = max_iter if max_iter is not None else int(tp.get('newton_max_iter', 5))
        tol    = tol_km   if tol_km   is not None else float(tp.get('newton_tol_km', 1.0))
        d_beta = 0.01   # central-difference step [deg]

        rt    = self.rt_params
        b_min = float(rt.get('beta_min', 5.0))
        b_max = float(rt.get('beta_max', 85.0))

        beta = float(np.clip(beta0, b_min, b_max))
        ray  = self._trace(tx_pos, beta)
        xl   = self._land_x(ray)

        for _ in range(max_it):
            err = xl - x_tgt
            if abs(err) <= tol:
                break
            r_p  = self._trace(tx_pos, beta + d_beta)
            r_m  = self._trace(tx_pos, beta - d_beta)
            dxdb = (self._land_x(r_p) - self._land_x(r_m)) / (2.0 * d_beta)
            if abs(dxdb) < 1e-10:
                break
            beta = float(np.clip(beta - err / dxdb, b_min, b_max))
            ray  = self._trace(tx_pos, beta)
            xl   = self._land_x(ray)

        ray['x_land'] = xl
        return ray, beta

    # ── Step 6: Main compute pipeline ────────────────────────────────────────

    def compute(self, tx_pos, x_tgt, Pt_W, Gt, Gr, sigma0, newton=True):
        """
        Full flux tube pipeline.

        Steps
        -----
        A. Shoot fan of rays from tx_pos
        B. Identify adjacent ray pairs (tubes) landing near x_tgt
        C. Optionally Newton-refine the center ray of each tube
        D. Compute Coleman backscatter power with pulse correction
        E. Deduplicate: same label + |dtau| < 0.05 ms -> keep max power

        Parameters
        ----------
        tx_pos  : (x0, z0) [km]
        x_tgt   : target horizontal distance [km]
        Pt_W    : transmit power [W]
        Gt, Gr  : antenna gains [linear]
        sigma0  : ground normalized RCS [linear, not dB]
        newton  : bool  enable Newton center-ray refinement

        Returns
        -------
        list of mode dicts (compatible with HybridPropagationModel)
        """
        # A: Fan shooting
        fan_rays = self.shoot_fan(tx_pos)

        # B: Form tubes
        tubes = self.compute_tubes(fan_rays, x_tgt)
        if not tubes:
            return []

        # C + D: Per-tube power computation
        raw = []
        for t in tubes:
            tau_sp = float(t['tau_spread_ms'])
            A_tube = float(t['A_tube_km2'])
            F_foc  = float(t['F_focus'])

            if newton:
                ref_ray, beta_ref = self.newton_refine(
                    tx_pos, t['beta_deg'], x_tgt)
                gp    = float(ref_ray['group_path_km'])
                tau   = float(ref_ray['tau_ms'])
                h_r   = float(ref_ray['h_reflect_km'])
                beta  = beta_ref
                pts   = np.array([[s[0], s[1]]
                                  for s in ref_ray['trajectory']], dtype=float)
                label = classify_mode(ref_ray)
            else:
                ray   = t['ray']
                gp    = float(ray['group_path_km'])
                tau   = float(ray['tau_ms'])
                h_r   = float(ray['h_reflect_km'])
                beta  = float(ray['beta_deg'])
                pts   = np.array([[s[0], s[1]]
                                  for s in ray['trajectory']], dtype=float)
                label = t['label']

            f_p  = self.pulse_correction(tau_sp)
            Pr_W = self.backscatter_power_W(Pt_W, Gt, Gr, gp, sigma0, A_tube)
            Pr_W *= f_p

            raw.append({
                'label'         : label,
                'tau_ms'        : tau,
                'tau_2way_ms'   : 2.0 * tau,
                'tau_spread_ms' : tau_sp,
                'delta_tau_ms'  : tau_sp,      # alias for build_pd_spectrum
                'Pr_W'          : Pr_W,
                'Pr_dBW'        : to_dBW(Pr_W),
                'h_reflect_km'  : h_r,
                'group_path_km' : gp,
                'beta_deg'      : beta,
                'phi_deg'       : 0.0,
                'wave_mode'     : 'iso',
                'F_focus'       : F_foc,
                'A_tube_km2'    : A_tube,
                'points'        : pts,
            })

        # E: Deduplicate
        return _dedup_modes(raw, tau_tol_ms=0.05)


# ── Module-level deduplication helper ────────────────────────────────────────

def _dedup_modes(results, tau_tol_ms=0.05):
    """
    Merge entries sharing the same label and nearly-equal tau_ms
    (within tau_tol_ms).  Keeps the highest-power entry per cluster.
    """
    if not results:
        return []

    sorted_r = sorted(results, key=lambda m: (m['label'], m['tau_ms']))

    merged = []
    i = 0
    while i < len(sorted_r):
        anchor = sorted_r[i]
        best   = anchor
        j = i + 1
        while (j < len(sorted_r)
               and sorted_r[j]['label'] == anchor['label']
               and abs(sorted_r[j]['tau_ms'] - anchor['tau_ms']) < tau_tol_ms):
            if sorted_r[j]['Pr_W'] > best['Pr_W']:
                best = sorted_r[j]
            j += 1
        merged.append(best)
        i = j

    return merged
