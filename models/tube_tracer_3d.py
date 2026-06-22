"""
models/tube_tracer_3d.py
3-D flux-tube ray tracer (Coleman 1997/1998 extended to 3-D elliptic cross-section).

Extends 2-D TubeRayTracer to 3-D by computing focusing in both principal planes:

  F_focus_el = R_eff * d_beta_rad / |d_x_land|     (elevation plane)
  F_focus_az = R_eff * d_az_rad  / |d_y_land_az|   (azimuth   plane)
  F_focus_3d = sqrt(F_focus_el * F_focus_az)        (geometric mean)
  A_tube_3d  = |d_x_land| * |d_y_land_az|          (elliptic footprint [km^2])

The elevation sensitivity comes from adjacent fan rays (beta_i, beta_{i+1}).
The azimuth sensitivity is measured by perturbing az by +/- delta_az_deg and
tracking the change in y_land: d_y / d_az.

For a free-space (isotropic) medium: F_focus_el = F_focus_az = 1.0  (unfocused).

Reference: Coleman (1997) Radio Sci. 32(5), Eq. 15-17;
           upgrade_plan_3d_scatter.txt Section A3.
"""
import numpy as np
from .ray_tracer_3d import trace_single_ray_3d
from config import FREQ_MHZ, C_MS, TUBE_3D, RT_3D
from utils import to_dBW


class TubeRayTracer3D:
    """
    3-D flux-tube ray tracer for OTH radar / HF point-to-point backscatter.

    Parameters
    ----------
    n3d         : RefractiveIndex3D   electron density + perturbation model
    freq_MHz    : float               operating frequency [MHz]
    tube_params : dict                TUBE_3D config dict (None -> config default)
    rt_params   : dict                RT_3D  config dict (None -> config default)
    """

    def __init__(self, n3d, freq_MHz=FREQ_MHZ, tube_params=None, rt_params=None):
        self.n3d         = n3d
        self.freq_MHz    = float(freq_MHz)
        self.tube_params = tube_params if tube_params is not None else TUBE_3D
        self.rt_params   = rt_params   if rt_params   is not None else RT_3D
        self.lambda_m    = C_MS / (self.freq_MHz * 1e6)

    # ── Private helpers ───────────────────────────────────────────────────────

    def _trace(self, tx_pos, beta_deg, az_deg=0.0):
        """Thin wrapper around trace_single_ray_3d."""
        return trace_single_ray_3d(
            tx_pos, float(beta_deg), float(az_deg),
            self.n3d, freq_MHz=self.freq_MHz, rt_params=self.rt_params)

    # ── Fan shooting ──────────────────────────────────────────────────────────

    def shoot_fan(self, tx_pos, az_deg=0.0):
        """
        Shoot elevation fan at fixed azimuth az_deg.
        Returns list of ray dicts (trace_single_ray_3d output).
        """
        tp    = self.tube_params
        rt    = self.rt_params
        db    = float(tp.get('delta_beta_deg', 0.5))
        n_max = int(tp.get('n_tube_rays', 80))
        b_min = float(rt.get('beta_min', 5.0))
        b_max = float(rt.get('beta_max', 85.0))

        betas = np.arange(b_min, b_max + db * 0.5, db)[:n_max]
        return [self._trace(tx_pos, float(b), az_deg) for b in betas]

    # ── Tube geometry ─────────────────────────────────────────────────────────

    def compute_tubes(self, fan_rays, x_tgt, tx_pos):
        """
        Form 3-D flux tubes from adjacent elevation fan rays near x_tgt.

        For each adjacent pair (r_i, r_{i+1}) with landing near x_tgt:
          - d_x = |x_land_{i+1} - x_land_i|         (elevation footprint width)
          - Shoot 2 extra rays at az = +/-delta_az to measure d_y_land/d_az
          - d_y = |y_land(+daz) - y_land(-daz)|       (azimuth footprint width)
          - A_tube_3d = d_x * d_y                      [km^2]

        Returns list of tube dicts.
        """
        tp      = self.tube_params
        tol_km  = float(tp.get('x_tgt_tol_km', 80.0))
        d_az    = float(tp.get('delta_az_deg', 0.5))
        T_pulse = float(tp.get('T_pulse_ms', 0.5))
        db_rad  = np.radians(float(tp.get('delta_beta_deg', 0.5)))
        daz_rad = np.radians(d_az)

        tubes = []
        for i in range(len(fan_rays) - 1):
            r0 = fan_rays[i]
            r1 = fan_rays[i + 1]
            x0 = float(r0['x_land_km'])
            x1 = float(r1['x_land_km'])

            if abs(x0 - x_tgt) > tol_km and abs(x1 - x_tgt) > tol_km:
                continue

            d_x = abs(x1 - x0)
            if d_x < 1e-6:
                continue

            R_eff_km   = float(r0['group_path_km'])
            F_focus_el = max(R_eff_km * db_rad / d_x, 1e-9)

            # Azimuth sensitivity: d(y_land) / d(az) at center elevation
            beta_c = float(r0['beta_deg_tx'])
            r_p    = self._trace(tx_pos, beta_c,  d_az)
            r_m    = self._trace(tx_pos, beta_c, -d_az)
            d_y    = abs(float(r_p['y_land_km']) - float(r_m['y_land_km']))

            if d_y < 1e-6:
                # Isotropic fallback: free-space limit dy/daz = R_eff
                d_y = R_eff_km * 2.0 * daz_rad

            dydaz_km_rad = d_y / (2.0 * daz_rad)   # [km/rad]
            F_focus_az   = max(R_eff_km / dydaz_km_rad, 1e-9)
            F_focus_3d   = float(np.sqrt(F_focus_el * F_focus_az))
            A_tube_3d    = d_x * d_y                # [km^2]

            tau_spread_ms = abs(float(r1['tau_ms']) - float(r0['tau_ms']))

            tubes.append({
                'beta_deg'      : beta_c,
                'x_land_km'     : 0.5 * (x0 + x1),
                'y_land_km'     : float(r0['y_land_km']),
                'tau_ms'        : float(r0['tau_ms']),
                'group_path_km' : R_eff_km,
                'h_reflect_km'  : float(r0['h_reflect_km']),
                'd_x_km'        : d_x,
                'd_y_km'        : d_y,
                'A_tube_km2'    : A_tube_3d,
                'F_focus_el'    : F_focus_el,
                'F_focus_az'    : F_focus_az,
                'F_focus_3d'    : F_focus_3d,
                'tau_spread_ms' : tau_spread_ms,
                'ray'           : r0,
            })

        return tubes

    # ── Backscatter power ─────────────────────────────────────────────────────

    def backscatter_power_W(self, Pt_W, Gt, Gr, P_one_km, sigma0, A_tube_km2):
        """
        Coleman (1997) radar equation with 3-D tube area.

        Pr = Pt*Gt*Gr*lambda^2*sigma0*A_tube / ((4*pi)^3 * P_one^4)
        """
        P_m  = float(P_one_km) * 1e3
        A_m2 = float(A_tube_km2) * 1e6
        if P_m <= 0.0 or A_m2 <= 0.0:
            return 0.0
        num   = Pt_W * Gt * Gr * (self.lambda_m**2) * sigma0 * A_m2
        denom = ((4.0 * np.pi)**3) * (P_m**4)
        return float(num / max(denom, 1e-300))

    # ── Main pipeline ─────────────────────────────────────────────────────────

    def compute(self, tx_pos, x_tgt, Pt_W, Gt, Gr, sigma0):
        """
        Full 3-D flux-tube backscatter pipeline.

        Parameters
        ----------
        tx_pos  : (x0, z0) [km]
        x_tgt   : target horizontal distance [km]
        Pt_W    : transmit power [W]
        Gt, Gr  : TX/RX antenna gains (linear)
        sigma0  : ground normalised RCS (linear)

        Returns list of mode dicts with fields:
          tau_ms, tau_2way_ms, tau_spread_ms, Pr_W, Pr_dBW,
          h_reflect_km, group_path_km, beta_deg,
          A_tube_km2, F_focus_el, F_focus_az, F_focus_3d
        """
        fan_rays = self.shoot_fan(tx_pos)
        tubes    = self.compute_tubes(fan_rays, x_tgt, tx_pos)

        results = []
        T_pulse = float(self.tube_params.get('T_pulse_ms', 0.5))

        for t in tubes:
            ts  = float(t['tau_spread_ms'])
            f_p = min(T_pulse / ts, 1.0) if ts > T_pulse and ts > 0.0 else 1.0

            Pr_W = self.backscatter_power_W(
                Pt_W, Gt, Gr,
                t['group_path_km'], sigma0, t['A_tube_km2'])
            Pr_W *= f_p

            results.append({
                'tau_ms'        : t['tau_ms'],
                'tau_2way_ms'   : 2.0 * t['tau_ms'],
                'tau_spread_ms' : t['tau_spread_ms'],
                'Pr_W'          : Pr_W,
                'Pr_dBW'        : to_dBW(Pr_W),
                'h_reflect_km'  : t['h_reflect_km'],
                'group_path_km' : t['group_path_km'],
                'beta_deg'      : t['beta_deg'],
                'A_tube_km2'    : t['A_tube_km2'],
                'F_focus_el'    : t['F_focus_el'],
                'F_focus_az'    : t['F_focus_az'],
                'F_focus_3d'    : t['F_focus_3d'],
            })

        return results
