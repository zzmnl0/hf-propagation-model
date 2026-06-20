"""
Phase 6 verification: 3-D ray tracing + Born/Rytov scatter theory.

Checks:
  1. A1 3-D ray: trace_single_ray_3d with az=0 matches 2-D tau within 0.1 ms
  2. A1 3-D TID: az=30 deg TID deflects ray (y_land != 0)
  3. A2 3-D P2P: find_all_rays_3d finds >= 1 mode with landing error < 5 km
  4. C1 Born:    sigma_v > 0 and physically in range 1e-20 to 1e-5 m^-1
  5. C2 Rytov:   S4 in [0, 0.6] for SCATTER params; sigma_tau_ms > 0
  Bonus C0:       scatter_phase_variance increases with Cs_rel

Run:
    conda run -n pytorch_cpu python tests/_verify_phase6.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from config import BG_X, BG_Z, TX_POS, RX_POS, FREQ_MHZ, SCATTER
from models.ionosphere_model import IonosphereModel
from models.ray_tracer import RefractiveIndex, trace_single_ray
from models.ray_tracer_3d import RefractiveIndex3D, trace_single_ray_3d
from models.point_to_point_3d import find_all_rays_3d
from models.point_to_point import find_all_rays_p2p
from models.scatter_selector import scatter_phase_variance, select_scatter_method
from models.scatter_born import born_sigma_v, born_scatter
from models.scatter_rytov import rytov_full


def run():
    n_fail = 0

    print("=" * 55)
    print("  Phase 6 verification: 3-D ray tracing + scatter")
    print("=" * 55)

    # ── Build baseline ionosphere ─────────────────────────────────────────────
    iono = IonosphereModel()
    Ne, _ = iono.build_Ne_field(BG_X, BG_Z)

    # 2-D reference
    nm_2d = RefractiveIndex(Ne, BG_X, BG_Z, FREQ_MHZ)
    # 3-D without Earth-flattening for tau comparison (same ionosphere as 2-D)
    nm_3d_noflat = RefractiveIndex3D(Ne, BG_X, BG_Z, FREQ_MHZ, earth_flat=False)
    # 3-D with Earth-flattening for all other 3-D checks
    nm_3d = RefractiveIndex3D(Ne, BG_X, BG_Z, FREQ_MHZ, earth_flat=True)

    # ── Check 1: A1 3-D az=0 ray tau reduces to 2-D single-ray tau ──────────
    # Both traces use the same ionosphere (no Earth-flattening) at the same
    # launch elevation, so 3-D with az=0 must reproduce the 2-D result.
    print("\nCheck 1: A1 3-D ray (az=0, noflat) tau == 2-D single-ray tau within 0.1 ms")

    beta_test = 30.0     # representative F-layer elevation
    ray2d_s   = trace_single_ray(TX_POS, beta_test, nm_2d, FREQ_MHZ)
    ray3d_s   = trace_single_ray_3d(TX_POS, beta_test, 0.0, nm_3d_noflat, FREQ_MHZ)
    tau_2d_s  = ray2d_s['tau_ms']
    tau_3d_s  = ray3d_s['tau_ms']
    diff      = abs(tau_3d_s - tau_2d_s)
    ok1       = diff < 0.1
    print("  2-D tau={:.3f}ms  3-D tau={:.3f}ms  diff={:.4f}ms  x2D={:.1f}km x3D={:.1f}km  {}".format(
        tau_2d_s, tau_3d_s, diff,
        ray2d_s['trajectory'][-1][0], ray3d_s['x_land_km'],
        'OK' if ok1 else 'FAIL'))
    if not ok1:
        n_fail += 1

    # ── Check 2: A1 3-D TID cross-path deflection ────────────────────────────
    print("\nCheck 2: A1 3-D TID az=30 deg -> y_land != 0")

    tid_3d = {
        **cfg.TID,
        'enable':           True,
        'n_components':     1,
        'az_deg_list':      [30.0],        # 30 deg from link -> cross-path component
        'amplitude_list':   [0.15],
        'period_s_list':    [2400.0],
        'lambda_h_km_list': [300.0],
        'link_bearing_deg': cfg.LINK_BEARING_DEG,
    }
    nm_3d_tid = RefractiveIndex3D(Ne, BG_X, BG_Z, FREQ_MHZ,
                                   tid_params=tid_3d, earth_flat=True)

    # Use baseline elevation (30 deg) with az=0; TID should deflect ray laterally
    ray_notid = trace_single_ray_3d(TX_POS, 30.0, 0.0, nm_3d, FREQ_MHZ)
    ray_tid   = trace_single_ray_3d(TX_POS, 30.0, 0.0, nm_3d_tid, FREQ_MHZ)
    y_notid   = abs(ray_notid['y_land_km'])
    y_tid     = abs(ray_tid['y_land_km'])
    ok2       = y_tid > 0.1 or y_notid > 0.1   # at least one has lateral component
    # A more specific check: the TID with az=30 changes y_land
    dy = abs(y_tid - y_notid)
    ok2 = dy > 0.05  # at least 50 m lateral shift
    print("  y_land (no TID)={:.3f}km  y_land (TID az=30)={:.3f}km  "
          "dy={:.3f}km  {}".format(y_notid, y_tid, dy, 'OK' if ok2 else 'FAIL'))
    if not ok2:
        n_fail += 1

    # ── Check 3: A2 3-D P2P finds modes (landing within 5 km) ───────────────
    print("\nCheck 3: A2 3-D P2P -> >= 1 mode, landing error < 5 km")

    p2p_3d = find_all_rays_3d(TX_POS, RX_POS, nm_3d, FREQ_MHZ,
                               p2p_params=cfg.P2P_3D)
    if not p2p_3d:
        print("  FAIL: no 3-D P2P modes found")
        n_fail += 1
    else:
        x_err = abs(p2p_3d[0]['x_land_km'] - float(RX_POS[0]))
        y_err = abs(p2p_3d[0]['y_land_km'])
        err   = np.hypot(x_err, y_err)
        ok3   = err < 5.0
        print("  n_modes={} best_tau={:.3f}ms landing_err={:.2f}km  {}".format(
            len(p2p_3d), p2p_3d[0]['tau_ms'], err, 'OK' if ok3 else 'FAIL'))
        for r in p2p_3d:
            print("    {} tau={:.3f}ms h={:.0f}km az_defl={:.2f}deg".format(
                r['label'], r['tau_ms'], r['h_reflect_km'],
                r['azimuth_deflect_deg']))
        if not ok3:
            n_fail += 1

    # ── Check 4: C1 Born sigma_v physically reasonable ───────────────────────
    print("\nCheck 4: C1 Born sigma_v in [1e-20, 1e-5] m^-1")

    Ne0   = 5e10   # moderate F2-layer [m^-3]
    sv    = born_sigma_v(FREQ_MHZ,
                         SCATTER['Cs_rel'],
                         SCATTER['p'],
                         SCATTER['L_outer_km'],
                         Ne0)
    ok4   = 1e-20 <= sv <= 1e-5
    print("  sigma_v={:.3e} m^-1  (Cs={}, p={}, L_outer={} km, Ne0={:.1e})  {}".format(
        sv, SCATTER['Cs_rel'], SCATTER['p'], SCATTER['L_outer_km'], Ne0,
        'OK' if ok4 else 'FAIL'))
    if not ok4:
        n_fail += 1

    # Also check Bragg wavenumber value at 10 MHz
    k0  = 2*np.pi * FREQ_MHZ*1e6 / 3e8
    kb  = 2*k0
    print("  Bragg wavenumber 2k0 = {:.4f} rad/m  (probed scale {:.1f} m)".format(
        kb, 2*np.pi/kb))

    # ── Check 5: C2 Rytov S4 and sigma_tau physically reasonable ─────────────
    print("\nCheck 5: C2 Rytov S4 in [0, 0.6] and sigma_tau > 0")

    ryt = rytov_full(FREQ_MHZ,
                     SCATTER['Cs_rel'],
                     SCATTER['p'],
                     SCATTER['L_outer_km'],
                     Ne0,
                     SCATTER['z_eff_km'],
                     SCATTER['DeltaL_km'])
    ok5a = 0.0 <= ryt['S4'] <= 0.6
    ok5b = ryt['sigma_tau_ms'] > 0.0
    print("  S4={:.4f}  sigma_phi={:.4f} rad  sigma_tau={:.4e} ms  "
          "r_F={:.2f} km  {}".format(
        ryt['S4'], ryt['sigma_phi_rad'], ryt['sigma_tau_ms'],
        ryt['fresnel_r_km'], 'OK' if (ok5a and ok5b) else 'FAIL'))
    if not (ok5a and ok5b):
        n_fail += 1

    # ── Bonus C0: sigma_phi^2 increases with Cs_rel ──────────────────────────
    print("\nBonus C0: scatter_phase_variance increases with Cs_rel")
    sp1 = scatter_phase_variance(Ne0, 0.01,  100.0, 10.0, FREQ_MHZ)
    sp2 = scatter_phase_variance(Ne0, 0.10,  100.0, 10.0, FREQ_MHZ)
    sp3 = scatter_phase_variance(Ne0, 0.30,  100.0, 10.0, FREQ_MHZ)
    okC = sp1 < sp2 < sp3
    print("  Cs=0.01->{:.3e}  Cs=0.10->{:.3e}  Cs=0.30->{:.3e}  {}".format(
        sp1, sp2, sp3, 'OK' if okC else 'FAIL'))
    print("  Methods: {}/{}/{}".format(
        select_scatter_method(sp1), select_scatter_method(sp2),
        select_scatter_method(sp3)))
    if not okC:
        n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
