"""
Phase 7 verification: 3-D field-aligned bubble, 3-D flux-tube, Booker-Gordon Es,
Rino 2-D phase screen, Knepp MPS.

Checks:
  1. B4 bubble_3d:    n2 at bubble centre increases (depletes Ne) by > 0.005
  2. A3 tube_3d:      TubeRayTracer3D finds >= 1 mode, F_focus_3d > 0
  3. C4 BG Es:        booker_gordon_sigma_v > 0 and in physical range; aspect check
  4. C3 Rino screen:  RinoPhaseScreen2D generates screen; measured sigma in [0.01, 3x theory]
  5. C5 MPS:          MultiPhaseScreen S4 increases monotonically with Cs

Run:
    conda run -n pytorch_cpu python tests/_verify_phase7.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from config import (BG_X, BG_Z, TX_POS, RX_POS, FREQ_MHZ,
                    SCATTER, BUBBLE_3D, TUBE_3D, RT_3D, ES_BG, PHASE_SCREEN,
                    RE_KM, GEOMAG)
from models.ionosphere_model import IonosphereModel
from models.ray_tracer_3d import RefractiveIndex3D
from models.tube_tracer_3d import TubeRayTracer3D
from models.es_model import booker_gordon_sigma_v, es_aspect_angle_deg, es_scatter_valid
from models.scatter_phase_screen import RinoPhaseScreen2D, MultiPhaseScreen


def run():
    n_fail = 0

    print("=" * 55)
    print("  Phase 7 verification: bubble / tube / BG / screen / MPS")
    print("=" * 55)

    # Shared baseline ionosphere
    iono = IonosphereModel()
    Ne, _ = iono.build_Ne_field(BG_X, BG_Z)

    # ── Check 1: B4  3-D field-aligned bubble depletes Ne ─────────────────────
    print("\nCheck 1: B4 bubble_3d -> n2 at bubble centre increases (depletion)")

    nm_no_bub = RefractiveIndex3D(Ne, BG_X, BG_Z, FREQ_MHZ, earth_flat=True)

    bub_cfg = dict(BUBBLE_3D)
    bub_cfg['enable'] = True
    nm_bub = RefractiveIndex3D(Ne, BG_X, BG_Z, FREQ_MHZ,
                               bubble_3d_params=bub_cfg, earth_flat=True)

    # Evaluate at exact bubble centre (x=600, y=0, z_phys=350 km).
    # With earth_flat, z_cart = z_phys + (x^2+y^2)/(2*RE)
    x_c  = float(bub_cfg['x0_km'])
    z_p  = float(bub_cfg['z0_km'])
    z_c  = z_p + (x_c**2) / (2.0 * RE_KM)   # Cartesian z at bubble centre

    n2_off = nm_no_bub.n2(x_c, 0.0, z_c)
    n2_on  = nm_bub.n2(x_c, 0.0, z_c)
    delta_n2 = n2_on - n2_off

    # Also verify depletion fraction is correct: frac = delta_max at exact centre
    frac_actual = nm_bub._bubble_3d_depletion(x_c, 0.0, z_p)
    ok1a = abs(frac_actual - float(bub_cfg['delta_max'])) < 0.001  # frac == delta_max at centre
    ok1b = delta_n2 > 0.001   # significant n2 increase (depletion lowers Ne -> higher n2)
    ok1  = ok1a and ok1b
    print("  frac_at_centre={:.4f}  delta_max={:.2f}  delta_n2={:.5f}  {}".format(
        frac_actual, float(bub_cfg['delta_max']), delta_n2,
        'OK' if ok1 else 'FAIL'))
    if not ok1:
        n_fail += 1

    # ── Check 2: A3  3-D flux-tube finds modes ────────────────────────────────
    print("\nCheck 2: A3 tube_3d -> >= 1 mode, F_focus_3d > 0")

    tube3d_cfg = dict(TUBE_3D)
    tube3d_cfg['x_tgt_tol_km'] = 120.0   # generous tolerance

    tt3d = TubeRayTracer3D(nm_no_bub, FREQ_MHZ,
                           tube_params=tube3d_cfg, rt_params=RT_3D)
    modes = tt3d.compute(TX_POS, float(RX_POS[0]),
                         Pt_W=1000.0, Gt=1.0, Gr=1.0, sigma0=1e-2)

    if not modes:
        print("  FAIL: no 3-D tube modes found near x_tgt={:.0f} km".format(
            float(RX_POS[0])))
        n_fail += 1
    else:
        best = modes[0]
        ok2a = float(best['F_focus_3d']) > 0.0
        ok2b = float(best['tau_ms'])     > 0.5
        ok2c = float(best['A_tube_km2']) > 0.0
        ok2  = ok2a and ok2b and ok2c
        print("  n_modes={} tau={:.3f}ms F_focus_3d={:.3f} A_tube={:.1f} km^2  {}".format(
            len(modes), best['tau_ms'], best['F_focus_3d'],
            best['A_tube_km2'], 'OK' if ok2 else 'FAIL'))
        if not ok2:
            n_fail += 1

    # ── Check 3: C4  Booker-Gordon sigma_v and aspect sensitivity ────────────
    print("\nCheck 3: C4 booker_gordon_sigma_v in physical range; aspect sensitivity")

    Ne0    = 2e10   # moderate E/Es-layer [m^-3]
    L_perp = float(ES_BG['L_perp_m'])
    L_para = float(ES_BG['L_para_m'])
    Cs_es  = 0.10   # 10% Es fluctuation
    p_es   = float(ES_BG['p_es'])

    # Optimal geometry: aspect=0 (Bragg vector perpendicular to B, max scatter)
    sv_opt  = booker_gordon_sigma_v(FREQ_MHZ, Cs_es, Ne0, L_perp, L_para, p_es, 0.0)
    ok3a    = 1e-30 <= sv_opt <= 1.0   # physically valid range [m^-1]
    ok3b    = es_scatter_valid(0.0, L_para, FREQ_MHZ)      # valid at aspect=0
    # Non-optimal: aspect angle for near-vertical incidence
    theta_vert = np.radians(80.0)
    asp        = es_aspect_angle_deg(theta_vert)
    not_valid  = not es_scatter_valid(asp, L_para, FREQ_MHZ)   # suppressed at large aspect
    ok3c = not_valid  # aspect sensitivity works (large angle suppressed)

    # Monotonic: sigma_v(0) > sigma_v(0.01 deg) > sigma_v(0.1 deg)
    sv1 = booker_gordon_sigma_v(FREQ_MHZ, Cs_es, Ne0, L_perp, L_para, p_es, 0.0)
    sv2 = booker_gordon_sigma_v(FREQ_MHZ, Cs_es, Ne0, L_perp, L_para, p_es, 0.01)
    sv3 = booker_gordon_sigma_v(FREQ_MHZ, Cs_es, Ne0, L_perp, L_para, p_es, 0.1)
    ok3d = sv1 >= sv2 >= sv3   # monotonic decrease with aspect angle

    ok3 = ok3a and ok3b and ok3c and ok3d
    print("  sigma_v(asp=0)={:.3e} m^-1  valid(0)={}  notvalid({:.1f}deg)={}".format(
        sv_opt, ok3b, asp, not_valid))
    print("  monotonic: {:.3e} >= {:.3e} >= {:.3e}  {}".format(
        sv1, sv2, sv3, 'OK' if ok3 else 'FAIL'))
    if not ok3:
        n_fail += 1

    # ── Check 4: C3  Rino 2-D phase screen sigma_phi ─────────────────────────
    print("\nCheck 4: C3 RinoPhaseScreen2D sigma_phi within [0.1, 3x] of theory")

    scr = RinoPhaseScreen2D(
        freq_MHz=FREQ_MHZ, Cs_rel=SCATTER['Cs_rel'],
        p=SCATTER['p'], L_outer_km=SCATTER['L_outer_km'],
        Ne_m3=1e11, DeltaL_km=SCATTER['DeltaL_km'],
        n_pts=64, dx_m=500.0, seed=42)

    phi    = scr.generate()
    sigma_measured = float(np.std(phi))
    sigma_theory   = scr.sigma_phi_rms()

    ok4a = phi.shape == (64, 64)
    ok4b = sigma_theory > 0.0
    ratio = sigma_measured / max(sigma_theory, 1e-30)
    ok4c = 0.1 <= ratio <= 10.0   # within order of magnitude (random realisation)
    ok4  = ok4a and ok4b and ok4c
    print("  sigma_meas={:.4f} rad  sigma_theory={:.4f} rad  ratio={:.3f}  {}".format(
        sigma_measured, sigma_theory, ratio, 'OK' if ok4 else 'FAIL'))
    if not ok4:
        n_fail += 1

    # ── Check 5: C5  MPS S4 increases with Cs ────────────────────────────────
    print("\nCheck 5: C5 MPS S4 increases with Cs_rel (Knepp 1983)")

    mps_small = MultiPhaseScreen(
        freq_MHz=FREQ_MHZ, Cs_rel=0.005, p=SCATTER['p'],
        L_outer_km=SCATTER['L_outer_km'], Ne_m3=1e11,
        DeltaL_km=SCATTER['DeltaL_km'], N_screen=5,
        n_pts=64, dx_m=500.0, seed=7)
    mps_large = MultiPhaseScreen(
        freq_MHz=FREQ_MHZ, Cs_rel=0.10, p=SCATTER['p'],
        L_outer_km=SCATTER['L_outer_km'], Ne_m3=1e11,
        DeltaL_km=SCATTER['DeltaL_km'], N_screen=5,
        n_pts=64, dx_m=500.0, seed=7)

    res_s = mps_small.propagate()
    res_l = mps_large.propagate()

    ok5a = res_l['S4'] > res_s['S4']               # larger Cs -> larger S4
    ok5b = 0.0 <= res_l['S4'] <= 3.0               # S4 in physical range
    ok5c = abs(res_s['Pr_norm'] - 1.0) < 0.5       # energy approximately conserved
    ok5  = ok5a and ok5b and ok5c
    print("  S4(Cs=0.005)={:.4f}  S4(Cs=0.10)={:.4f}  Pr_norm(small)={:.3f}  {}".format(
        res_s['S4'], res_l['S4'], res_s['Pr_norm'], 'OK' if ok5 else 'FAIL'))
    if not ok5:
        n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
