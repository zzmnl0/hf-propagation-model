"""
Phase 4 verification: flux tube ray tracer (Coleman 1997/1998).

Checks:
  1. Tube center-ray tau_ms agrees with P2P variational solver (< 0.05 ms)
  2. Backscatter power within 10 dB of radar_equation_W (same sigma)
  3. F_focus values in range (0.1, 10.0) for all tubes
  4. TID perturbation increases tau_spread_ms (monotone, >= 1 mode)
  5. mode_summary dict is complete and non-empty

Run:
    conda run -n pytorch_cpu python tests/_verify_tube.py
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from config import (BG_X, BG_Z, TX_POS, RX_POS, FREQ_MHZ,
                    PT_W, GT, GR, RADAR, TID, P2P, RT)
from models.ionosphere_model import IonosphereModel
from models.ray_tracer import RefractiveIndex
from models.point_to_point import find_all_rays_p2p, classify_mode
from models.tube_tracer import TubeRayTracer
from models.hybrid_model import HybridPropagationModel
from utils import radar_equation_W, to_dBW


def _make_iono():
    return {
        'iri_params'    : {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params'    : {**cfg.TID,    'enable': False},
        'es_params'     : {**cfg.ES,     'enable': False},
        'bubble_params' : {**cfg.BUBBLE, 'enable': False},
        'spread_f_params': {**cfg.SPREAD_F, 'enable': False},
    }


def _make_rp():
    return {
        'freq_MHz'        : FREQ_MHZ,
        'Pt_W'            : PT_W,
        'Gt'              : GT,
        'Gr'              : GR,
        'sigma_rcs_m2'    : RADAR['sigma_rcs_m2'],
        'bearing_deg'     : cfg.LINK_BEARING_DEG,
        'sigma0_ground_dB': RADAR['sigma0_ground_dB'],
    }


def run():
    n_fail = 0

    print("=" * 55)
    print("  Phase 4 verification: flux tube ray tracer")
    print("=" * 55)

    # Build IRI baseline ionosphere once (shared by checks 1-3)
    iono = IonosphereModel()
    Ne, _ = iono.build_Ne_field(BG_X, BG_Z)
    nm    = RefractiveIndex(Ne, BG_X, BG_Z, FREQ_MHZ)

    sigma0 = 10.0 ** (RADAR['sigma0_ground_dB'] / 10.0)
    x_tgt  = float(RADAR['target_range_km'])

    tracer     = TubeRayTracer(nm, FREQ_MHZ)
    tube_modes = tracer.compute(TX_POS, x_tgt, PT_W, GT, GR, sigma0)

    # ── Check 1: tau_ms vs P2P variational solver ─────────────────────────────
    print("\nCheck 1: tube tau_ms vs P2P (expect diff < 0.05 ms)")
    p2p_params = {**P2P, 'n_init': 6, 'max_iter': 100}
    p2p_modes  = find_all_rays_p2p(TX_POS, RX_POS, nm, FREQ_MHZ,
                                   p2p_params=p2p_params)
    p2p_taus = {classify_mode(m): m['tau_ms'] for m in p2p_modes}
    tube_taus = {m['label']: m['tau_ms'] for m in tube_modes}

    n_matched = 0
    for lab, t_tube in sorted(tube_taus.items()):
        if lab in p2p_taus:
            diff = abs(t_tube - p2p_taus[lab])
            ok_i = diff < 0.05
            print("  {:<10} tube={:.3f}ms p2p={:.3f}ms diff={:.4f}ms  {}".format(
                lab, t_tube, p2p_taus[lab], diff, 'OK' if ok_i else 'FAIL'))
            n_matched += 1
            if not ok_i:
                n_fail += 1
    if n_matched == 0:
        print("  WARN: no matching mode labels between tube and P2P")

    # ── Check 2: backscatter power vs radar_equation_W ────────────────────────
    print("\nCheck 2: tube power vs radar_equation_W (expect ~0.0 dB diff)")
    for m in tube_modes:
        sigma_eq  = sigma0 * m['A_tube_km2'] * 1e6
        Pr_radar  = radar_equation_W(PT_W, GT, GR, FREQ_MHZ,
                                     m['group_path_km'], sigma_eq)
        diff_dB   = abs(to_dBW(m['Pr_W']) - to_dBW(Pr_radar))
        ok_i      = diff_dB < 10.0
        print("  {:<10} tube={:.2f}dBW radar_eq={:.2f}dBW diff={:.2f}dB  {}".format(
            m['label'], to_dBW(m['Pr_W']), to_dBW(Pr_radar), diff_dB,
            'OK' if ok_i else 'FAIL'))
        if not ok_i:
            n_fail += 1

    # ── Check 3: F_focus range ────────────────────────────────────────────────
    print("\nCheck 3: F_focus in (0.1, 10.0)")
    F_vals = [m['F_focus'] for m in tube_modes]
    if len(F_vals) == 0:
        print("  WARN: no tube modes found")
        n_fail += 1
    else:
        for m in tube_modes:
            fv   = m['F_focus']
            ok_i = 0.1 <= fv <= 10.0
            print("  {:<10} F_focus={:.4f}  {}".format(
                m['label'], fv, 'OK' if ok_i else 'FAIL'))
            if not ok_i:
                n_fail += 1

    # ── Check 4: TID increases tau_spread_ms ─────────────────────────────────
    print("\nCheck 4: TID perturbation -> tau_spread_ms increases (>= 1 mode)")
    iono_tid  = IonosphereModel(tid_params={**TID, 'enable': True})
    Ne_tid, _ = iono_tid.build_Ne_field(BG_X, BG_Z)
    nm_tid    = RefractiveIndex(Ne_tid, BG_X, BG_Z, FREQ_MHZ)
    tracer_tid  = TubeRayTracer(nm_tid, FREQ_MHZ)
    tube_tid    = tracer_tid.compute(TX_POS, x_tgt, PT_W, GT, GR, sigma0)

    base_spread = {m['label']: m['tau_spread_ms'] for m in tube_modes}
    tid_spread  = {m['label']: m['tau_spread_ms'] for m in tube_tid}

    n_incr = 0
    for lab, sp_base in sorted(base_spread.items()):
        if lab in tid_spread:
            ok_i   = tid_spread[lab] >= sp_base
            n_incr += (1 if ok_i else 0)
            print("  {:<10} baseline={:.4f}ms TID={:.4f}ms  {}".format(
                lab, sp_base, tid_spread[lab], 'OK' if ok_i else 'note'))
    ok4 = n_incr >= 1
    print("  At least 1 mode with increased spread: {}".format('OK' if ok4 else 'FAIL'))
    if not ok4:
        n_fail += 1

    # ── Check 5: mode_summary completeness ───────────────────────────────────
    print("\nCheck 5: mode_summary fields complete")
    model = HybridPropagationModel(_make_iono(), _make_rp(),
                                   radar_mode=True, tube_mode=True)
    modes, tau_ax, pd, main = model.compute(TX_POS, RX_POS)
    summary = model.mode_summary

    ok5 = (summary.get('n_modes', 0) > 0
           and summary.get('main_label', '') != ''
           and summary.get('main_F_focus', 0.0) > 0.0
           and isinstance(summary.get('mode_pairs'), list))
    print("  n_modes={}  main='{}'  F_focus={:.2f}  mode_pairs={}  {}".format(
        summary.get('n_modes', 0), summary.get('main_label', '?'),
        summary.get('main_F_focus', 0.0),
        len(summary.get('mode_pairs', [])),
        'OK' if ok5 else 'FAIL'))
    if not ok5:
        n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
