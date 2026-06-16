"""
Phase 1 verification: OTH radar geometry adaptation.

Checks:
  1. tau_2way_ms = 2 * tau_ms  (error < 0.001 ms)
  2. tau_2way_ms is None in communication mode
  3. Radar Pr << Friis Pr  (difference > 40 dB expected at ~1000 km, 5 m^2)
  4. phi_deg == LINK_BEARING_DEG in radar mode
  5. Mode count identical between comm and radar

Run:
    conda run -n pytorch_cpu python tests/_verify_radar.py
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import config as cfg
from models.hybrid_model import HybridPropagationModel


def _make_iono():
    return {
        'iri_params':    {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':    {**cfg.TID,    'enable': False},
        'es_params':     {**cfg.ES,     'enable': False},
        'bubble_params': {**cfg.BUBBLE, 'enable': False},
    }


def _make_comm_rp():
    return {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W,
            'Gt': cfg.GT, 'Gr': cfg.GR}


def _make_radar_rp():
    return {**_make_comm_rp(),
            'sigma_rcs_m2': cfg.RADAR['sigma_rcs_m2'],
            'bearing_deg':  cfg.LINK_BEARING_DEG}


def run():
    # Use reduced P2P settings for speed in the test
    p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
    n_fail = 0

    print("Building communication model (Friis) ...")
    comm = HybridPropagationModel(_make_iono(), _make_comm_rp(), radar_mode=False)
    comm_modes, _, _, _ = comm.compute(cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)
    print("  {} modes found".format(len(comm_modes)))

    print("Building radar model (radar equation) ...")
    radar = HybridPropagationModel(_make_iono(), _make_radar_rp(), radar_mode=True)
    radar_modes, _, _, _ = radar.compute(cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)
    print("  {} modes found".format(len(radar_modes)))

    # ── Check 1: tau_2way_ms = 2 * tau_ms ────────────────────────────────────
    print("\n[1] tau_2way_ms = 2 * tau_ms (tolerance 0.001 ms)")
    for m in radar_modes:
        diff = abs(m['tau_2way_ms'] - 2.0 * m['tau_ms'])
        ok = diff < 0.001
        print("    {:16s}  tau={:.3f}  tau_2way={:.3f}  diff={:.5f}  {}".format(
            m['label'], m['tau_ms'], m['tau_2way_ms'], diff,
            'OK' if ok else 'FAIL'))
        if not ok:
            n_fail += 1

    # ── Check 2: tau_2way_ms is None in comm mode ─────────────────────────────
    print("\n[2] tau_2way_ms is None in communication mode")
    for m in comm_modes:
        ok = m.get('tau_2way_ms') is None
        print("    {:16s}  tau_2way={}  {}".format(
            m['label'], m.get('tau_2way_ms'), 'OK' if ok else 'FAIL'))
        if not ok:
            n_fail += 1

    # ── Check 3: radar power << Friis power ───────────────────────────────────
    print("\n[3] Radar Pr << Friis Pr (expect > 40 dB difference)")
    for rm in radar_modes:
        cm_match = min(comm_modes, key=lambda m: abs(m['tau_ms'] - rm['tau_ms']))
        if abs(cm_match['tau_ms'] - rm['tau_ms']) > 0.5:
            continue
        diff_dB = (cm_match['Pr_dBW'] - rm['Pr_dBW'])   # positive = radar weaker
        ok = diff_dB > 40.0
        print("    {:16s}  Friis={:.1f} dBW  Radar={:.1f} dBW  diff={:.1f} dB  {}".format(
            rm['label'], cm_match['Pr_dBW'], rm['Pr_dBW'], diff_dB,
            'OK' if ok else 'FAIL'))
        if not ok:
            n_fail += 1

    # ── Check 4: phi_deg = LINK_BEARING_DEG ──────────────────────────────────
    print("\n[4] phi_deg = LINK_BEARING_DEG ({})".format(cfg.LINK_BEARING_DEG))
    for m in radar_modes:
        ok = abs(m['phi_deg'] - cfg.LINK_BEARING_DEG) < 0.01
        print("    {:16s}  phi={:.2f}  {}".format(
            m['label'], m['phi_deg'], 'OK' if ok else 'FAIL'))
        if not ok:
            n_fail += 1

    # ── Check 5: mode count identical ────────────────────────────────────────
    print("\n[5] Mode count  comm={:d}  radar={:d}  {}".format(
        len(comm_modes), len(radar_modes),
        'OK' if len(comm_modes) == len(radar_modes) else 'WARN (expected equal)'))

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    if n_fail == 0:
        print("Phase 1 verification PASSED ({} checks).".format(
            len(radar_modes) * 2 + len(comm_modes) + len(radar_modes)))
    else:
        print("{} check(s) FAILED.".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
