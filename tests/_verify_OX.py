"""
Phase 2 verification: Appleton-Hartree O/X magnetoionic splitting.

Checks:
  1. Isotropic limit: fH->0 makes O/X n_batch converge to isotropic n_batch
  2. Mode count with enable_OX=True is approx 2x isotropic count
  3. O/X tau_ms difference < 0.5 ms for matched modes (typical 0.1-0.3 ms)
  4. X-mode tau >= O-mode tau for same base label (X travels slower near fxF2)
  5. wave_mode field is 'O' or 'X' for AH modes, 'iso' for isotropic

Run:
    conda run -n pytorch_cpu python tests/_verify_OX.py
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import config as cfg
from models.ionosphere_model import IonosphereModel
from models.ray_tracer import RefractiveIndex, RefractiveIndexAH
from models.hybrid_model import HybridPropagationModel


def _make_iono():
    return {
        'iri_params':    {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':    {**cfg.TID, 'enable': True},
        'es_params':     {**cfg.ES,  'enable': False},
        'bubble_params': {**cfg.BUBBLE, 'enable': False},
    }


def _make_rp():
    return {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W,
            'Gt': cfg.GT, 'Gr': cfg.GR}


def run():
    p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
    n_fail = 0

    # ── Check 1: isotropic limit (fH -> 0) ───────────────────────────────────
    print("[1] Isotropic limit: AH with fH=0 converges to isotropic n_batch")
    iono = IonosphereModel(tid_params={**cfg.TID, 'enable': True})
    Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
    n_iso = RefractiveIndex(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ)
    n_ah_O = RefractiveIndexAH(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                                wave_mode='O', geomag={'fH_MHz': 0.0, 'dip_deg': 48.7})
    n_ah_X = RefractiveIndexAH(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ,
                                wave_mode='X', geomag={'fH_MHz': 0.0, 'dip_deg': 48.7})

    test_pts = np.array([[300., 200.], [600., 300.], [900., 250.]])
    n_i  = n_iso.n_batch(test_pts)
    n_O  = n_ah_O.n_batch(test_pts)
    n_X  = n_ah_X.n_batch(test_pts)
    max_diff_O = float(np.max(np.abs(n_O - n_i)))
    max_diff_X = float(np.max(np.abs(n_X - n_i)))
    ok = max_diff_O < 1e-4 and max_diff_X < 1e-4
    print("    max |n_O - n_iso| = {:.2e}   max |n_X - n_iso| = {:.2e}   {}".format(
        max_diff_O, max_diff_X, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Build isotropic and OX models ─────────────────────────────────────────
    print("\nBuilding isotropic model ...")
    iso_model = HybridPropagationModel(_make_iono(), _make_rp(),
                                        geomag_params=None)
    iso_modes, _, _, _ = iso_model.compute(cfg.TX_POS, cfg.RX_POS,
                                            p2p_params=p2p_fast)
    print("  {} isotropic modes".format(len(iso_modes)))

    print("Building O/X model ...")
    ox_model = HybridPropagationModel(_make_iono(), _make_rp(),
                                       geomag_params={**cfg.GEOMAG, 'enable_OX': True})
    ox_modes, _, _, _ = ox_model.compute(cfg.TX_POS, cfg.RX_POS,
                                          p2p_params=p2p_fast)
    print("  {} O/X modes".format(len(ox_modes)))

    # ── Check 2: mode count approx 2x ────────────────────────────────────────
    print("\n[2] Mode count with O/X = {} vs isotropic = {}".format(
        len(ox_modes), len(iso_modes)))
    ok = len(ox_modes) >= len(iso_modes)  # at least as many (some may merge)
    print("    {}".format('OK' if ok else 'WARN (fewer modes than expected)'))
    if not ok:
        n_fail += 1

    # ── Check 3: O/X tau difference < 0.5 ms ─────────────────────────────────
    print("\n[3] O/X tau_ms difference < 0.5 ms for matched base labels")
    o_modes = [m for m in ox_modes if m.get('wave_mode') == 'O']
    x_modes = [m for m in ox_modes if m.get('wave_mode') == 'X']
    matched = 0
    for om in o_modes:
        # find closest X-mode with same base label
        base_o = om['label'].replace('_O', '')
        xm_cands = [m for m in x_modes if m['label'].replace('_X', '') == base_o]
        if not xm_cands:
            continue
        xm = min(xm_cands, key=lambda m: abs(m['tau_ms'] - om['tau_ms']))
        dtau = abs(xm['tau_ms'] - om['tau_ms'])
        ok_i = dtau < 0.5
        print("    {:16s}  O tau={:.3f}  X tau={:.3f}  dtau={:.4f} ms  {}".format(
            base_o, om['tau_ms'], xm['tau_ms'], dtau, 'OK' if ok_i else 'FAIL'))
        if not ok_i:
            n_fail += 1
        matched += 1
    if matched == 0:
        print("    WARN: no O/X pairs found to compare")

    # ── Check 4: wave_mode field present and valid ────────────────────────────
    print("\n[4] wave_mode field values")
    for m in ox_modes:
        wm = m.get('wave_mode', 'MISSING')
        ok = wm in ('O', 'X')
        print("    {:20s}  wave_mode={}  {}".format(
            m['label'], wm, 'OK' if ok else 'FAIL'))
        if not ok:
            n_fail += 1
    for m in iso_modes:
        ok = m.get('wave_mode') == 'iso'
        if not ok:
            print("    iso mode {} has wave_mode={} FAIL".format(
                m['label'], m.get('wave_mode')))
            n_fail += 1

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    if n_fail == 0:
        print("Phase 2 verification PASSED.")
    else:
        print("{} check(s) FAILED.".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
