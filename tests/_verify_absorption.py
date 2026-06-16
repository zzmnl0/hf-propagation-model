"""
Phase 3 verification: D-layer absorption.

Checks:
  1. A_dB > 0 for typical daytime parameters
  2. A_dB = 0 at nighttime (chi = 90 deg)
  3. A_dB decreases with higher frequency  (1/f^2 dependence)
  4. A_dB increases at lower elevation  (1/sin(beta) dependence)
  5. Absorption reduces Pr_W in HybridPropagationModel
  6. Radar mode applies 2x absorption vs comm mode

Run:
    conda run -n pytorch_cpu python tests/_verify_absorption.py
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from utils import d_layer_absorption_dB
import config as cfg


def _make_iono():
    return {
        'iri_params':    {'dt': cfg.IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON},
        'tid_params':    {**cfg.TID, 'enable': False},
        'es_params':     {**cfg.ES,  'enable': False},
        'bubble_params': {**cfg.BUBBLE, 'enable': False},
    }


def _make_rp():
    return {'freq_MHz': cfg.FREQ_MHZ, 'Pt_W': cfg.PT_W,
            'Gt': cfg.GT, 'Gr': cfg.GR}


def run():
    p2p_fast = {**cfg.P2P, 'n_init': 6, 'max_iter': 100}
    n_fail   = 0

    # ── Check 1: A_dB > 0 for daytime ────────────────────────────────────────
    A = d_layer_absorption_dB(10.0, 30.0, chi_deg=60.0, A0=500.0)
    ok = A > 0
    print("[1] A_dB(f=10, beta=30, chi=60): {:.3f} dB  {}".format(A, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 2: A_dB = 0 at nighttime (chi = 90) ────────────────────────────
    A_night = d_layer_absorption_dB(10.0, 30.0, chi_deg=90.0, A0=500.0)
    ok = abs(A_night) < 1e-10
    print("[2] A_dB at chi=90 (night): {:.2e} dB  {}".format(A_night, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 3: A_dB decreases with higher frequency ─────────────────────────
    A_10 = d_layer_absorption_dB(10.0, 30.0)
    A_20 = d_layer_absorption_dB(20.0, 30.0)
    ok = A_20 < A_10
    print("[3] A(10MHz)={:.3f}  A(20MHz)={:.3f}  (higher f -> lower absorption)  {}".format(
        A_10, A_20, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 4: A_dB increases at lower elevation ────────────────────────────
    A_hi = d_layer_absorption_dB(10.0, 60.0)   # high elevation
    A_lo = d_layer_absorption_dB(10.0, 20.0)   # low elevation
    ok = A_lo > A_hi
    print("[4] A(beta=20)={:.3f}  A(beta=60)={:.3f}  (low elev -> more absorption)  {}".format(
        A_lo, A_hi, 'OK' if ok else 'FAIL'))
    if not ok:
        n_fail += 1

    # ── Check 5: absorption reduces Pr_W in HybridPropagationModel ───────────
    print("\n[5] Absorption reduces Pr_W:")
    from models.hybrid_model import HybridPropagationModel

    m_clean = HybridPropagationModel(
        _make_iono(), _make_rp(),
        absorption_params={'enable': False})
    modes_clean, _, _, _ = m_clean.compute(
        cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)

    m_abs = HybridPropagationModel(
        _make_iono(), _make_rp(),
        absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 60.0})
    modes_abs, _, _, _ = m_abs.compute(
        cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)

    clean_d = {m['label']: m for m in modes_clean}
    n_checked = 0
    for ma in modes_abs:
        mc = clean_d.get(ma['label'])
        if mc is None:
            continue
        ok_i = ma['Pr_W'] < mc['Pr_W']
        print("    {:16s}  clean={:.2e}W  abs={:.2e}W  {}".format(
            ma['label'], mc['Pr_W'], ma['Pr_W'], 'OK' if ok_i else 'FAIL'))
        if not ok_i:
            n_fail += 1
        n_checked += 1
    if n_checked == 0:
        print("    WARN: no matching modes to compare")

    # ── Check 6: radar mode applies 2x absorption ─────────────────────────────
    print("\n[6] Radar mode applies 2x absorption vs comm mode:")
    rp_radar = {**_make_rp(), 'sigma_rcs_m2': cfg.RADAR['sigma_rcs_m2']}

    m_comm = HybridPropagationModel(
        _make_iono(), _make_rp(),
        absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 60.0})
    modes_comm, _, _, _ = m_comm.compute(
        cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)

    m_radar = HybridPropagationModel(
        _make_iono(), rp_radar, radar_mode=True,
        absorption_params={'enable': True, 'A0': 500.0, 'chi_deg': 60.0})
    modes_radar, _, _, _ = m_radar.compute(
        cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)

    # For radar mode: base Pr is already much lower (sigma/R^4 vs Friis),
    # but the ratio of (with abs)/(without abs) should differ by 2x in dB.
    m_radar_noabs = HybridPropagationModel(
        _make_iono(), rp_radar, radar_mode=True,
        absorption_params={'enable': False})
    modes_rna, _, _, _ = m_radar_noabs.compute(
        cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)
    m_comm_noabs = HybridPropagationModel(
        _make_iono(), _make_rp(),
        absorption_params={'enable': False})
    modes_cna, _, _, _ = m_comm_noabs.compute(
        cfg.TX_POS, cfg.RX_POS, p2p_params=p2p_fast)

    radar_d = {m['label']: m for m in modes_radar}
    rna_d   = {m['label']: m for m in modes_rna}
    comm_d  = {m['label']: m for m in modes_comm}
    cna_d   = {m['label']: m for m in modes_cna}

    n_2x = 0
    n_2x_ok = 0
    for lab in sorted(radar_d.keys()):
        if lab not in (rna_d | comm_d | cna_d):
            continue
        if lab not in rna_d or lab not in comm_d or lab not in cna_d:
            continue
        dB_radar = 10*np.log10(max(radar_d[lab]['Pr_W'], 1e-30)) - \
                   10*np.log10(max(rna_d[lab]['Pr_W'],   1e-30))
        dB_comm  = 10*np.log10(max(comm_d[lab]['Pr_W'],  1e-30)) - \
                   10*np.log10(max(cna_d[lab]['Pr_W'],   1e-30))
        ratio = dB_radar / dB_comm if abs(dB_comm) > 0.01 else 0.0
        ok_i = abs(ratio - 2.0) < 0.05  # ratio should be ~2 (2x absorption)
        print("    {:16s}  radar_loss={:.2f} dB  comm_loss={:.2f} dB  ratio={:.2f}  {}".format(
            lab, -dB_radar, -dB_comm, ratio, 'OK' if ok_i else 'FAIL'))
        n_2x += 1
        if ok_i:
            n_2x_ok += 1
        else:
            n_fail += 1
    if n_2x == 0:
        print("    WARN: no modes found for 2x check")

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 50)
    if n_fail == 0:
        print("Phase 3 absorption verification PASSED.")
    else:
        print("{} check(s) FAILED.".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
