"""
Phase 5 verification: 3-D ionospheric background upgrades.

Checks:
  1. B1 lateral IRI: NmF2 varies along the link (TX lat != RX lat)
  2. B2 IGRF 3-D:    fH_MHz and dip_deg vary with horizontal position
  3. B3 multi-TID:   2-component superposition produces two spatial periods
  4. Backward compat: lateral off + n_components=1 gives same Ne as legacy path

Run:
    conda run -n pytorch_cpu python tests/_verify_phase5.py
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import config as cfg
from config import BG_X, BG_Z, TX_LAT, TX_LON, LINK_BEARING_DEG, IRI_DT
from models.ionosphere_model import IonosphereModel
from models.geomag_field import GeomagnField2D
from utils import destination_point


def run():
    n_fail = 0

    print("=" * 55)
    print("  Phase 5 verification: 3-D ionospheric background")
    print("=" * 55)

    # ── Check 1: B1 lateral IRI -- NmF2 varies along link ─────────────────────
    print("\nCheck 1: B1 lateral IRI - NmF2 varies TX->RX")

    lat_rx, lon_rx = destination_point(TX_LAT, TX_LON, LINK_BEARING_DEG,
                                       float(cfg.RX_RANGE))
    lat_mid = (TX_LAT + lat_rx) / 2.0

    iri_params_lat = {
        'dt':  IRI_DT,
        'lat': lat_mid,
        'lon': TX_LON,
    }
    lat_iri = {
        'enable':      True,
        'spacing_km':  50.0,
        'tx_lat':      TX_LAT,
        'tx_lon':      TX_LON,
        'bearing_deg': LINK_BEARING_DEG,
    }

    print("  Building lateral IRI (28 IRI calls, ~10s) ...")
    model_lat = IonosphereModel(iri_params=iri_params_lat,
                                lateral_iri_params=lat_iri)
    Ne_lat, _ = model_lat.build_Ne_field(BG_X, BG_Z)

    # NmF2 at TX end (x~0) vs RX end (x~1169 km)
    ix_tx = np.argmin(np.abs(BG_X - 0.0))
    ix_rx = np.argmin(np.abs(BG_X - cfg.RX_RANGE))
    NmF2_tx = Ne_lat[ix_tx, :].max()
    NmF2_rx = Ne_lat[ix_rx, :].max()
    diff_pct = abs(NmF2_tx - NmF2_rx) / max(NmF2_tx, 1.0) * 100.0
    ok1 = diff_pct > 1.0   # expect at least 1% difference over 1169 km at different latitudes
    print("  NmF2 TX={:.3e} m^-3  RX={:.3e} m^-3  diff={:.1f}%  {}".format(
        NmF2_tx, NmF2_rx, diff_pct, 'OK' if ok1 else 'FAIL'))
    if not ok1:
        n_fail += 1

    # ── Check 2: B2 IGRF 3-D -- fH and dip vary horizontally ──────────────────
    print("\nCheck 2: B2 IGRF 3-D - fH and dip vary along link")

    # Use a coarse x grid for speed (full BG_X would be 271 x 271 ppigrf queries)
    x_coarse = np.arange(0.0, cfg.RX_RANGE + 100.0, 100.0)
    gf = GeomagnField2D(x_coarse, BG_Z, TX_LAT, TX_LON, LINK_BEARING_DEG, IRI_DT)

    pts_tx = np.array([[0.0,    300.0]])
    pts_rx = np.array([[cfg.RX_RANGE, 300.0]])
    fH_tx  = float(gf.fH_MHz_batch(pts_tx)[0])
    fH_rx  = float(gf.fH_MHz_batch(pts_rx)[0])
    dip_tx = float(gf.dip_deg_batch(pts_tx)[0])
    dip_rx = float(gf.dip_deg_batch(pts_rx)[0])

    dfH  = abs(fH_tx  - fH_rx)
    ddip = abs(dip_tx - dip_rx)
    ok2a = dfH  > 0.005   # > 5 kHz difference expected over 1169 km N-S path
    ok2b = ddip > 0.5     # > 0.5 deg dip difference expected
    print("  fH  TX={:.3f} MHz  RX={:.3f} MHz  diff={:.3f} MHz  {}".format(
        fH_tx, fH_rx, dfH, 'OK' if ok2a else 'FAIL'))
    print("  dip TX={:.1f} deg  RX={:.1f} deg  diff={:.1f} deg  {}".format(
        dip_tx, dip_rx, ddip, 'OK' if ok2b else 'FAIL'))
    if not ok2a:
        n_fail += 1
    if not ok2b:
        n_fail += 1

    # ── Check 3: B3 multi-TID -- 2 components create 2 spatial periods ────────
    print("\nCheck 3: B3 multi-TID - 2-component superposition")

    # Single-component reference
    tid_single = {**cfg.TID,
                  'enable':      True,
                  'n_components': 1,
                  'amplitude':   0.10,
                  'lambda_h_km': 300.0,
                  'T_s':         2400.0}
    m1 = IonosphereModel(tid_params=tid_single)
    Ne1, _ = m1.build_Ne_field(BG_X, BG_Z)

    # Two-component: 300 km (0 deg) + 200 km (30 deg)
    tid_multi = {**cfg.TID,
                 'enable':           True,
                 'n_components':     2,
                 'az_deg_list':      [0.0,   30.0],
                 'amplitude_list':   [0.08,  0.06],
                 'period_s_list':    [2400.0, 1800.0],
                 'lambda_h_km_list': [300.0,  200.0],
                 'link_bearing_deg': LINK_BEARING_DEG}
    m2 = IonosphereModel(tid_params=tid_multi)
    Ne2, _ = m2.build_Ne_field(BG_X, BG_Z)

    # Ne fields should differ (multi-component adds a second spatial frequency)
    # Check at F2-peak altitude
    iz_f2 = np.argmax(Ne1[len(BG_X)//2, :])
    ne1_slice = Ne1[:, iz_f2]
    ne2_slice = Ne2[:, iz_f2]
    rms_diff = np.sqrt(np.mean((ne2_slice - ne1_slice)**2))
    rel_diff  = rms_diff / max(ne1_slice.mean(), 1.0)
    ok3 = rel_diff > 0.001   # multi-component should differ by >0.1% rms
    print("  1-comp vs 2-comp Ne at F2 peak (z={:.0f}km): rms_diff={:.2f}% {}".format(
        BG_Z[iz_f2], rel_diff * 100.0, 'OK' if ok3 else 'FAIL'))
    if not ok3:
        n_fail += 1

    # ── Check 4: backward compatibility ──────────────────────────────────────
    print("\nCheck 4: backward compat - lateral off + n_comp=1 = legacy result")

    # Legacy model (original behavior)
    m_legacy = IonosphereModel()
    Ne_leg, _ = m_legacy.build_Ne_field(BG_X, BG_Z)

    # New model with all Phase 5 features disabled
    iri_compat = {'dt': IRI_DT, 'lat': cfg.IRI_LAT, 'lon': cfg.IRI_LON}
    lat_off    = {'enable': False}
    m_new      = IonosphereModel(iri_params=iri_compat, lateral_iri_params=lat_off)
    Ne_new, _  = m_new.build_Ne_field(BG_X, BG_Z)

    max_absdiff = np.max(np.abs(Ne_new - Ne_leg))
    ok4 = max_absdiff < 1.0   # should be identical (floating-point exact)
    print("  max |Ne_new - Ne_leg| = {:.2e} m^-3  {}".format(
        max_absdiff, 'OK' if ok4 else 'FAIL'))
    if not ok4:
        n_fail += 1

    print("\n{} fail(s)\n".format(n_fail))
    return n_fail


if __name__ == '__main__':
    sys.exit(run())
