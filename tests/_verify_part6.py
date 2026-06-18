"""
Verify Part 6: Hybrid propagation model (M5-M7 synthesis).

Checks:
  [1] build_pd_spectrum: Gaussian spreading, tau axis coverage, peak positions
  [2] identify_main_mode: highest P-D peak maps to highest-power mode
  [3] HybridPropagationModel.__init__: constructs without error (IRI-only)
  [4] HybridPropagationModel.compute: pipeline runs, finds modes, returns valid outputs
  [5] mode_results dict: all required keys present
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.signal import find_peaks
import config as cfg
from models.hybrid_model import (build_pd_spectrum, identify_main_mode,
                                  HybridPropagationModel)

print("=" * 55)
print("  Part 6 verification")
print("=" * 55)

# ── [1] build_pd_spectrum ─────────────────────────────────────────────────────
print("\n[1] build_pd_spectrum: Gaussian spreading and peak positions")

modes_syn = [
    {'tau_ms': 5.0, 'Pr_W': 1.0, 'delta_tau_ms': 0.10},
    {'tau_ms': 7.0, 'Pr_W': 2.0, 'delta_tau_ms': 0.05},
]
tau_axis, pd_W = build_pd_spectrum(modes_syn, tau_res_ms=0.01)

assert tau_axis[0]  <= 3.0, "tau_axis starts too late: {:.3f}".format(tau_axis[0])
assert tau_axis[-1] >= 8.9, "tau_axis ends too early: {:.3f}".format(tau_axis[-1])
print("  tau_axis: [{:.2f}, {:.2f}] ms  ({} pts)".format(
    tau_axis[0], tau_axis[-1], len(tau_axis)))

for m in modes_syn:
    idx     = int(np.argmin(np.abs(tau_axis - m['tau_ms'])))
    val     = float(pd_W[idx])
    assert val > 0.0, "P-D zero at tau={:.1f}ms".format(m['tau_ms'])
    assert val > 0.1 * m['Pr_W'], \
        "P-D too low at tau={:.1f}ms: {:.3e}".format(m['tau_ms'], val)
    print("  P-D at tau={:.1f}ms: {:.4e} W".format(m['tau_ms'], val))

assert np.all(pd_W >= 0.0), "Negative P-D values"

peaks_idx, _ = find_peaks(pd_W, height=pd_W.max() * 0.01, distance=50)
assert len(peaks_idx) >= 2, "Expected >= 2 peaks, found {}".format(len(peaks_idx))
print("  Peaks detected: {}  (expect >= 2)".format(len(peaks_idx)))
print("  [PASS]")

# ── [2] identify_main_mode ────────────────────────────────────────────────────
print("\n[2] identify_main_mode: highest P-D peak maps to highest-power mode")

main_mode, ranked, _ = identify_main_mode(tau_axis, pd_W, modes_syn)

assert main_mode is not None,  "main_mode is None"
assert len(ranked) >= 1,       "ranked list empty"
# strongest mode is tau=7.0ms (Pr_W=2.0)
assert abs(main_mode['tau_ms'] - 7.0) < 0.5, \
    "main_mode tau {:.3f} not near 7.0ms".format(main_mode['tau_ms'])
print("  main_mode: tau={:.3f}ms  Pr_W={:.3f}W  (expect tau~7.0, Pr~2.0)".format(
    main_mode['tau_ms'], main_mode.get('Pr_W', float('nan'))))
print("  ranked modes: {}".format(len(ranked)))
print("  [PASS]")

# ── [3] HybridPropagationModel.__init__ (IRI only) ───────────────────────────
print("\n[3] HybridPropagationModel.__init__: IRI-only construction")

iono_params = {
    'iri_params': {
        'dt':  cfg.IRI_DT,
        'lat': cfg.IRI_LAT,
        'lon': cfg.IRI_LON,
    },
}
radar_params = {
    'Pt_W':     cfg.PT_W,
    'Gt':       cfg.GT,
    'Gr':       cfg.GR,
    'freq_MHz': cfg.FREQ_MHZ,
}

model = HybridPropagationModel(iono_params, radar_params)
assert model.freq == cfg.FREQ_MHZ, "freq mismatch"
assert model.iono is not None,     "iono model not constructed"
assert model.es   is None,         "es should be None (not enabled)"
assert model.pe   is not None,     "pe model not constructed"
print("  freq_MHz : {:.1f}".format(model.freq))
print("  iono     : {}".format(type(model.iono).__name__))
print("  es       : {}".format(model.es))
print("  pe       : {}".format(type(model.pe).__name__))
print("  [PASS]")

# ── [4] HybridPropagationModel.compute (IRI only, small n_init) ──────────────
print("\n[4] HybridPropagationModel.compute: full pipeline (IRI, n_init=6)")

p2p_fast = dict(cfg.P2P)
p2p_fast['n_init']   = 6
p2p_fast['max_iter'] = 200

mode_results, tau_axis_c, pd_W_c, main_mode_c = model.compute(
    tx_km      = cfg.TX_POS,
    rx_km      = cfg.RX_POS,
    t          = 0.0,
    p2p_params = p2p_fast,
)

assert len(mode_results) > 0,     "No propagation modes found"
assert len(tau_axis_c)  > 1,      "tau_axis is trivial"
assert np.all(pd_W_c >= 0.0),     "Negative P-D values in output"
assert main_mode_c is not None,   "main_mode is None after compute"

tau_list = [m['tau_ms'] for m in mode_results]
print("  Modes found : {}".format(len(mode_results)))
print("  tau range   : [{:.2f}, {:.2f}] ms".format(min(tau_list), max(tau_list)))
print("  pd_W max    : {:.4e} W".format(float(pd_W_c.max())))
print("  main_mode   : label='{}' tau={:.3f}ms Pr={:.3e}W".format(
    main_mode_c['label'], main_mode_c['tau_ms'], main_mode_c['Pr_W']))
print("  [PASS]")

# ── [5] Required output keys ──────────────────────────────────────────────────
print("\n[5] mode_results dict: all required keys present")

REQUIRED = {'label', 'tau_ms', 'delta_tau_ms', 'Pr_W', 'Pr_dBW',
            'h_reflect_km', 'group_path_km', 'beta_deg', 'phi_deg'}

for i, m in enumerate(mode_results):
    missing = REQUIRED - set(m.keys())
    assert not missing, "Mode {} missing keys: {}".format(i, missing)
    assert np.isfinite(m['tau_ms']),        "tau_ms not finite"
    assert m['tau_ms'] > 0.0,              "tau_ms <= 0"
    assert np.isfinite(m['Pr_W']),         "Pr_W not finite"
    assert m['Pr_W'] >= 0.0,              "Pr_W negative"
    assert np.isfinite(m['h_reflect_km']), "h_reflect_km not finite"
    assert m['h_reflect_km'] > 0.0,       "h_reflect_km <= 0"

print("  All {} modes have required keys.".format(len(mode_results)))
for m in mode_results:
    print("    {} | tau={:.2f}ms | h={:.0f}km | Pr={:.1f}dBW".format(
        m['label'], m['tau_ms'], m['h_reflect_km'], m['Pr_dBW']))
print("  [PASS]")

print()
print("=" * 55)
print("  All Part 6 checks PASSED.")
print("=" * 55)
