"""
analysis/mode_analyzer.py
Post-processing analysis tools for HybridPropagationModel mode results.

Functions:
  analyze_ox_pairs(mode_results)              -> list[dict]
  analyze_mode_features(mode_results)         -> dict
  build_mode_summary(mode_results, main_mode) -> dict
  print_mode_report(mode_results, summary)

CLI:
  python analysis/mode_analyzer.py --csv output/modes_<scenario>.csv
"""
import argparse
import os
import sys
import numpy as np


# ── Public API ────────────────────────────────────────────────────────────────

def analyze_ox_pairs(mode_results):
    """
    Match O/X wave mode pairs by base label and delay proximity.

    Returns list of dicts with keys:
      label_O, label_X, tau_O, tau_X, delta_tau_OX_ms,
      Pr_O_dBW, Pr_X_dBW, delta_Pr_OX_dB
    Returns [] for tube_mode results (wave_mode='iso').
    """
    o_modes = [m for m in mode_results if m.get('wave_mode') == 'O']
    x_modes = [m for m in mode_results if m.get('wave_mode') == 'X']
    if not o_modes or not x_modes:
        return []

    def _base(label):
        return label.replace('_O', '').replace('_X', '')

    pairs = []
    for om in o_modes:
        base_o = _base(om['label'])
        best_x = None
        best_dt = 0.5
        for xm in x_modes:
            if _base(xm['label']) != base_o:
                continue
            dt = abs(xm['tau_ms'] - om['tau_ms'])
            if dt < best_dt:
                best_dt = dt
                best_x  = xm
        if best_x is not None:
            pairs.append({
                'label_O'        : om['label'],
                'label_X'        : best_x['label'],
                'tau_O'          : om['tau_ms'],
                'tau_X'          : best_x['tau_ms'],
                'delta_tau_OX_ms': best_x['tau_ms'] - om['tau_ms'],
                'Pr_O_dBW'       : om.get('Pr_dBW', -999.0),
                'Pr_X_dBW'       : best_x.get('Pr_dBW', -999.0),
                'delta_Pr_OX_dB' : om.get('Pr_dBW', -999.0) - best_x.get('Pr_dBW', -999.0),
            })
    return pairs


def analyze_mode_features(mode_results):
    """
    Compute inter-mode statistics for delay and power.

    Returns dict:
      labels, tau_ms, Pr_dBW, delta_tau_ij, delta_P_ij_dB,
      F_focus, tau_spread_ms
    Modes are sorted by Pr_W descending (highest power first).
    """
    if not mode_results:
        return {
            'labels': [], 'tau_ms': [], 'Pr_dBW': [],
            'delta_tau_ij': np.zeros((0, 0)),
            'delta_P_ij_dB': np.zeros((0, 0)),
            'F_focus': [], 'tau_spread_ms': [],
        }

    sorted_m = sorted(mode_results,
                      key=lambda m: m.get('Pr_W', 0.0), reverse=True)
    labels   = [m.get('label', '') for m in sorted_m]
    taus     = [m['tau_ms']             for m in sorted_m]
    prs      = [m.get('Pr_dBW', -999.0) for m in sorted_m]
    ffocus   = [m.get('F_focus', 1.0)   for m in sorted_m]
    tspread  = [m.get('tau_spread_ms', 0.0) for m in sorted_m]

    n = len(sorted_m)
    dt_ij = np.zeros((n, n))
    dp_ij = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dt_ij[i, j] = taus[i] - taus[j]
            dp_ij[i, j] = prs[i]  - prs[j]

    return {
        'labels'       : labels,
        'tau_ms'       : taus,
        'Pr_dBW'       : prs,
        'delta_tau_ij' : dt_ij,
        'delta_P_ij_dB': dp_ij,
        'F_focus'      : ffocus,
        'tau_spread_ms': tspread,
    }


def build_mode_summary(mode_results, main_mode=None):
    """
    Build mode_summary dict independently of HybridPropagationModel.

    Returns:
      n_modes, main_label, main_tau_ms, main_Pr_dBW,
      main_F_focus, tau_spread_main, mode_pairs
    """
    _empty = {
        'n_modes': 0, 'main_label': '', 'main_tau_ms': 0.0,
        'main_Pr_dBW': -999.0, 'main_F_focus': 0.0,
        'tau_spread_main': 0.0, 'mode_pairs': [],
    }
    if not mode_results:
        return _empty

    if main_mode is None:
        main_mode = max(mode_results, key=lambda m: m.get('Pr_W', 0.0))

    return {
        'n_modes'        : len(mode_results),
        'main_label'     : main_mode.get('label', ''),
        'main_tau_ms'    : main_mode.get('tau_ms', 0.0),
        'main_Pr_dBW'    : main_mode.get('Pr_dBW', -999.0),
        'main_F_focus'   : main_mode.get('F_focus', 1.0),
        'tau_spread_main': main_mode.get('tau_spread_ms', 0.0),
        'mode_pairs'     : analyze_ox_pairs(mode_results),
    }


def print_mode_report(mode_results, summary=None):
    """
    Print formatted mode table with F_focus and tau_spread columns,
    followed by main-mode summary and O/X pair info.
    """
    if not mode_results:
        print("  (no modes)")
        return

    if summary is None:
        summary = build_mode_summary(mode_results)

    sorted_m = sorted(mode_results,
                      key=lambda m: m.get('Pr_W', 0.0), reverse=True)

    print("")
    print("  {:<12} {:>8} {:>10} {:>8} {:>12}".format(
        "label", "tau(ms)", "Pr(dBW)", "F_focus", "tauSprd(ms)"))
    print("  " + "-" * 56)
    for m in sorted_m:
        print("  {:<12} {:>8.3f} {:>10.1f} {:>8.3f} {:>12.4f}".format(
            m.get('label', '?'),
            m.get('tau_ms', 0.0),
            m.get('Pr_dBW', -999.0),
            m.get('F_focus', 1.0),
            m.get('tau_spread_ms', 0.0)))
    print("")

    print("  Summary: n_modes={}  main={}  tau={:.3f}ms  Pr={:.1f}dBW  F_focus={:.2f}".format(
        summary['n_modes'], summary['main_label'],
        summary['main_tau_ms'], summary['main_Pr_dBW'],
        summary['main_F_focus']))

    pairs = summary.get('mode_pairs', [])
    if pairs:
        print("  O/X pairs:")
        for p in pairs:
            print("    {}/{}: dTau={:.3f}ms  dPr={:.1f}dB".format(
                p['label_O'], p['label_X'],
                p['delta_tau_OX_ms'], p['delta_Pr_OX_dB']))


# ── CLI entry ─────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Analyze propagation mode features from saved CSV')
    parser.add_argument('--csv', required=True,
        help='Mode CSV file (output/modes_<scenario>.csv)')
    parser.add_argument('--scenario', default=None,
        help='Scenario label (overrides CSV filename)')
    args = parser.parse_args()

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from analysis.compare_pd import load_model_modes

    scenario = (args.scenario
                or os.path.basename(args.csv)
                   .replace('modes_', '').replace('.csv', ''))
    out_dir  = os.path.dirname(args.csv)
    modes    = load_model_modes(scenario, out_dir)

    features = analyze_mode_features(modes)
    ox_pairs = analyze_ox_pairs(modes)
    summary  = build_mode_summary(modes)

    print_mode_report(modes, summary)

    if ox_pairs:
        print("\nO/X pairs:")
        for p in ox_pairs:
            print("  {}/{}: dTau={:.3f}ms  dPr={:.1f}dB".format(
                p['label_O'], p['label_X'],
                p['delta_tau_OX_ms'], p['delta_Pr_OX_dB']))
