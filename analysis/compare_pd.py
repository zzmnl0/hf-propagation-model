"""
compare_pd.py - Compare measured vs. modeled P-D spectrum.

Usage:
    python analysis/compare_pd.py --measured data/pd_measured_YYYYMMDD.csv
                                  --scenario baseline
                                  --output output/compare_baseline.png

Measured CSV format (required columns: tau_ms, power_dBW):
    tau_ms, power_dBW, freq_MHz, datetime, bearing_deg
    3.90, -85.2, 10.0, 2024-01-15T04:00:00Z, 0.0
"""
import argparse
import csv
import os
import sys
import numpy as np


# ── Data loading ─────────────────────────────────────────────────────────────

def load_measured_pd(csv_path: str) -> tuple:
    """
    Load measured P-D spectrum from CSV.

    Required columns: tau_ms, power_dBW
    Optional columns: freq_MHz, datetime, bearing_deg

    Returns
    -------
    tau_ms  : (N,) array  group delays [ms]
    pd_dBW  : (N,) array  measured power [dBW]
    meta    : dict  optional metadata from first data row
    """
    tau_ms = []
    pd_dBW = []
    meta   = {}

    with open(csv_path, newline='', encoding='utf-8') as f:
        lines = [l for l in f if l.strip() and not l.strip().startswith('#')]

    if len(lines) < 2:
        return np.array([]), np.array([]), {}

    reader = csv.DictReader(lines)
    for i, row in enumerate(reader):
        try:
            tau_ms.append(float(row['tau_ms']))
            pd_dBW.append(float(row['power_dBW']))
        except (KeyError, ValueError):
            continue
        if i == 0:
            for key in ('freq_MHz', 'datetime', 'bearing_deg'):
                if key in row and row[key].strip():
                    try:
                        meta[key] = float(row[key]) if key != 'datetime' else row[key]
                    except ValueError:
                        meta[key] = row[key]

    return np.array(tau_ms), np.array(pd_dBW), meta


def load_model_modes(scenario: str, out_dir: str) -> list:
    """
    Load mode results from a scenario CSV saved by save_modes_csv.

    Returns list of dicts with at least 'label', 'tau_ms', 'Pr_dBW'.
    """
    csv_path = os.path.join(out_dir, 'modes_{}.csv'.format(scenario))
    if not os.path.exists(csv_path):
        raise FileNotFoundError('Model CSV not found: {}'.format(csv_path))

    modes = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                m = {
                    'label':  row.get('label', ''),
                    'tau_ms': float(row.get('tau_ms', 0)),
                    'Pr_dBW': float(row.get('Pr_dBW', -999)),
                    'Pr_W':   float(row.get('Pr_W', 0)),
                }
                modes.append(m)
            except (KeyError, ValueError):
                continue
    return modes


# ── Visualization ─────────────────────────────────────────────────────────────

def overlay_plot(tau_model,  pd_model_W,  modes_model,
                 tau_meas=None, pd_meas_dBW=None,
                 title='Model vs Measured', save_path=None):
    """
    Two-panel comparison plot.

    Top    : P-D spectrum overlay (model blue, measured orange)
    Bottom : mode markers at predicted delays

    Returns matplotlib Figure.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True)

    if len(pd_model_W) > 0 and len(tau_model) > 0:
        pd_dBW_model = 10.0 * np.log10(np.maximum(pd_model_W, 1e-30))
        ax1.plot(tau_model, pd_dBW_model, 'b-', lw=1.5, label='Model')

    if (tau_meas is not None and pd_meas_dBW is not None
            and len(tau_meas) > 0):
        ax1.plot(tau_meas, pd_meas_dBW, 'o-', color='orange',
                 lw=1.2, ms=4, label='Measured')

    ax1.set_ylabel('Power [dBW]')
    ax1.set_title(title)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    for m in modes_model:
        ax2.axvline(m['tau_ms'], color='b', alpha=0.55, lw=1.0)
        ax2.text(m['tau_ms'], 0.5, m['label'],
                 rotation=90, fontsize=6, ha='right', va='center',
                 transform=ax2.get_xaxis_transform())

    ax2.set_xlabel('Group delay [ms]')
    ax2.set_ylabel('Modes')
    ax2.set_yticks([])
    ax2.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print('  Saved: {}'.format(save_path))

    return fig


# ── Scoring ──────────────────────────────────────────────────────────────────

def score_mode_match(modes_model: list, peaks_meas: list,
                     tau_tol_ms: float = 0.5) -> dict:
    """
    Compute hit rate for main-mode identification.

    Parameters
    ----------
    modes_model : list of mode dicts (need 'tau_ms'; optionally 'Pr_dBW')
    peaks_meas  : list of dicts with 'tau_ms' (optionally 'power_dBW')
    tau_tol_ms  : matching window [ms]

    Returns
    -------
    dict: n_hit, n_total, hit_rate, tau_error_ms (list), Pr_error_dB (list)
    """
    n_total   = len(peaks_meas)
    n_hit     = 0
    tau_errors = []
    Pr_errors  = []

    for pk in peaks_meas:
        tau_pk = float(pk['tau_ms'])
        if not modes_model:
            continue
        closest = min(modes_model, key=lambda m: abs(m['tau_ms'] - tau_pk))
        dtau = abs(closest['tau_ms'] - tau_pk)
        if dtau <= tau_tol_ms:
            n_hit += 1
            tau_errors.append(closest['tau_ms'] - tau_pk)
            if 'power_dBW' in pk and 'Pr_dBW' in closest:
                Pr_errors.append(closest['Pr_dBW'] - float(pk['power_dBW']))

    return {
        'n_hit':        n_hit,
        'n_total':      n_total,
        'hit_rate':     n_hit / max(n_total, 1),
        'tau_error_ms': tau_errors,
        'Pr_error_dB':  Pr_errors,
    }


# ── CLI entry point ───────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare measured vs. modeled P-D spectrum')
    parser.add_argument('--measured',  required=True,
                        help='Measured P-D CSV file (tau_ms, power_dBW)')
    parser.add_argument('--scenario',  default='baseline',
                        help='Model scenario label (default: baseline)')
    parser.add_argument('--output',    default=None,
                        help='Output plot path (default: output/compare_<scenario>.png)')
    parser.add_argument('--tau_tol',   type=float, default=0.5,
                        help='Mode match tolerance [ms] (default: 0.5)')
    args = parser.parse_args()

    # Locate project root (one level above analysis/)
    _root   = os.path.join(os.path.dirname(__file__), '..')
    _out    = os.path.join(_root, 'output')
    sys.path.insert(0, _root)

    # Load measured data
    tau_m, pd_m, meta = load_measured_pd(args.measured)
    print('Loaded {} measured points  meta={}'.format(len(tau_m), meta))

    # Load model modes
    try:
        modes = load_model_modes(args.scenario, _out)
        print('Loaded {} model modes from output/modes_{}.csv'.format(
            len(modes), args.scenario))
    except FileNotFoundError as e:
        print(str(e))
        sys.exit(1)

    # Score
    peaks_meas = [{'tau_ms': t, 'power_dBW': p}
                  for t, p in zip(tau_m, pd_m)]
    score = score_mode_match(modes, peaks_meas, tau_tol_ms=args.tau_tol)
    print('Hit rate: {}/{} = {:.1%}  tol={} ms'.format(
        score['n_hit'], score['n_total'],
        score['hit_rate'], args.tau_tol))
    if score['tau_error_ms']:
        print('Mean tau error: {:.3f} ms  rms: {:.3f} ms'.format(
            float(np.mean(score['tau_error_ms'])),
            float(np.std(score['tau_error_ms']))))
    if score['Pr_error_dB']:
        print('Mean Pr error: {:.1f} dB'.format(
            float(np.mean(score['Pr_error_dB']))))

    # Plot
    save_path = args.output or os.path.join(
        _out, 'compare_{}.png'.format(args.scenario))
    tau_model = np.array([m['tau_ms'] for m in modes])
    pd_model_W = np.array([m['Pr_W'] for m in modes])

    overlay_plot(
        tau_model, pd_model_W, modes,
        tau_meas=tau_m, pd_meas_dBW=pd_m,
        title='Model ({}) vs Measured'.format(args.scenario),
        save_path=save_path)
