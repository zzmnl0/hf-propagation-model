"""Power-delay spectrum construction and main-mode identification."""

from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks

from config import UPDATE2


def build_pd_spectrum(mode_results: list,
                      tau_axis: np.ndarray | None = None,
                      tau_res_ms: float = 0.02,
                      method: str = "gaussian") -> tuple[np.ndarray, np.ndarray]:
    """Construct a Gaussian, tube-reference, or scatter-coupled P-D spectrum."""
    valid = [
        mode for mode in mode_results
        if mode.get("convergence_status", "passed") == "passed"
        and not mode.get("under_resolved_Ne", False)
        and mode.get("branch_status") not in ("branch_lost", "constraint_failed")
    ]
    if not valid:
        return np.array([0.0]), np.array([0.0])
    if method == "gaussian" and UPDATE2.get("SCATTER_PD", {}).get("apply_to_pd", False):
        method = "scatter_coupled"
    if method == "scatter_coupled":
        from .pd_scatter_kernel import apply_scatter_kernel
        return apply_scatter_kernel(valid, tau_axis, tau_res_ms)

    taus = [mode["tau_ms"] for mode in valid]
    if tau_axis is None:
        tau_axis = np.arange(max(0.0, min(taus) - 2.0), max(taus) + 2.0,
                             tau_res_ms)
    if method == "tube":
        pd = np.zeros(len(tau_axis))
        for mode in valid:
            tau_c = float(mode["tau_ms"])
            power = float(mode.get("Pr_W", 0.0))
            spread = max(float(mode.get("tau_spread_ms", 0.0)),
                         float(mode.get("delta_tau_ms", 0.0)), tau_res_ms)
            mask = ((tau_axis >= tau_c - spread / 2.0)
                    & (tau_axis <= tau_c + spread / 2.0))
            if np.any(mask):
                pd[mask] += power / int(mask.sum())
        return tau_axis, pd
    if method != "gaussian":
        raise ValueError(f"unsupported P-D method: {method}")

    pd = np.zeros(len(tau_axis))
    for mode in valid:
        tau_c = float(mode["tau_ms"])
        sigma = max(float(mode.get("delta_tau_ms", tau_res_ms)), tau_res_ms)
        pd += float(mode.get("Pr_W", 0.0)) * np.exp(
            -(tau_axis - tau_c) ** 2 / (2.0 * sigma ** 2))
    return tau_axis, pd


def _find_ox_pairs(mode_results: list) -> list:
    pairs = []
    ordinary = {mode["label"].replace("_O", ""): mode
                for mode in mode_results if mode.get("wave_mode") == "O"}
    extraordinary = {mode["label"].replace("_X", ""): mode
                     for mode in mode_results if mode.get("wave_mode") == "X"}
    for base, mode_o in ordinary.items():
        mode_x = extraordinary.get(base)
        if mode_x is None:
            continue
        delay = abs(mode_o["tau_ms"] - mode_x["tau_ms"])
        if delay < 0.5:
            pairs.append({
                "label_O": mode_o["label"], "label_X": mode_x["label"],
                "tau_O": mode_o["tau_ms"], "tau_X": mode_x["tau_ms"],
                "delta_tau_OX_ms": delay,
                "Pr_O_dBW": mode_o.get("Pr_dBW", -999.0),
                "Pr_X_dBW": mode_x.get("Pr_dBW", -999.0),
                "delta_Pr_OX_dB": (mode_o.get("Pr_dBW", -999.0)
                                    - mode_x.get("Pr_dBW", -999.0)),
            })
    return pairs


def identify_main_mode(tau_axis: np.ndarray, pd_W: np.ndarray,
                       mode_results: list) -> tuple[dict | None, list, dict]:
    """Rank valid propagation modes by peaks on the supplied P-D spectrum."""
    empty = {
        "n_modes": 0, "main_label": "", "main_tau_ms": 0.0,
        "main_Pr_dBW": -999.0, "main_F_focus": 0.0,
        "tau_spread_main": 0.0, "mode_pairs": [],
    }
    valid = [
        mode for mode in mode_results
        if mode.get("convergence_status", "passed") == "passed"
        and not mode.get("under_resolved_Ne", False)
        and mode.get("branch_status") not in ("branch_lost", "constraint_failed")
    ]
    if len(pd_W) == 0 or not valid:
        return None, [], empty
    threshold = max(float(np.max(pd_W)) * 0.01, 1.0e-30)
    peaks, _ = find_peaks(pd_W, height=threshold, distance=5)
    if len(peaks) == 0:
        peaks = np.array([int(np.argmax(pd_W))])
    ordered = peaks[np.argsort(-pd_W[peaks])]
    ranked = [min(valid, key=lambda mode: abs(mode["tau_ms"] - tau_axis[index]))
              for index in ordered]
    main = ranked[0]
    return main, ranked, {
        "n_modes": len(valid), "main_label": main.get("label", ""),
        "main_tau_ms": main["tau_ms"],
        "main_Pr_dBW": main.get("Pr_dBW", -999.0),
        "main_F_focus": main.get("F_focus", 1.0),
        "tau_spread_main": main.get("tau_spread_ms", 0.0),
        "mode_pairs": _find_ox_pairs(valid),
    }
