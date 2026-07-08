"""Generate paper-oriented result summaries and figures.

Run from code/:
    python analysis/report_results.py
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

SCENARIO_ORDER = [
    "baseline",
    "tid",
    "es",
    "bubble",
    "full",
    "OX",
    "radar",
    "spreadf",
    "tube_radar",
]

SCENARIO_LABELS = {
    "baseline": "Baseline",
    "tid": "TID",
    "es": "Es",
    "bubble": "Bubble",
    "full": "Full",
    "OX": "O/X",
    "radar": "Radar",
    "spreadf": "Spread-F",
    "tube_radar": "Tube radar",
}

NUMERIC_COLUMNS = [
    "freq_MHz",
    "tau_ms",
    "tau_2way_ms",
    "delta_tau_ms",
    "Pr_W",
    "Pr_dBW",
    "pd_power_W",
    "h_reflect_km",
    "group_path_km",
    "beta_deg",
    "phi_deg",
    "geom_path_km",
    "x_land_km",
    "y_land_km",
    "z_land_km",
    "az_deg_tx",
    "az_deg_rx",
    "azimuth_deflect_deg",
    "az_deflect_deg",
    "abs_dB",
    "abs_total_dB",
    "OX_gp_ratio",
    "theta_Es_deg",
    "aspect_Es_deg",
    "es_loss_dB",
    "sigma_phi_sq",
    "sigma_phi_rad",
    "S4",
    "sigma_tau_ms",
    "F_focus",
    "F_focus_el",
    "F_focus_az",
    "F_focus_3d",
    "A_tube_km2",
    "A_tube_3d_km2",
    "tau_spread_ms",
]

DETAIL_COLUMNS = [
    "scenario",
    "label",
    "wave_mode",
    "freq_MHz",
    "tau_ms",
    "tau_2way_ms",
    "group_path_km",
    "geom_path_km",
    "h_reflect_km",
    "Pr_dBW",
    "beta_deg",
    "y_land_km",
    "az_deflect_deg",
    "abs_total_dB",
    "scatter_method",
    "sigma_phi_sq",
    "S4",
    "sigma_tau_ms",
    "F_focus_3d",
    "A_tube_km2",
    "tau_spread_ms",
]


def _read_modes() -> pd.DataFrame:
    frames = []
    for path in sorted(OUT.glob("modes_*.csv")):
        scenario = path.stem.replace("modes_", "")
        df = pd.read_csv(path)
        if "scenario" not in df.columns:
            df["scenario"] = scenario
        df["scenario"] = scenario
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    data = pd.concat(frames, ignore_index=True, sort=False)
    for col in NUMERIC_COLUMNS:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors="coerce")
    for col in DETAIL_COLUMNS:
        if col not in data.columns:
            data[col] = np.nan
    data["scenario_label"] = data["scenario"].map(SCENARIO_LABELS).fillna(data["scenario"])
    return data


def _main_modes(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scenario in SCENARIO_ORDER:
        sdf = data[data["scenario"] == scenario]
        if sdf.empty:
            continue
        if "Pr_W" in sdf and sdf["Pr_W"].notna().any():
            idx = sdf["Pr_W"].idxmax()
        else:
            idx = sdf["Pr_dBW"].idxmax()
        row = sdf.loc[idx].copy()
        row["n_modes"] = int(len(sdf))
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    summary = pd.DataFrame(rows)
    base = summary[summary["scenario"] == "baseline"]
    if not base.empty:
        b = base.iloc[0]
        summary["delta_tau_us"] = (summary["tau_ms"] - b.get("tau_ms", np.nan)) * 1000.0
        summary["delta_group_path_km"] = summary["group_path_km"] - b.get("group_path_km", np.nan)
        summary["delta_Pr_dB"] = summary["Pr_dBW"] - b.get("Pr_dBW", np.nan)
    else:
        summary["delta_tau_us"] = np.nan
        summary["delta_group_path_km"] = np.nan
        summary["delta_Pr_dB"] = np.nan
    summary["scenario_sort"] = summary["scenario"].map(
        {name: i for i, name in enumerate(SCENARIO_ORDER)}
    )
    return summary.sort_values("scenario_sort")


def _add_ox_metrics(data: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    summary = summary.copy()
    summary["delta_tau_OX_us"] = np.nan
    summary["delta_group_path_OX_km"] = np.nan
    ox = data[data["scenario"] == "OX"].copy()
    if ox.empty:
        return summary
    if "wave_mode" not in ox or ox["wave_mode"].isna().all():
        ox["wave_mode"] = np.where(ox["label"].astype(str).str.endswith("_X"), "X", "O")
    o = ox[ox["wave_mode"].astype(str).str.upper() == "O"]
    x = ox[ox["wave_mode"].astype(str).str.upper() == "X"]
    if o.empty or x.empty:
        return summary
    o0 = o.sort_values("tau_ms").iloc[0]
    x0 = x.iloc[(x["tau_ms"] - o0["tau_ms"]).abs().argmin()]
    idx = summary.index[summary["scenario"] == "OX"]
    if len(idx):
        summary.loc[idx, "delta_tau_OX_us"] = (x0["tau_ms"] - o0["tau_ms"]) * 1000.0
        summary.loc[idx, "delta_group_path_OX_km"] = (
            x0["group_path_km"] - o0["group_path_km"]
        )
    return summary


def _write_tables(data: pd.DataFrame, summary: pd.DataFrame) -> None:
    detail = data[[c for c in DETAIL_COLUMNS if c in data.columns]].copy()
    detail.to_csv(OUT / "report_mode_details.csv", index=False, encoding="utf-8")
    keep = [
        "scenario",
        "scenario_label",
        "n_modes",
        "label",
        "wave_mode",
        "tau_ms",
        "tau_2way_ms",
        "group_path_km",
        "h_reflect_km",
        "Pr_dBW",
        "beta_deg",
        "delta_tau_us",
        "delta_group_path_km",
        "delta_Pr_dB",
        "delta_tau_OX_us",
        "delta_group_path_OX_km",
        "abs_total_dB",
        "scatter_method",
        "sigma_phi_sq",
        "S4",
        "F_focus_3d",
        "A_tube_km2",
    ]
    summary[[c for c in keep if c in summary.columns]].to_csv(
        OUT / "report_scenario_summary.csv", index=False, encoding="utf-8"
    )


def _style_axes(ax):
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _save(fig, name: str) -> None:
    fig.savefig(OUT / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_scenario_comparison(summary: pd.DataFrame) -> None:
    labels = summary["scenario_label"].tolist()
    x = np.arange(len(labels))
    metrics = [
        ("tau_ms", "Delay (ms)"),
        ("group_path_km", "Group path (km)"),
        ("Pr_dBW", "Power (dBW)"),
        ("h_reflect_km", "Reflection height (km)"),
    ]
    colors = ["#0072B2", "#009E73", "#D55E00", "#CC79A7"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.5))
    for ax, (col, ylabel), color in zip(axes.ravel(), metrics, colors):
        ax.bar(x, summary[col], color=color, alpha=0.85)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        _style_axes(ax)
    _save(fig, "report_fig_scenario_comparison.png")


def _plot_delta(summary: pd.DataFrame) -> None:
    rows = summary[summary["scenario"] != "baseline"].copy()
    labels = rows["scenario_label"].tolist()
    x = np.arange(len(labels))
    metrics = [
        ("delta_tau_us", "Delta delay (us)", "#0072B2"),
        ("delta_group_path_km", "Delta group path (km)", "#009E73"),
        ("delta_Pr_dB", "Delta power (dB)", "#D55E00"),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 7.0), sharex=True)
    for ax, (col, ylabel, color) in zip(axes, metrics):
        ax.axhline(0.0, color="0.25", linewidth=0.8)
        ax.bar(x, rows[col], color=color, alpha=0.85)
        ax.set_ylabel(ylabel)
        _style_axes(ax)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, rotation=35, ha="right")
    _save(fig, "report_fig_delta_vs_baseline.png")


def _plot_ox(data: pd.DataFrame) -> None:
    ox = data[data["scenario"] == "OX"].copy()
    if ox.empty:
        return
    if "wave_mode" not in ox or ox["wave_mode"].isna().all():
        ox["wave_mode"] = np.where(ox["label"].astype(str).str.endswith("_X"), "X", "O")
    ox = ox.sort_values(["label", "tau_ms"])
    labels = ox["label"].astype(str).tolist()
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))
    colors = ["#0072B2" if str(w).upper() == "O" else "#D55E00" for w in ox["wave_mode"]]
    axes[0].bar(x, ox["tau_ms"], color=colors, alpha=0.85)
    axes[0].set_ylabel("Delay (ms)")
    axes[1].bar(x, ox["group_path_km"], color=colors, alpha=0.85)
    axes[1].set_ylabel("Group path (km)")
    if len(ox) >= 2:
        o = ox[ox["wave_mode"].astype(str).str.upper() == "O"]
        xx = ox[ox["wave_mode"].astype(str).str.upper() == "X"]
        if not o.empty and not xx.empty:
            d_tau = (xx.iloc[0]["tau_ms"] - o.iloc[0]["tau_ms"]) * 1000.0
            d_gp = xx.iloc[0]["group_path_km"] - o.iloc[0]["group_path_km"]
            axes[2].bar([0, 1], [d_tau, d_gp], color=["#0072B2", "#009E73"], alpha=0.85)
            axes[2].set_xticks([0, 1])
            axes[2].set_xticklabels(["Delay us", "Path km"])
            axes[2].set_ylabel("X - O")
    for ax in axes[:2]:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
    for ax in axes:
        _style_axes(ax)
    _save(fig, "report_fig_ox_split.png")


def _plot_tube(data: pd.DataFrame) -> None:
    tube = data[data["scenario"] == "tube_radar"].sort_values("beta_deg")
    if tube.empty:
        return
    fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.0), sharex=True)
    ycols = [
        ("tau_2way_ms", "Two-way delay (ms)", "#0072B2"),
        ("Pr_dBW", "Power (dBW)", "#D55E00"),
        ("group_path_km", "Group path (km)", "#009E73"),
    ]
    for ax, (col, ylabel, color) in zip(axes, ycols):
        ax.plot(tube["beta_deg"], tube[col], marker="o", color=color, linewidth=1.5)
        ax.set_ylabel(ylabel)
        _style_axes(ax)
    axes[-1].set_xlabel("Elevation angle (deg)")
    _save(fig, "report_fig_tube_radar_scan.png")


def _plot_matrix(summary: pd.DataFrame) -> None:
    rows = summary.copy()
    metrics = ["delta_tau_us", "delta_group_path_km", "delta_Pr_dB", "sigma_phi_sq", "S4"]
    mat = rows[metrics].copy()
    for col in metrics:
        vals = mat[col].astype(float).replace([np.inf, -np.inf], np.nan)
        scale = np.nanmax(np.abs(vals.to_numpy()))
        mat[col] = 0.0 if not scale or math.isnan(scale) else vals / scale
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows["scenario_label"])
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels(["Delay", "Path", "Power", "Sigma phi2", "S4"], rotation=30, ha="right")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if pd.notna(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.82, label="Normalized effect")
    _save(fig, "report_fig_result_matrix.png")


def _write_markdown(data: pd.DataFrame, summary: pd.DataFrame) -> None:
    def md_table(df: pd.DataFrame, cols: list[str]) -> str:
        out = []
        out.append("| " + " | ".join(cols) + " |")
        out.append("| " + " | ".join(["---"] * len(cols)) + " |")
        for _, row in df[cols].iterrows():
            vals = []
            for col in cols:
                val = row[col]
                if pd.isna(val):
                    vals.append("")
                elif isinstance(val, (float, np.floating)):
                    vals.append(f"{float(val):.3f}")
                else:
                    vals.append(str(val))
            out.append("| " + " | ".join(vals) + " |")
        return "\n".join(out)

    lines = [
        "# Enhanced result summary",
        "",
        "This file is generated from `output/modes_*.csv` by `analysis/report_results.py`.",
        "Missing diagnostics are left blank in the CSV tables and are not interpreted as zero.",
        "",
        "## Scenario main-mode summary",
        "",
    ]
    md_cols = [
        "scenario_label",
        "n_modes",
        "label",
        "tau_ms",
        "group_path_km",
        "Pr_dBW",
        "delta_tau_us",
        "delta_group_path_km",
        "delta_Pr_dB",
    ]
    lines.append(md_table(summary, md_cols))
    lines.extend(["", "## O/X split", ""])
    ox = summary[summary["scenario"] == "OX"]
    if ox.empty or pd.isna(ox.iloc[0].get("delta_tau_OX_us", np.nan)):
        lines.append("No complete O/X pair was available.")
    else:
        row = ox.iloc[0]
        lines.append(
            "The X-O delay difference is {:.3f} us and the group-path difference is {:.3f} km.".format(
                row["delta_tau_OX_us"], row["delta_group_path_OX_km"]
            )
        )
    lines.extend(["", "## Tube radar scan", ""])
    tube = data[data["scenario"] == "tube_radar"]
    if tube.empty:
        lines.append("No tube-radar rows were available.")
    else:
        lines.append(
            "Tube radar covers {:.1f}-{:.1f} deg elevation, {:.3f}-{:.3f} ms two-way delay, and {:.2f}-{:.2f} dBW received power.".format(
                tube["beta_deg"].min(),
                tube["beta_deg"].max(),
                tube["tau_2way_ms"].min(),
                tube["tau_2way_ms"].max(),
                tube["Pr_dBW"].min(),
                tube["Pr_dBW"].max(),
            )
        )
    lines.extend(
        [
            "",
            "## Generated figures",
            "",
            "- `report_fig_scenario_comparison.png`",
            "- `report_fig_delta_vs_baseline.png`",
            "- `report_fig_ox_split.png`",
            "- `report_fig_tube_radar_scan.png`",
            "- `report_fig_result_matrix.png`",
            "",
            "## Physical notes",
            "",
            "`group_path_km` is the cold-plasma group path integral `integral(ds/n)`, not the geometric path length. O/X splitting is post-processed along 3-D isotropic trajectories and is not a full 3-D Appleton-Hartree ray trace. Scatter coupling diagnostics attach Born/Rytov/MPS regime labels without changing `Pr_W` by default.",
        ]
    )
    (OUT / "report_results_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    data = _read_modes()
    if data.empty:
        raise SystemExit("No modes_*.csv files found in output")
    summary = _main_modes(data)
    summary = _add_ox_metrics(data, summary)
    _write_tables(data, summary)
    _plot_scenario_comparison(summary)
    _plot_delta(summary)
    _plot_ox(data)
    _plot_tube(data)
    _plot_matrix(summary)
    _write_markdown(data, summary)
    print("OK report_results")


if __name__ == "__main__":
    main()
