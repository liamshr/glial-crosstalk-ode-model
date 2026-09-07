#!/usr/bin/env python3
"""Analyze steady-state model sensitivity to ROS.

This script extracts the notebook analysis that produces the
"Sensitivity of DVs to ROS" graph.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from glial_crosstalk.core import STATE_NAMES, simulate
from glial_crosstalk.parameters import DEFAULT_PARAMETERS


def summarize_solution(solution, dv_names):
    """Summarize final-state, peak, area, and peak-time metrics."""
    t = solution.t
    y = solution.y
    last_n = max(1, int(np.ceil(len(t) * 0.05)))
    last_indices = np.arange(len(t) - last_n, len(t))

    summary = {}
    for index, name in enumerate(dv_names):
        series = y[index, :]
        summary[name] = {
            "steady_state": np.mean(series[last_indices]),
            "peak": np.max(series),
            "auc": np.trapezoid(series, t),
            "time_to_peak": t[np.argmax(series)],
        }
    return summary


def run_sensitivity_analysis(
    *,
    ros_grid=None,
    min_change_tol=0.10,
    sat_tol=0.01,
    max_ros_allowed=1e6,
):
    """Run the ROS sweep and return metrics, relative changes, and thresholds."""
    if ros_grid is None:
        ros_grid = np.unique(
            np.concatenate((np.array([0.0]), np.logspace(-2, 5, 60)))
        )
    ros_grid = np.asarray(ros_grid, dtype=float)
    ros_grid = ros_grid[ros_grid <= max_ros_allowed]

    dv_names = list(STATE_NAMES)
    metrics = []
    for ros in ros_grid:
        params = DEFAULT_PARAMETERS.updated(ROS=float(ros))
        solution = simulate(params=params)
        if not solution.success:
            raise RuntimeError(f"Solver failed for ROS={ros}: {solution.message}")
        metrics.append(
            {
                "ROS": ros,
                "summary": summarize_solution(solution, dv_names),
            }
        )

    rows = []
    for record in metrics:
        row = {"ROS": record["ROS"]}
        for dv in dv_names:
            summary = record["summary"][dv]
            row[f"{dv}_ss"] = summary["steady_state"]
            row[f"{dv}_peak"] = summary["peak"]
            row[f"{dv}_auc"] = summary["auc"]
            row[f"{dv}_ttp"] = summary["time_to_peak"]
        rows.append(row)
    df_metrics = pd.DataFrame(rows).sort_values("ROS").reset_index(drop=True)

    control_rows = df_metrics.loc[df_metrics["ROS"] == 0.0]
    if control_rows.empty:
        raise RuntimeError("Control (ROS=0) not found in metrics.")
    control = control_rows.iloc[0]

    eps = 1e-12
    relative_change = pd.DataFrame({"ROS": df_metrics["ROS"]})
    for dv in dv_names:
        base = control[f"{dv}_ss"]
        denominator = base if abs(base) > eps else abs(base) + eps
        relative_change[f"{dv}_ss_rel"] = (
            df_metrics[f"{dv}_ss"] - base
        ) / denominator
    relative_columns = [f"{dv}_ss_rel" for dv in dv_names]
    relative_change["max_abs_ss_rel"] = relative_change[relative_columns].abs().max(axis=1)

    candidates = relative_change[
        relative_change["max_abs_ss_rel"] >= min_change_tol
    ]
    min_effective = (
        float(candidates["ROS"].iloc[0]) if not candidates.empty else None
    )

    saturation_ros = None
    for index in range(1, len(df_metrics)):
        changes = []
        start = df_metrics.iloc[index]
        for later_index in range(index, len(df_metrics)):
            later = df_metrics.iloc[later_index]
            dv_changes = []
            for dv in dv_names:
                start_value = start[f"{dv}_ss"]
                denominator = (
                    start_value if abs(start_value) > eps else start_value + eps
                )
                dv_changes.append(
                    abs((later[f"{dv}_ss"] - start_value) / denominator)
                )
            changes.append(max(dv_changes))
        if max(changes) < sat_tol:
            saturation_ros = float(start["ROS"])
            break

    if (
        min_effective is not None
        and saturation_ros is not None
        and saturation_ros > min_effective
    ):
        mid_ros = float(np.sqrt(min_effective * saturation_ros))
    else:
        nonzero = df_metrics.loc[df_metrics["ROS"] > 0, "ROS"].to_numpy()
        mid_ros = float(nonzero[len(nonzero) // 2]) if len(nonzero) else None

    return (
        df_metrics,
        relative_change,
        min_effective,
        saturation_ros,
        mid_ros,
    )


def plot_sensitivity(relative_change, min_effective, saturation_ros, output=None):
    """Create the notebook's sensitivity graph and optionally save it."""
    figure = plt.figure(figsize=(8, 5))
    plt.semilogx(
        relative_change["ROS"].replace(0, np.nan).fillna(1e-6),
        relative_change["max_abs_ss_rel"],
        marker="o",
    )
    min_change_tol = 0.10
    plt.axhline(
        min_change_tol,
        color="grey",
        linestyle="--",
        label=f"{min_change_tol * 100:.0f}% threshold",
    )
    if min_effective is not None:
        plt.axvline(
            min_effective,
            color="green",
            linestyle=":",
            label=f"min_effective={min_effective}",
        )
    if saturation_ros is not None:
        plt.axvline(
            saturation_ros,
            color="red",
            linestyle=":",
            label=f"saturation={saturation_ros}",
        )
    plt.xlabel("ROS (log scale)")
    plt.ylabel("max absolute relative change in steady-state vs control")
    plt.title("Sensitivity of DVs to ROS")
    plt.legend()
    plt.tight_layout()
    if output is not None:
        figure.savefig(output, dpi=300, bbox_inches="tight")
    return figure


def main():
    parser = argparse.ArgumentParser(
        description="Generate the model's Sensitivity of DVs to ROS graph."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("figures/sensitivity_DVs_to_ROS.png"),
        help="Path for the generated PNG.",
    )
    args = parser.parse_args()

    _, relative_change, min_effective, saturation_ros, mid_ros = (
        run_sensitivity_analysis()
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot_sensitivity(relative_change, min_effective, saturation_ros, args.output)
    plt.close("all")

    print(f"min_effective ROS: {min_effective}")
    print(f"saturation ROS: {saturation_ros}")
    print(f"mid ROS: {mid_ros}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
