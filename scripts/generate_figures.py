"""Generate reproducible publication figures for the model."""

from pathlib import Path

import matplotlib.pyplot as plt

from glial_crosstalk.core import simulate
from glial_crosstalk.parameters import DEFAULT_PARAMETERS
from glial_crosstalk.plotting import trajectory_figure
from euler_microglial_polarization import run_euler_simulation
from sensitivity_analysis import plot_sensitivity, run_sensitivity_analysis

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
ROS_VALUES = (0.00, 0.01, 31.62, 100000.00)


def save_euler_figure(output):
    """Save the standalone Euler microglial polarization example."""
    time, m1_trace, m2_trace = run_euler_simulation()
    figure, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    ax.plot(time, m1_trace, label="M1 (pro-inflammatory)")
    ax.plot(time, m2_trace, label="M2 (anti-inflammatory)")
    ax.set_xlabel("Time (arbitrary units)")
    ax.set_ylabel("Concentration (arbitrary units)")
    ax.set_title("Euler simulation of microglial polarization (ROS=2.0)")
    ax.legend(frameon=False)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    return figure


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for ros in ROS_VALUES:
        sol = simulate(params=DEFAULT_PARAMETERS.updated(ROS=ros))
        if not sol.success:
            raise RuntimeError(f"Solver failed for ROS={ros}: {sol.message}")
        stem = f"ros_{ros:.2f}_trajectories"
        title = f"Model trajectories with constant ROS input = {ros:.2f}"
        trajectory_figure(sol, OUT / f"{stem}.png", title=title)
        trajectory_figure(sol, OUT / f"{stem}.pdf", title=title)

    _, relative_change, min_effective, saturation_ros, _ = run_sensitivity_analysis()
    for suffix in ("png", "pdf"):
        figure = plot_sensitivity(
            relative_change,
            min_effective,
            saturation_ros,
            OUT / f"sensitivity_DVs_to_ROS.{suffix}",
        )
        plt.close(figure)

    for suffix in ("png", "pdf"):
        figure = save_euler_figure(OUT / f"euler_microglial_polarization.{suffix}")
        plt.close(figure)

    plt.close("all")


if __name__ == "__main__":
    main()
