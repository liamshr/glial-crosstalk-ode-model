"""Generate baseline figures: python scripts/generate_figures.py."""

from pathlib import Path

import matplotlib.pyplot as plt

from glial_crosstalk.core import simulate
from glial_crosstalk.parameters import DEFAULT_PARAMETERS
from glial_crosstalk.plotting import phase_portrait, trajectory_figure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"


def main():
    solution = simulate(params=DEFAULT_PARAMETERS)
    if not solution.success:
        raise RuntimeError(solution.message)
    trajectory_figure(solution, OUT / "baseline_trajectories.png")
    trajectory_figure(solution, OUT / "baseline_trajectories.pdf")
    phase_portrait(solution, output=OUT / "baseline_M1_D_phase_portrait.png")
    phase_portrait(solution, output=OUT / "baseline_M1_D_phase_portrait.pdf")
    for ros in (0.0, 1.0, 100.0):
        sol = simulate(params=DEFAULT_PARAMETERS.updated(ROS=ros))
        phase_portrait(sol, output=OUT / f"ros_{ros:g}_M1_D_phase_portrait.png")
    plt.close("all")


if __name__ == "__main__":
    main()
