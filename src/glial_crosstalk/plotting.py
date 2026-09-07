"""Publication-oriented plots for model trajectories and phase portraits."""

from pathlib import Path

import matplotlib.pyplot as plt

from .core import STATE_NAMES, rates


def trajectory_figure(solution, output=None, *, title=None, dpi=300):
    """Plot all biological groups and auxiliary recruitment rates."""
    fig, axes = plt.subplots(4, 1, figsize=(9, 12), sharex=True, constrained_layout=True)
    groups = [
        (("M1", "M2"), "Microglia", ("#c0392b", "#2980b9")),
        (("IL1", "IL12", "IL10", "IL4"), "Cytokines", ("#e67e22", "#f39c12", "#3498db", "#1abc9c")),
        (("Aq", "Ap"), "Astroglia", ("#16a085", "#8e44ad")),
        (("D",), "Damage signal", ("#2c3e50",)),
    ]
    for ax, (names, title, colors) in zip(axes, groups):
        for name, color in zip(names, colors):
            ax.plot(solution.t, solution.y[STATE_NAMES.index(name)], label=name, lw=2, color=color)
        ax.set_ylabel("a.u.")
        ax.set_title(title)
        ax.legend(ncol=len(names), frameon=False)
    axes[-1].set_xlabel("Time (a.u.)")
    if title is not None:
        fig.suptitle(title)
        fig.set_constrained_layout_pads(h_pad=0.25)
    return _save_or_return(fig, output, dpi)


def phase_portrait(solution, x="M1", y="D", output=None, *, dpi=300):
    """Plot a two-state phase portrait from one trajectory."""
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    ax.plot(solution.y[STATE_NAMES.index(x)], solution.y[STATE_NAMES.index(y)], lw=2)
    ax.scatter(solution.y[STATE_NAMES.index(x), 0], solution.y[STATE_NAMES.index(y), 0],
               color="tab:green", label="initial")
    ax.scatter(solution.y[STATE_NAMES.index(x), -1], solution.y[STATE_NAMES.index(y), -1],
               color="tab:red", label="final")
    ax.set(xlabel=f"{x} (a.u.)", ylabel=f"{y} (a.u.)", title=f"{x}–{y} phase portrait")
    ax.legend(frameon=False)
    return _save_or_return(fig, output, dpi)


def _save_or_return(fig, output, dpi):
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
    return fig
