"""One two-panel figure per concentration."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"     # simulated
ORANGE = "#eb6834"   # exact
BAND = "#e1e0d9"

RATIO_YMAX = 2.5
MAX_POINTS = 4000


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def plot_concentration(path: Path, *, label: str, realized_label: str, n_cells: int,
                       trajectories: list[np.ndarray], sim: np.ndarray, exact: np.ndarray,
                       band: tuple[float, float], confidence: float, runs: int,
                       first: int | None, stable: int | None) -> None:
    """trajectories: estimate / realized truth for n = 1..n_cells, one array per run."""
    x_max = n_cells if stable is None else min(n_cells, 2 * stable)
    n = np.unique(np.linspace(1, x_max, MAX_POINTS).round().astype(np.int64))
    idx = n - 1

    fig, (top, bot) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, facecolor=SURFACE,
                                   gridspec_kw={"height_ratios": [1, 1], "hspace": 0.2})

    _style(top)
    top.axhspan(band[0], band[1], color=BAND, alpha=0.7, linewidth=0)
    top.axhline(1.0, color=MUTED, linewidth=0.8)
    for t in trajectories:
        top.plot(n, np.minimum(t[idx], RATIO_YMAX * 1.5), color=BLUE, linewidth=0.8, alpha=0.45)
    top.set_ylim(0, RATIO_YMAX)
    top.set_ylabel("Estimate ÷ true", color=INK_2, fontsize=10)
    top.text(x_max, band[1], f" {band[0]:g}–{band[1]:g} band", va="bottom", ha="right",
             color=INK_2, fontsize=9)
    top.set_title(f"{len(trajectories)} individual simulated runs (clipped at {RATIO_YMAX:g}×)",
                  loc="left", color=INK_2, fontsize=10)

    _style(bot)
    bot.plot(n, sim[idx], color=BLUE, linewidth=2, label=f"Simulated ({runs:,} runs)")
    bot.plot(n, exact[idx], color=ORANGE, linewidth=1.2, label="Exact (hypergeometric)")
    bot.axhline(confidence, color=MUTED, linewidth=0.8, linestyle=(0, (4, 3)))
    bot.annotate(f"{confidence:.0%}", (0, confidence), xytext=(5, 3),
                 textcoords="offset points", color=INK_2, fontsize=9)
    same = first is not None and first == stable
    if first is not None and not same:
        bot.axvline(first, color=MUTED, linewidth=0.8, linestyle=(0, (1, 2)))
        bot.annotate(f"first crossing\n{first:,}", (first, 0.05), xytext=(-5, 0),
                     textcoords="offset points", ha="right", va="bottom",
                     color=INK_2, fontsize=9)
    if stable is not None:
        for ax in (top, bot):
            ax.axvline(stable, color=INK_2, linewidth=1)
        place_left = stable > 0.75 * x_max
        text = f"required (stable) n\n{stable:,}" + (" = first crossing" if same else "")
        bot.annotate(text, (stable, 0.05), xytext=(-5 if place_left else 5, 0),
                     textcoords="offset points", ha="right" if place_left else "left",
                     va="bottom", color=INK, fontsize=9, fontweight="bold")
    bot.set_ylim(0, 1.02)
    bot.set_xlim(0, x_max)
    bot.set_ylabel("P(estimate in band)", color=INK_2, fontsize=10)
    bot.set_xlabel("Number of zaps", color=INK_2, fontsize=10)
    bot.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
    bot.legend(loc="lower left", bbox_to_anchor=(-0.01, 0.98), ncol=2, frameon=False,
               fontsize=9, labelcolor=INK_2, borderaxespad=0)

    fig.suptitle(f"True [Ir] {label} nominal ({realized_label} realized)", x=0.125, ha="left",
                 color=INK, fontsize=13, fontweight="bold")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
