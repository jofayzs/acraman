"""Cross-check: the simulated success curve must overlay the exact one."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OverlayCheck:
    passed: bool
    k_se: float
    n_checked: int
    max_gap_se: float      # largest |sim - exact| in standard-error units
    max_gap_se_at: int
    max_gap_pp: float      # largest |sim - exact| in percentage points
    max_gap_pp_at: int


def check_points(n_cells: int, anchors: list[int | None], log_points: int,
                 window_fraction: float) -> np.ndarray:
    """~log_points log-spaced n over the slide plus every n within +/- window of each anchor."""
    pts = [np.unique(np.round(np.geomspace(1, n_cells, log_points)).astype(np.int64))]
    for a in anchors:
        if a is not None:
            lo = max(1, int(a * (1 - window_fraction)))
            hi = min(n_cells, int(np.ceil(a * (1 + window_fraction))))
            pts.append(np.arange(lo, hi + 1, dtype=np.int64))
    return np.unique(np.concatenate(pts))


def overlay_check(sim: np.ndarray, exact: np.ndarray, runs: int, n_values: np.ndarray,
                  k_se: float) -> OverlayCheck:
    """Pass iff |sim - exact| <= k_se * SE at every checked n.

    SE = sqrt(max(p(1-p), 1/R) / R); the 1/R floor avoids false failures where p is ~0 or ~1.
    """
    idx = n_values - 1
    p = exact[idx]
    gap = np.abs(sim[idx] - p)
    se = np.sqrt(np.maximum(p * (1 - p), 1 / runs) / runs)
    gap_se = gap / se
    i_se = int(np.argmax(gap_se))
    i_pp = int(np.argmax(gap))
    return OverlayCheck(
        passed=bool(np.all(gap_se <= k_se)), k_se=k_se, n_checked=int(n_values.size),
        max_gap_se=float(gap_se[i_se]), max_gap_se_at=int(n_values[i_se]),
        max_gap_pp=float(gap[i_pp] * 100), max_gap_pp_at=int(n_values[i_pp]),
    )
