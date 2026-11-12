"""Monte Carlo simulation of zapping runs.

Two paths give the same hits(n):

* reference path: the literal model (build a slide, permute the cells, read
  every zap). Slow, used for plotted trajectories and as the place v2
  features plug in.
* fast path: a run is fully described by the zap numbers at which its
  spherules get hit, so only those K numbers are drawn.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np

from . import model


# ---------- reference path ----------

def reference_hits(slide: np.ndarray, order: np.ndarray, hit_threshold_ppb: float) -> np.ndarray:
    """Cumulative hits after n = 1..len(order) zaps, going through the model step by step."""
    readings = model.read_zaps(slide, order)
    return np.cumsum(model.is_hit(readings, hit_threshold_ppb))


def reference_run(n_cells: int, n_spherules: int, grade_ppb: float, hit_threshold_ppb: float,
                  rng: np.random.Generator) -> np.ndarray:
    slide = model.build_slide(n_cells, n_spherules, grade_ppb, rng)
    order = model.zap_order(n_cells, rng)
    return reference_hits(slide, order, hit_threshold_ppb)


# ---------- fast path ----------

def hit_positions(slide: np.ndarray, order: np.ndarray, hit_threshold_ppb: float) -> np.ndarray:
    """Sorted 1-based zap numbers at which a hit occurs, for a given slide and zap order."""
    return np.flatnonzero(model.is_hit(model.read_zaps(slide, order), hit_threshold_ppb)) + 1


def hits_from_positions(positions: np.ndarray, n_cells: int) -> np.ndarray:
    """Cumulative hits after n = 1..n_cells zaps, from sorted hit positions."""
    return np.searchsorted(positions, np.arange(1, n_cells + 1), side="right")


def draw_positions(n_cells: int, n_spherules: int, runs: int,
                   rng: np.random.Generator) -> np.ndarray:
    """(runs, K) sorted zap numbers of the spherule hits; each row is one run."""
    pos = np.empty((runs, n_spherules), dtype=np.int64)
    for r in range(runs):
        pos[r] = rng.choice(n_cells, size=n_spherules, replace=False)
    pos += 1
    pos.sort(axis=1)
    return pos


def _band_n_ranges(n_cells: int, n_spherules: int,
                   band: tuple[Fraction, Fraction]) -> tuple[np.ndarray, np.ndarray]:
    """For each hit count h = 0..K, the inclusive range of n at which h hits is in band."""
    lo, hi = band
    h = np.arange(n_spherules + 1, dtype=np.int64)
    # hi.den*h*N <= hi.num*K*n  ->  n >= ceil(hi.den*h*N / (hi.num*K))
    n_lo = -((-hi.denominator * h * n_cells) // (hi.numerator * n_spherules))
    # lo.num*K*n <= lo.den*h*N  ->  n <= floor(lo.den*h*N / (lo.num*K))
    n_hi = (lo.denominator * h * n_cells) // (lo.numerator * n_spherules)
    return np.maximum(n_lo, 1), np.minimum(n_hi, n_cells)


def success_counts(positions: np.ndarray, n_cells: int, n_spherules: int,
                   band: tuple[Fraction, Fraction]) -> np.ndarray:
    """Number of runs in band at n = 1..n_cells, given (runs, K) sorted hit positions."""
    runs = positions.shape[0]
    n_lo, n_hi = _band_n_ranges(n_cells, n_spherules, band)
    # A run holds h hits for n in [seg_start[h], seg_end[h]].
    seg_start = np.concatenate([np.ones((runs, 1), dtype=np.int64), positions], axis=1)
    seg_end = np.concatenate(
        [positions - 1, np.full((runs, 1), n_cells, dtype=np.int64)], axis=1)
    start = np.maximum(seg_start, n_lo)
    end = np.minimum(seg_end, n_hi)
    valid = start <= end
    diff = (np.bincount(start[valid], minlength=n_cells + 2)
            - np.bincount(end[valid] + 1, minlength=n_cells + 2))
    return np.cumsum(diff)[1:n_cells + 1]


def simulate_success_fraction(n_cells: int, n_spherules: int, band: tuple[Fraction, Fraction],
                              runs: int, rng: np.random.Generator,
                              hits_detectable: bool = True, chunk: int = 500) -> np.ndarray:
    """Fraction of simulated runs in band at n = 1..n_cells. Index i holds n = i + 1."""
    if n_spherules == 0 or not hits_detectable:
        return np.zeros(n_cells)
    counts = np.zeros(n_cells, dtype=np.int64)
    done = 0
    while done < runs:
        size = min(chunk, runs - done)
        counts += success_counts(draw_positions(n_cells, n_spherules, size, rng),
                                 n_cells, n_spherules, band)
        done += size
    return counts / runs
