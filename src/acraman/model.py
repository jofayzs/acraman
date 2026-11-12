"""The v1 physical model as small swappable steps.

Later versions replace individual steps: `build_slide` for spherule sizes and
background Ir, `zap_order` for targeted sampling or real circle geometry,
`read_zaps` for instrument noise and detection limits, `estimate` for bias
corrections. The reference simulation path in `simulate.py` runs through
exactly these functions.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np


def build_slide(n_cells: int, n_spherules: int, grade_ppb: float,
                rng: np.random.Generator) -> np.ndarray:
    """Per-cell Ir concentration: spherules in random cells at `grade_ppb`, the rest 0."""
    slide = np.zeros(n_cells)
    slide[rng.choice(n_cells, size=n_spherules, replace=False)] = grade_ppb
    return slide


def zap_order(n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Cells in the order they are zapped: uniform random, never the same cell twice."""
    return rng.permutation(n_cells)


def read_zaps(slide: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Instrument reading for each zapped cell. v1: the cell's true value, no noise."""
    return slide[cells]


def is_hit(readings: np.ndarray, hit_threshold_ppb: float) -> np.ndarray:
    return readings >= hit_threshold_ppb


def estimate(hits: np.ndarray, n_zaps: np.ndarray, grade_ppb: float) -> np.ndarray:
    """Bulk concentration estimate after n zaps: (hits / n) x grade."""
    return hits / n_zaps * grade_ppb


def in_band(hits: int, n_zaps: int, n_spherules: int, n_cells: int,
            band: tuple[Fraction, Fraction]) -> bool:
    """Is estimate / realized truth inside the band (inclusive)?

    estimate / truth = (hits / n) / (K / N) = hits*N / (K*n), compared in
    integers so ties on the band edge are decided exactly.
    """
    if n_spherules == 0:
        return False
    lo, hi = band
    lhs = hits * n_cells
    rhs = n_spherules * n_zaps
    return lo.numerator * rhs <= lo.denominator * lhs and hi.denominator * lhs <= hi.numerator * rhs


def hit_bounds(n_zaps: np.ndarray, n_spherules: int, n_cells: int,
               band: tuple[Fraction, Fraction]) -> tuple[np.ndarray, np.ndarray]:
    """Smallest and largest hit counts that are in band at each n (vectorised `in_band`)."""
    lo, hi = band
    n = np.asarray(n_zaps, dtype=np.int64)
    kmin = -((-lo.numerator * n_spherules * n) // (lo.denominator * n_cells))
    kmax = (hi.numerator * n_spherules * n) // (hi.denominator * n_cells)
    return kmin, kmax
