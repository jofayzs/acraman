"""Exact solution: hits after n zaps without replacement are hypergeometric."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy.stats import hypergeom

from .model import hit_bounds


def success_probability(n_cells: int, n_spherules: int, band: tuple[Fraction, Fraction],
                        hits_detectable: bool = True) -> np.ndarray:
    """P(estimate in band) for n = 1..n_cells. Index i holds n = i + 1."""
    if n_spherules == 0 or not hits_detectable:
        return np.zeros(n_cells)
    n = np.arange(1, n_cells + 1, dtype=np.int64)
    kmin, kmax = hit_bounds(n, n_spherules, n_cells, band)
    kmax = np.minimum(kmax, np.minimum(n, n_spherules))
    p = (hypergeom.cdf(kmax, n_cells, n_spherules, n)
         - hypergeom.cdf(kmin - 1, n_cells, n_spherules, n))
    p[kmax < kmin] = 0.0
    return np.clip(p, 0.0, 1.0)


def first_crossing(p: np.ndarray, confidence: float) -> int | None:
    """Smallest n with p >= confidence."""
    ok = np.flatnonzero(p >= confidence)
    return int(ok[0]) + 1 if ok.size else None


def stable_n(p: np.ndarray, confidence: float) -> int | None:
    """Smallest n with p >= confidence at n and at every larger n up to the whole slide."""
    if p.size == 0 or p[-1] < confidence:
        return None
    below = np.flatnonzero(p < confidence)
    return int(below[-1]) + 2 if below.size else 1
