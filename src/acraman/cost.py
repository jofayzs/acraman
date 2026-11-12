"""Turn a zap count into slide area, machine time, money and a feasibility verdict."""

from __future__ import annotations

from dataclasses import dataclass

from .config import Config


@dataclass(frozen=True)
class Cost:
    zaps: int
    slide_fraction: float
    hours: float
    usd: float


def cost_of(zaps: int, cfg: Config) -> Cost:
    hours = zaps / cfg.zaps_per_hour
    return Cost(zaps=zaps, slide_fraction=zaps / cfg.total_cells, hours=hours,
                usd=hours * cfg.usd_per_hour)


@dataclass(frozen=True)
class Verdict:
    tier: str      # Feasible | Costly | Infeasible | Impossible
    reason: str    # which constraint triggered the tier


def verdict_for(zaps: int | None, cfg: Config) -> Verdict:
    """Feasibility tier for the headline (stable) zap count; None means no n on the slide works."""
    if zaps is None:
        return Verdict("Impossible",
                       "no number of zaps up to the whole slide reaches the confidence level")
    c = cost_of(zaps, cfg)
    reasons = []
    if c.hours > cfg.infeasible_min_hours:
        reasons.append(f"{c.hours:,.0f} machine hours exceeds the {cfg.infeasible_min_hours:g} h limit")
    if c.slide_fraction > cfg.infeasible_slide_fraction:
        reasons.append(f"{c.slide_fraction:.1%} of the slide exceeds the "
                       f"{cfg.infeasible_slide_fraction:.0%} limit")
    if reasons:
        return Verdict("Infeasible", " and ".join(reasons))
    if c.hours > cfg.feasible_max_hours:
        return Verdict("Costly", f"{c.hours:,.0f} machine hours is between the "
                                 f"{cfg.feasible_max_hours:g} h and {cfg.infeasible_min_hours:g} h limits")
    return Verdict("Feasible", f"{c.hours:,.0f} machine hours is within the "
                               f"{cfg.feasible_max_hours:g} h limit")
