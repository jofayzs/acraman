"""Load and validate run parameters."""

from __future__ import annotations

import math
import tomllib
import warnings
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path


def _frac(x) -> Fraction:
    """Exact fraction from a decimal as written in the config (0.8 -> 4/5, not the binary float)."""
    return Fraction(str(x))


@dataclass(frozen=True)
class Config:
    slide_width_cm: float = 3.0
    slide_height_cm: float = 2.0
    spot_diameter_um: float = 60.0
    grade_ppb: float = 100.0
    hit_threshold_ppb: float = 100.0
    concentrations_ppb: tuple[float, ...] = (1.0, 0.5, 0.001)
    band_low: float = 0.8
    band_high: float = 1.2
    confidence: float = 0.90
    zaps_per_hour: float = 15.0
    usd_per_hour: float = 100.0
    feasible_max_hours: float = 100.0
    infeasible_min_hours: float = 500.0
    infeasible_slide_fraction: float = 0.5
    runs: int = 10_000
    seed: int = 20261111
    trajectories: int = 20
    k_se: float = 5.0
    log_points: int = 200
    window_fraction: float = 0.10

    def __post_init__(self) -> None:
        if not 0 < self.band_low <= 1 <= self.band_high:
            raise ValueError("band must satisfy 0 < band_low <= 1 <= band_high")
        if not 0 < self.confidence < 1:
            raise ValueError("confidence must be between 0 and 1")
        if self.grade_ppb <= 0 or self.runs <= 0:
            raise ValueError("grade_ppb and runs must be positive")
        if self.hit_threshold_ppb > self.grade_ppb:
            warnings.warn(
                f"hit_threshold_ppb ({self.hit_threshold_ppb}) is above grade_ppb "
                f"({self.grade_ppb}): no zap can ever count as a hit.",
                stacklevel=2,
            )

    @property
    def total_cells(self) -> int:
        slide_area_cm2 = self.slide_width_cm * self.slide_height_cm
        radius_cm = self.spot_diameter_um / 2 * 1e-4
        return round(slide_area_cm2 / (math.pi * radius_cm**2))

    @property
    def band(self) -> tuple[Fraction, Fraction]:
        return _frac(self.band_low), _frac(self.band_high)

    @property
    def hits_detectable(self) -> bool:
        return self.hit_threshold_ppb <= self.grade_ppb

    def n_spherules(self, concentration_ppb: float) -> int:
        return round(_frac(concentration_ppb) / _frac(self.grade_ppb) * self.total_cells)

    def realized_ppb(self, concentration_ppb: float) -> float:
        return self.n_spherules(concentration_ppb) * self.grade_ppb / self.total_cells

    def as_dict(self) -> dict:
        d = asdict(self)
        d["concentrations_ppb"] = list(self.concentrations_ppb)
        return d


_SECTIONS = {
    "slide": {"width_cm": "slide_width_cm", "height_cm": "slide_height_cm",
              "spot_diameter_um": "spot_diameter_um"},
    "spherules": {"grade_ppb": "grade_ppb"},
    "measurement": {"hit_threshold_ppb": "hit_threshold_ppb"},
    "target": {"concentrations_ppb": "concentrations_ppb"},
    "success": {"band_low": "band_low", "band_high": "band_high", "confidence": "confidence"},
    "cost": {"zaps_per_hour": "zaps_per_hour", "usd_per_hour": "usd_per_hour"},
    "feasibility": {"feasible_max_hours": "feasible_max_hours",
                    "infeasible_min_hours": "infeasible_min_hours",
                    "infeasible_slide_fraction": "infeasible_slide_fraction"},
    "simulation": {"runs": "runs", "seed": "seed", "trajectories": "trajectories"},
    "validation": {"k_se": "k_se", "log_points": "log_points",
                   "window_fraction": "window_fraction"},
}


def load_config(path: str | Path) -> Config:
    with open(path, "rb") as f:
        raw = tomllib.load(f)
    kwargs = {}
    for section, body in raw.items():
        if section not in _SECTIONS:
            raise ValueError(f"unknown config section [{section}]")
        for key, value in body.items():
            if key not in _SECTIONS[section]:
                raise ValueError(f"unknown config key {section}.{key}")
            kwargs[_SECTIONS[section][key]] = value
    if "concentrations_ppb" in kwargs:
        kwargs["concentrations_ppb"] = tuple(kwargs["concentrations_ppb"])
    return Config(**kwargs)


def concentration_label(ppb: float) -> str:
    """'1 ppb', '0.5 ppb', '1 ppt'."""
    if ppb >= 0.1:
        return f"{ppb:g} ppb"
    return f"{ppb * 1000:g} ppt"
