"""Write summary.csv, report.md and params.json."""

from __future__ import annotations

import csv
import json
import platform
from dataclasses import asdict, dataclass
from importlib.metadata import version
from pathlib import Path

from .config import Config
from .cost import Cost, Verdict, cost_of
from .validate import OverlayCheck


@dataclass(frozen=True)
class Result:
    nominal_ppb: float
    label: str
    realized_ppb: float
    n_spherules: int
    exact_first: int | None
    exact_stable: int | None
    sim_first: int | None
    sim_stable: int | None
    verdict: Verdict
    check: OverlayCheck
    figure: str

    def cost(self, cfg: Config) -> Cost | None:
        return None if self.exact_stable is None else cost_of(self.exact_stable, cfg)


def conditions(cfg: Config) -> str:
    return (f"blind random sampling with {cfg.grade_ppb:g} ppb one-cell spherules")


def best_case_statement(cfg: Config) -> str:
    return ("Best case for blind (uniform random) sampling, assuming every spherule is exactly "
            f"one spot in size at {cfg.grade_ppb:g} ppb and the rest of the rock is 0 ppb.")


def _n(x: int | None) -> str:
    return "none" if x is None else f"{x:,}"


def _ppb(x: float) -> str:
    return f"{x:.6g} ppb"


def write_csv(path: Path, results: list[Result], cfg: Config) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nominal_ppb", "realized_ppb", "spherules", "required_zaps_stable_exact",
                    "first_crossing_exact", "pct_of_slide", "machine_hours", "cost_usd",
                    "verdict", "verdict_trigger", "sim_stable_informational",
                    "sim_first_crossing_informational", "overlay_check_passed",
                    "max_gap_se", "max_gap_pct_points"])
        for r in results:
            c = r.cost(cfg)
            w.writerow([r.nominal_ppb, f"{r.realized_ppb:.9g}", r.n_spherules,
                        r.exact_stable, r.exact_first,
                        "" if c is None else f"{c.slide_fraction * 100:.2f}",
                        "" if c is None else f"{c.hours:.1f}",
                        "" if c is None else f"{c.usd:.0f}",
                        r.verdict.tier, r.verdict.reason, r.sim_stable, r.sim_first,
                        r.check.passed, f"{r.check.max_gap_se:.2f}", f"{r.check.max_gap_pp:.2f}"])


def write_params(path: Path, cfg: Config, config_path: str) -> None:
    record = {
        "config_file": config_path,
        "parameters": cfg.as_dict(),
        "derived": {
            "total_cells": cfg.total_cells,
            "spherules": {str(c): cfg.n_spherules(c) for c in cfg.concentrations_ppb},
            "realized_ppb": {str(c): cfg.realized_ppb(c) for c in cfg.concentrations_ppb},
        },
        "versions": {"python": platform.python_version(),
                     **{p: version(p) for p in ("acraman", "numpy", "scipy", "matplotlib")}},
    }
    path.write_text(json.dumps(record, indent=2) + "\n")


def write_report(path: Path, results: list[Result], cfg: Config) -> None:
    all_passed = all(r.check.passed for r in results)
    L: list[str] = []
    L += ["# Ir zap simulator v1: results", "",
          f"**{best_case_statement(cfg)}**", "",
          f"Slide: {cfg.slide_width_cm:g} cm × {cfg.slide_height_cm:g} cm, "
          f"{cfg.spot_diameter_um:g} µm spot, treated as {cfg.total_cells:,} spot-sized cells. "
          f"A run succeeds at n zaps when its estimate is within {cfg.band_low:g}–{cfg.band_high:g}× "
          "the slide's realized concentration (edges included). "
          f"\"Required zaps\" is the stable n: the smallest n at which at least {cfg.confidence:.0%} "
          "of runs succeed at that n and at every larger n. All headline numbers are exact "
          "(hypergeometric); the simulation is a cross-check.", ""]

    L += ["## Summary", "",
          "| Nominal [Ir] | Realized [Ir] | Spherules | Required zaps (stable n) | First crossing "
          "| % of slide | Machine hours | Cost | Verdict |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        c = r.cost(cfg)
        cells = (["n/a"] * 3 if c is None else
                 [f"{c.slide_fraction:.1%}", f"{c.hours:,.0f}", f"${c.usd:,.0f}"])
        L.append(f"| {r.label} | {_ppb(r.realized_ppb)} | {r.n_spherules:,} | "
                 f"{_n(r.exact_stable)} | {_n(r.exact_first)} | " + " | ".join(cells)
                 + f" | **{r.verdict.tier}** |")
    L += ["", f"Cost basis: {cfg.zaps_per_hour:g} zaps per hour at ${cfg.usd_per_hour:g} per hour. "
          "% of slide, hours and cost are for the stable n.", ""]

    L += ["## Feasibility verdicts", "",
          "**Placeholder limits pending a real budget/instrument allocation from the team:** "
          f"Feasible up to {cfg.feasible_max_hours:g} machine hours; Costly up to "
          f"{cfg.infeasible_min_hours:g} hours; Infeasible above {cfg.infeasible_min_hours:g} hours "
          f"or above {cfg.infeasible_slide_fraction:.0%} of the slide; Impossible when no number "
          "of zaps on the slide is enough.", ""]
    for r in results:
        L.append(f"- **{r.label}: {r.verdict.tier}** under {conditions(cfg)}: {r.verdict.reason}.")
    L += ["", "Notes:", "",
          "- Targeted sampling (for example Fe-rich areas only, planned for v2) could reduce the "
          "required zaps, so an Infeasible verdict here does not rule out targeted strategies.",
          "- Required zaps scale roughly with spherule grade ÷ bulk concentration: richer, rarer "
          "spherules (for example 1 ppm) would need roughly proportionally more zaps.", ""]
    if not cfg.hits_detectable:
        L += [f"**Warning:** the hit threshold ({cfg.hit_threshold_ppb:g} ppb) is above the "
              f"spherule grade ({cfg.grade_ppb:g} ppb), so no zap can count as a hit.", ""]

    L += ["## Simulation cross-check", "",
          f"Overall: **{'PASS' if all_passed else 'FAIL'}**. {cfg.runs:,} runs per concentration, "
          f"seed {cfg.seed}. The simulated success fraction must stay within {cfg.k_se:g} standard "
          "errors of the exact probability at every checked n (log-spaced n across the slide plus "
          f"every n within ±{cfg.window_fraction:.0%} of the exact first crossing and stable n).", "",
          "| Nominal [Ir] | Check | n values checked | Largest gap (SE) | at n | "
          "Largest gap (points) | at n | Simulated stable n* | Simulated first crossing* |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        k = r.check
        L.append(f"| {r.label} | {'PASS' if k.passed else 'FAIL'} | {k.n_checked:,} | "
                 f"{k.max_gap_se:.2f} | {k.max_gap_se_at:,} | {k.max_gap_pp:.2f} | "
                 f"{k.max_gap_pp_at:,} | {_n(r.sim_stable)} | {_n(r.sim_first)} |")
    L += ["", "\\* Simulated, noisy, informational only. These never decide pass or fail.", ""]

    L += ["## Figures", ""]
    for r in results:
        L += [f"### {r.label}", "", f"![{r.label}]({r.figure})", ""]

    L += ["## Parameters", "", "Full record in `params.json`.", "", "```json",
          json.dumps(cfg.as_dict(), indent=2), "```", ""]
    path.write_text("\n".join(L))
