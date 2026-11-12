"""Command line entry point: run every concentration and write the outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from . import exact, model, plots, report, simulate, validate
from .config import Config, concentration_label, load_config
from .cost import verdict_for


def run(cfg: Config, out: Path, config_path: str = "") -> list[report.Result]:
    out.mkdir(parents=True, exist_ok=True)
    N = cfg.total_cells
    band = cfg.band
    n_all = np.arange(1, N + 1)
    results = []
    seeds = np.random.SeedSequence(cfg.seed).spawn(len(cfg.concentrations_ppb))

    for conc, seed in zip(cfg.concentrations_ppb, seeds):
        label = concentration_label(conc)
        K = cfg.n_spherules(conc)
        realized = cfg.realized_ppb(conc)
        sim_seed, traj_seed = seed.spawn(2)

        p_exact = exact.success_probability(N, K, band, cfg.hits_detectable)
        first = exact.first_crossing(p_exact, cfg.confidence)
        stable = exact.stable_n(p_exact, cfg.confidence)

        p_sim = simulate.simulate_success_fraction(
            N, K, band, cfg.runs, np.random.default_rng(sim_seed), cfg.hits_detectable)
        points = validate.check_points(N, [first, stable], cfg.log_points, cfg.window_fraction)
        check = validate.overlay_check(p_sim, p_exact, cfg.runs, points, cfg.k_se)

        traj_rng = np.random.default_rng(traj_seed)
        trajectories = []
        for _ in range(cfg.trajectories if K > 0 else 0):
            hits = simulate.reference_run(N, K, cfg.grade_ppb, cfg.hit_threshold_ppb, traj_rng)
            trajectories.append(model.estimate(hits, n_all, cfg.grade_ppb) / realized)

        figure = f"success_{label.replace(' ', '')}.png"
        plots.plot_concentration(
            out / figure, label=label, realized_label=f"{realized:.6g} ppb", n_cells=N,
            trajectories=trajectories, sim=p_sim, exact=p_exact,
            band=(cfg.band_low, cfg.band_high), confidence=cfg.confidence, runs=cfg.runs,
            first=first, stable=stable)

        results.append(report.Result(
            nominal_ppb=conc, label=label, realized_ppb=realized, n_spherules=K,
            exact_first=first, exact_stable=stable,
            sim_first=exact.first_crossing(p_sim, cfg.confidence),
            sim_stable=exact.stable_n(p_sim, cfg.confidence),
            verdict=verdict_for(stable, cfg), check=check, figure=figure))

    report.write_csv(out / "summary.csv", results, cfg)
    report.write_report(out / "report.md", results, cfg)
    report.write_params(out / "params.json", cfg, config_path)
    return results


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/v1_default.toml", help="TOML parameter file")
    ap.add_argument("--out", default="outputs/v1", help="output folder")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    results = run(cfg, Path(args.out), args.config)

    print(report.best_case_statement(cfg))
    for r in results:
        c = r.cost(cfg)
        detail = "" if c is None else (f"{c.zaps:,} zaps (first crossing {r.exact_first:,}), "
                                       f"{c.slide_fraction:.1%} of slide, {c.hours:,.0f} h, "
                                       f"${c.usd:,.0f}; ")
        print(f"  {r.label:>8}: {detail}{r.verdict.tier}; "
              f"sim check {'PASS' if r.check.passed else 'FAIL'} "
              f"(max gap {r.check.max_gap_se:.2f} SE)")
    print(f"Outputs written to {args.out}/")
    if not all(r.check.passed for r in results):
        print("ERROR: simulation does not agree with the exact calculation.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
