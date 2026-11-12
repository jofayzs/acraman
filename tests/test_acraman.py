from fractions import Fraction

import numpy as np
import pytest

from acraman import exact, model, simulate, validate
from acraman.config import Config, concentration_label, load_config
from acraman.cost import cost_of, verdict_for

CFG = Config()
N = CFG.total_cells
BAND = CFG.band


def test_grid_and_spherule_counts():
    assert N == 212_207
    assert [CFG.n_spherules(c) for c in (1.0, 0.5, 0.001)] == [2122, 1061, 2]
    assert CFG.realized_ppb(0.001) == pytest.approx(0.000942, rel=1e-3)
    assert BAND == (Fraction(4, 5), Fraction(6, 5))


@pytest.mark.parametrize("conc, first, stable", [
    (1.0, 6_334, 6_667),
    (0.5, 12_168, 12_834),
    (0.001, 201_318, 201_318),
])
def test_exact_validation_table(conc, first, stable):
    p = exact.success_probability(N, CFG.n_spherules(conc), BAND)
    assert exact.first_crossing(p, 0.9) == first
    assert exact.stable_n(p, 0.9) == stable


def test_in_band_matches_agreed_integer_rule_and_includes_edges():
    K = 2122
    for h, n in [(60, 6000), (72, 6000), (48, 6000), (80, 6667), (0, 10), (5, 5)]:
        assert model.in_band(h, n, K, N, BAND) == (4 * K * n <= 5 * h * N <= 6 * K * n)
    # exact ties on both edges count as success
    assert model.in_band(4, 100, 5, 100, BAND)       # ratio exactly 0.8
    assert model.in_band(6, 100, 5, 100, BAND)       # ratio exactly 1.2
    assert not model.in_band(7, 100, 5, 100, BAND)
    assert not model.in_band(3, 100, 0, 100, BAND)   # no spherules: never in band


def test_hit_bounds_agree_with_in_band():
    K, n_cells = 37, 1500
    n = np.arange(1, n_cells + 1)
    kmin, kmax = model.hit_bounds(n, K, n_cells, BAND)
    for i in range(0, n_cells, 7):
        ok = [h for h in range(K + 1) if model.in_band(h, int(n[i]), K, n_cells, BAND)]
        expected = [h for h in range(int(kmin[i]), int(kmax[i]) + 1) if h <= K]
        assert ok == expected


def test_fast_and_reference_paths_identical_on_shared_input():
    rng = np.random.default_rng(1)
    n_cells, K = 20_000, 150
    slide = model.build_slide(n_cells, K, 100.0, rng)
    order = model.zap_order(n_cells, rng)
    ref = simulate.reference_hits(slide, order, 100.0)
    positions = simulate.hit_positions(slide, order, 100.0)
    assert positions.size == K
    fast = simulate.hits_from_positions(positions, n_cells)
    assert np.array_equal(ref, fast)

    # success-interval bookkeeping agrees with testing every n directly
    counts = simulate.success_counts(positions[None, :], n_cells, K, BAND)
    direct = [model.in_band(int(h), n, K, n_cells, BAND) for n, h in enumerate(ref, start=1)]
    assert np.array_equal(counts.astype(bool), np.array(direct))


def test_reference_path_matches_exact_on_small_slide():
    n_cells, K, runs, k_se = 2_000, 40, 3_000, 5.0
    rng = np.random.default_rng(7)
    n = np.arange(1, n_cells + 1)
    kmin, kmax = model.hit_bounds(n, K, n_cells, BAND)
    counts = np.zeros(n_cells)
    for _ in range(runs):
        hits = simulate.reference_run(n_cells, K, 100.0, 100.0, rng)
        counts += (hits >= kmin) & (hits <= kmax)
    p = exact.success_probability(n_cells, K, BAND)
    check = validate.overlay_check(counts / runs, p, runs, n, k_se)
    assert check.passed, check


def test_fast_simulation_overlays_exact_on_default_seed():
    K, runs = CFG.n_spherules(1.0), 2_000
    rng = np.random.default_rng(np.random.SeedSequence(CFG.seed))
    p = exact.success_probability(N, K, BAND)
    sim = simulate.simulate_success_fraction(N, K, BAND, runs, rng)
    pts = validate.check_points(N, [6_334, 6_667], CFG.log_points, CFG.window_fraction)
    assert validate.overlay_check(sim, p, runs, pts, CFG.k_se).passed


def test_overlay_check_detects_a_wrong_curve():
    p = np.full(1000, 0.5)
    n = np.arange(1, 1001)
    assert validate.overlay_check(p.copy(), p, 10_000, n, 5.0).passed
    bad = validate.overlay_check(p + 0.05, p, 10_000, n, 5.0)
    assert not bad.passed and bad.max_gap_pp == pytest.approx(5.0)


def test_stable_and_first_crossing_edge_cases():
    assert exact.stable_n(np.array([0.95, 0.5, 0.95, 0.95]), 0.9) == 3
    assert exact.first_crossing(np.array([0.95, 0.5, 0.95, 0.95]), 0.9) == 1
    assert exact.stable_n(np.array([0.95, 0.95]), 0.9) == 1
    assert exact.stable_n(np.array([0.95, 0.5]), 0.9) is None
    assert exact.first_crossing(np.zeros(5), 0.9) is None


def test_cost_arithmetic():
    c = cost_of(6_667, CFG)
    assert c.hours == pytest.approx(444.47, abs=0.01)
    assert c.usd == pytest.approx(44_446.67, abs=0.01)
    assert c.slide_fraction == pytest.approx(0.0314, abs=1e-4)


def test_verdict_tiers_and_triggers():
    assert verdict_for(1_500, CFG).tier == "Feasible"          # exactly 100 h
    assert verdict_for(1_501, CFG).tier == "Costly"
    assert verdict_for(6_667, CFG).tier == "Costly"
    assert verdict_for(7_500, CFG).tier == "Costly"            # exactly 500 h
    v = verdict_for(12_834, CFG)
    assert v.tier == "Infeasible" and "hours" in v.reason and "slide" not in v.reason
    v = verdict_for(201_318, CFG)
    assert v.tier == "Infeasible" and "hours" in v.reason and "slide" in v.reason
    assert verdict_for(None, CFG).tier == "Impossible"
    # slide-area trigger alone
    cheap = Config(zaps_per_hour=1e6)
    v = verdict_for(150_000, cheap)
    assert v.tier == "Infeasible" and "slide" in v.reason and "hours" not in v.reason


def test_zero_spherules_is_impossible():
    assert CFG.n_spherules(0.0002) == 0
    p = exact.success_probability(N, 0, BAND)
    assert exact.stable_n(p, 0.9) is None
    assert verdict_for(exact.stable_n(p, 0.9), CFG).tier == "Impossible"
    assert not simulate.simulate_success_fraction(N, 0, BAND, 10, np.random.default_rng(0)).any()


def test_threshold_above_grade_warns_and_never_succeeds():
    with pytest.warns(UserWarning, match="no zap can ever count"):
        cfg = Config(hit_threshold_ppb=150.0)
    p = exact.success_probability(N, 2122, BAND, cfg.hits_detectable)
    assert exact.stable_n(p, 0.9) is None


def test_default_config_file_matches_code_defaults():
    assert load_config("config/v1_default.toml") == CFG
    with pytest.raises(ValueError):
        Config(band_low=1.1)


def test_labels():
    assert [concentration_label(c) for c in (1.0, 0.5, 0.001)] == ["1 ppb", "0.5 ppb", "1 ppt"]
