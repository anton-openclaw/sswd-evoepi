"""Tests for sswd_evoepi.results and sswd_evoepi.metrics modules."""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from sswd_evoepi.results import (
    REGION_ORDER,
    SCORED_REGIONS,
    SeedResult,
    RunResult,
    load_seed,
    load_run,
    load_sweep,
)
from sswd_evoepi.metrics import (
    RECOVERY_TARGETS,
    ARRIVAL_TARGETS,
    rmse_log,
    score_regions,
    score_run,
    score_seed,
    compare_runs,
    format_comparison_table,
    quick_score,
)

# ---------------------------------------------------------------------------
# Fake result JSON factory
# ---------------------------------------------------------------------------

def _make_fake_result(seed: int = 42, rmse: float = 0.663) -> dict:
    """Create a minimal fake result_seed JSON dict."""
    return {
        "seed": seed,
        "wall_time_seconds": 100.0,
        "scoring": {
            "rmse_log": rmse,
            "per_region": {},
            "within_2x": 3,
            "within_5x": 6,
            "n_targets": 8,
        },
        "region_recovery": {
            "AK-PWS": 0.031, "AK-FN": 0.048, "AK-FS": 0.105,
            "BC-N": 0.152, "SS-S": 0.183, "JDF": 0.201,
            "OR": 0.035, "CA-N": 0.013,
            "AK-AL": 0.0, "AK-WG": 0.0, "AK-OC": 0.001,
            "AK-EG": 0.002, "BC-C": 0.10, "SS-N": 0.12,
            "WA-O": 0.08, "CA-C": 0.005, "CA-S": 0.001, "BJ": 0.0,
        },
        "region_details": {
            "AK-PWS": {
                "n_nodes": 50, "peak_pop": 10000, "final_pop": 310,
                "recovery_frac": 0.031, "crash_pct": 96.9,
                "yearly_totals": [10000, 9000, 5000, 310],
                "yearly_recruits": [0, 100, 50, 20],
                "yearly_disease_deaths": [0, 3000, 2000, 100],
                "yearly_mean_resistance": [0.15, 0.16, 0.17, 0.18],
                "yearly_mean_tolerance": [0.10, 0.10, 0.11, 0.11],
                "yearly_mean_recovery": [0.02, 0.02, 0.03, 0.03],
                "yearly_va_resistance": [0.01, 0.01, 0.009, 0.008],
                "yearly_va_tolerance": [0.005, 0.005, 0.004, 0.004],
                "yearly_va_recovery": [0.002, 0.002, 0.002, 0.001],
            }
        },
        "overall": {"pop_crash_pct": 85.0, "final_pop_frac": 0.15},
        "arrival_timing": {},
    }


def _write_result(directory: Path, seed: int, rmse: float = 0.663) -> Path:
    """Write a fake result JSON file and return the path."""
    data = _make_fake_result(seed=seed, rmse=rmse)
    path = directory / f"result_seed{seed}.json"
    path.write_text(json.dumps(data))
    return path


@pytest.fixture
def sweep_dir(tmp_path):
    """Create a fake sweep directory with 2 configs × 2 seeds each."""
    w1 = tmp_path / "W1"
    w1.mkdir()
    _write_result(w1, seed=1, rmse=0.663)
    _write_result(w1, seed=2, rmse=0.700)

    w10 = tmp_path / "W10"
    w10.mkdir()
    _write_result(w10, seed=1, rmse=0.800)
    _write_result(w10, seed=2, rmse=0.850)

    return tmp_path


# ---------------------------------------------------------------------------
# results.py tests
# ---------------------------------------------------------------------------

class TestRegionConstants:
    def test_region_order_length(self):
        assert len(REGION_ORDER) == 18

    def test_scored_regions_length(self):
        assert len(SCORED_REGIONS) == 8

    def test_scored_regions_subset_of_order(self):
        assert all(r in REGION_ORDER for r in SCORED_REGIONS)


class TestSeedResult:
    def test_from_json(self):
        data = _make_fake_result(seed=42)
        sr = SeedResult.from_json(data)
        assert sr.seed == 42
        assert sr.wall_time == 100.0
        assert sr.rmse_log == 0.663
        assert sr.recovery("AK-PWS") == pytest.approx(0.031)

    def test_from_file(self, tmp_path):
        path = _write_result(tmp_path, seed=7)
        sr = SeedResult.from_file(path)
        assert sr.seed == 7

    def test_yearly_pop(self):
        data = _make_fake_result()
        sr = SeedResult.from_json(data)
        pop = sr.yearly_pop("AK-PWS")
        assert pop == [10000, 9000, 5000, 310]

    def test_yearly_pop_missing_region(self):
        sr = SeedResult.from_json(_make_fake_result())
        assert sr.yearly_pop("NONEXISTENT") == []

    def test_yearly_trait(self):
        sr = SeedResult.from_json(_make_fake_result())
        assert sr.yearly_trait("AK-PWS", "resistance") == [0.15, 0.16, 0.17, 0.18]
        assert sr.yearly_trait("AK-PWS", "tolerance") == [0.10, 0.10, 0.11, 0.11]
        assert sr.yearly_trait("AK-PWS", "recovery") == [0.02, 0.02, 0.03, 0.03]

    def test_yearly_trait_invalid(self):
        sr = SeedResult.from_json(_make_fake_result())
        with pytest.raises(ValueError):
            sr.yearly_trait("AK-PWS", "speed")

    def test_raw_preserved(self):
        data = _make_fake_result()
        sr = SeedResult.from_json(data)
        assert sr.raw == data


class TestRunResult:
    def test_from_dir(self, sweep_dir):
        run = RunResult.from_dir(sweep_dir / "W1")
        assert run.config_name == "W1"
        assert run.n_seeds() == 2

    def test_mean_rmse(self, sweep_dir):
        run = load_run(sweep_dir / "W1")
        expected = (0.663 + 0.700) / 2
        assert run.mean_rmse() == pytest.approx(expected)

    def test_mean_recovery(self, sweep_dir):
        run = load_run(sweep_dir / "W1")
        # Both seeds have the same recovery values
        assert run.mean_recovery("AK-PWS") == pytest.approx(0.031)

    def test_recovery_std_single_seed(self, tmp_path):
        d = tmp_path / "single"
        d.mkdir()
        _write_result(d, seed=1)
        run = load_run(d)
        assert run.recovery_std("AK-PWS") == 0.0

    def test_empty_dir(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        run = load_run(d)
        assert run.n_seeds() == 0
        assert run.mean_rmse() == float("inf")
        assert run.mean_recovery("AK-PWS") == 0.0


class TestLoadSweep:
    def test_sweep_loads_all(self, sweep_dir):
        runs = load_sweep(sweep_dir)
        assert len(runs) == 2
        assert "W1" in runs
        assert "W10" in runs

    def test_sweep_sorted_by_rmse(self, sweep_dir):
        runs = load_sweep(sweep_dir)
        names = list(runs.keys())
        # W1 has lower mean RMSE than W10
        assert names[0] == "W1"
        assert names[1] == "W10"

    def test_sweep_specific_configs(self, sweep_dir):
        runs = load_sweep(sweep_dir, configs=["W10"])
        assert len(runs) == 1
        assert "W10" in runs

    def test_sweep_empty(self, tmp_path):
        runs = load_sweep(tmp_path)
        assert len(runs) == 0


class TestLoadSeed:
    def test_convenience(self, tmp_path):
        path = _write_result(tmp_path, seed=99)
        sr = load_seed(path)
        assert sr.seed == 99


# ---------------------------------------------------------------------------
# metrics.py tests
# ---------------------------------------------------------------------------

class TestRecoveryTargets:
    def test_values(self):
        assert RECOVERY_TARGETS["AK-PWS"] == 0.30
        assert RECOVERY_TARGETS["AK-FN"] == 0.30
        assert RECOVERY_TARGETS["AK-FS"] == 0.20
        assert RECOVERY_TARGETS["BC-N"] == 0.20
        assert RECOVERY_TARGETS["SS-S"] == 0.05
        assert RECOVERY_TARGETS["JDF"] == 0.02
        assert RECOVERY_TARGETS["OR"] == 0.0025
        assert RECOVERY_TARGETS["CA-N"] == 0.001

    def test_count(self):
        assert len(RECOVERY_TARGETS) == 8

    def test_all_scored_regions_present(self):
        for r in SCORED_REGIONS:
            assert r in RECOVERY_TARGETS


class TestArrivalTargets:
    def test_count(self):
        assert len(ARRIVAL_TARGETS) == 17

    def test_origin(self):
        assert ARRIVAL_TARGETS["CA-S"] == 0


class TestRmseLog:
    def test_perfect_match(self):
        actual = dict(RECOVERY_TARGETS)
        assert rmse_log(actual) == pytest.approx(0.0, abs=1e-10)

    def test_positive_rmse(self):
        actual = {r: 0.001 for r in RECOVERY_TARGETS}
        result = rmse_log(actual)
        assert result > 0
        assert math.isfinite(result)

    def test_zero_actual_uses_eps(self):
        actual = {r: 0.0 for r in RECOVERY_TARGETS}
        result = rmse_log(actual)
        assert math.isfinite(result)
        assert result > 0

    def test_custom_targets(self):
        targets = {"A": 0.5, "B": 0.5}
        actual = {"A": 0.5, "B": 0.5}
        assert rmse_log(actual, targets) == pytest.approx(0.0, abs=1e-10)

    def test_empty_targets(self):
        assert rmse_log({}, {}) == float("inf")


class TestScoreRegions:
    def test_perfect_match(self):
        actual = dict(RECOVERY_TARGETS)
        scored = score_regions(actual)
        for region, info in scored.items():
            assert info["within_2x"] is True
            assert info["grade_str"] == "✓ 2×"

    def test_within_5x(self):
        # AK-PWS target is 0.30; 0.08 gives fold ~3.75 → within 5× but not 2×
        actual = dict(RECOVERY_TARGETS)
        actual["AK-PWS"] = 0.08
        scored = score_regions(actual)
        assert scored["AK-PWS"]["within_2x"] is False
        assert scored["AK-PWS"]["within_5x"] is True
        assert scored["AK-PWS"]["grade_str"] == "~ 5×"

    def test_beyond_5x(self):
        actual = dict(RECOVERY_TARGETS)
        actual["AK-PWS"] = 0.001  # 300× off from 0.30
        scored = score_regions(actual)
        assert scored["AK-PWS"]["within_2x"] is False
        assert scored["AK-PWS"]["within_5x"] is False
        assert "✗" in scored["AK-PWS"]["grade_str"]


class TestScoreRun:
    def test_score_run(self, sweep_dir):
        run = load_run(sweep_dir / "W1")
        scored = score_run(run)
        assert "rmse_log" in scored
        assert "per_region" in scored
        assert "within_2x" in scored
        assert scored["n_targets"] == 8
        assert math.isfinite(scored["rmse_log"])

    def test_score_empty_run(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        run = load_run(d)
        scored = score_run(run)
        assert scored["rmse_log"] > 0


class TestScoreSeed:
    def test_score_seed(self):
        sr = SeedResult.from_json(_make_fake_result())
        scored = score_seed(sr)
        assert "rmse_log" in scored
        assert math.isfinite(scored["rmse_log"])
        assert scored["n_targets"] == 8


class TestCompareRuns:
    def test_compare(self, sweep_dir):
        runs = load_sweep(sweep_dir)
        comparisons = compare_runs(runs)
        assert len(comparisons) == 2
        # Should be sorted by RMSE (lowest first)
        assert comparisons[0]["rmse"] <= comparisons[1]["rmse"]

    def test_compare_by_within_2x(self, sweep_dir):
        runs = load_sweep(sweep_dir)
        comparisons = compare_runs(runs, sort_by="within_2x")
        assert len(comparisons) == 2


class TestFormatComparisonTable:
    def test_format(self, sweep_dir):
        runs = load_sweep(sweep_dir)
        table = format_comparison_table(runs)
        assert "Config" in table
        assert "RMSE" in table
        assert "W1" in table
        assert "W10" in table

    def test_empty(self):
        table = format_comparison_table({})
        assert "no runs" in table


class TestQuickScore:
    def test_quick_score(self, sweep_dir):
        result = quick_score(sweep_dir / "W1")
        assert "W1" in result
        assert "RMSE=" in result
        assert "AK-PWS" in result
        assert "within 2×" in result
