"""Tests for sswd_evoepi.config — configuration loading and validation."""

import tempfile
from pathlib import Path

import pytest
import yaml

from sswd_evoepi.config import (
    ConservationSection,
    DiseaseSection,
    GeneticsSection,
    OutputSection,
    PopulationSection,
    SimulationConfig,
    SimulationSection,
    SpatialSection,
    deep_merge,
    default_config,
    load_config,
    validate_config,
)


# ── deep_merge tests ──────────────────────────────────────────────────

class TestDeepMerge:
    def test_simple_override(self):
        base = {'a': 1, 'b': 2}
        override = {'b': 3}
        result = deep_merge(base, override)
        assert result == {'a': 1, 'b': 3}

    def test_nested_merge(self):
        base = {'x': {'a': 1, 'b': 2}, 'y': 10}
        override = {'x': {'b': 3, 'c': 4}}
        result = deep_merge(base, override)
        assert result == {'x': {'a': 1, 'b': 3, 'c': 4}, 'y': 10}

    def test_new_key(self):
        base = {'a': 1}
        override = {'b': 2}
        result = deep_merge(base, override)
        assert result == {'a': 1, 'b': 2}

    def test_override_dict_with_scalar(self):
        base = {'a': {'nested': 1}}
        override = {'a': 'replaced'}
        result = deep_merge(base, override)
        assert result == {'a': 'replaced'}

    def test_empty_override(self):
        base = {'a': 1, 'b': 2}
        result = deep_merge(base, {})
        assert result == {'a': 1, 'b': 2}


# ── default_config tests ─────────────────────────────────────────────

class TestDefaultConfig:
    def test_creates_valid_config(self):
        config = default_config()
        assert isinstance(config, SimulationConfig)

    def test_default_values(self):
        config = default_config()
        assert config.simulation.start_year == 1910
        assert config.simulation.end_year == 2100
        assert config.simulation.seed == 42
        assert config.disease.scenario == "ubiquitous"
        assert config.genetics.n_loci == 52
        assert config.genetics.n_additive == 51

    def test_no_cost_of_resistance(self):
        """CE-1: cost_resistance must not appear in config."""
        config = default_config()
        assert not hasattr(config.genetics, 'cost_resistance')

    def test_both_scenarios_configurable(self):
        """CE-2: disease.scenario accepts both values."""
        config = default_config()
        assert config.disease.scenario in {"ubiquitous", "invasion"}

    def test_annual_survival_length(self):
        config = default_config()
        assert len(config.population.annual_survival) == 5


# ── YAML loading tests ───────────────────────────────────────────────

class TestLoadConfig:
    def test_load_from_yaml(self, tmp_path):
        """Load a minimal YAML config."""
        yaml_content = {
            'simulation': {'seed': 99, 'start_year': 1900, 'end_year': 2050},
            'disease': {'scenario': 'ubiquitous'},
        }
        config_path = tmp_path / "test.yaml"
        with open(config_path, 'w') as f:
            yaml.dump(yaml_content, f)

        config = load_config(config_path)
        assert config.simulation.seed == 99
        assert config.simulation.start_year == 1900
        # Unspecified sections get defaults
        assert config.genetics.n_loci == 52
        assert config.population.L_inf == 1000.0

    def test_load_with_scenario_override(self, tmp_path):
        """Scenario YAML overrides base."""
        base = {
            'simulation': {'seed': 42, 'start_year': 1910, 'end_year': 2100},
            'disease': {'scenario': 'ubiquitous', 'K_half': 87000.0},
        }
        scenario = {
            'disease': {'K_half': 50000.0},
        }
        base_path = tmp_path / "base.yaml"
        scen_path = tmp_path / "scenario.yaml"
        with open(base_path, 'w') as f:
            yaml.dump(base, f)
        with open(scen_path, 'w') as f:
            yaml.dump(scenario, f)

        config = load_config(base_path, scenario_path=scen_path)
        assert config.disease.K_half == 50000.0
        assert config.disease.scenario == "ubiquitous"  # unchanged

    def test_load_with_sweep_overrides(self, tmp_path):
        """Dict-based sweep overrides are applied last."""
        base = {
            'simulation': {'seed': 42, 'start_year': 1910, 'end_year': 2100},
        }
        base_path = tmp_path / "base.yaml"
        with open(base_path, 'w') as f:
            yaml.dump(base, f)

        overrides = {'simulation': {'seed': 123}}
        config = load_config(base_path, sweep_overrides=overrides)
        assert config.simulation.seed == 123

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_config(tmp_path / "nonexistent.yaml")

    def test_load_real_default_yaml(self):
        """Load the actual configs/default.yaml from the project."""
        default_path = Path(__file__).parent.parent / "configs" / "default.yaml"
        if default_path.exists():
            config = load_config(default_path)
            assert config.simulation.seed == 42
            assert config.disease.scenario == "ubiquitous"
            assert config.disease.Ea_I2D == 2000.0  # ERRATA E1
            assert config.disease.sigma_1_eff == 5.0  # ERRATA E2
            assert config.disease.sigma_D == 150.0    # ERRATA E14
            assert config.genetics.n_loci == 52


# ── Validation tests ──────────────────────────────────────────────────

class TestValidation:
    def test_invalid_scenario(self):
        config = default_config()
        config.disease.scenario = "invalid"
        with pytest.raises(ValueError, match="disease.scenario"):
            validate_config(config)

    def test_invasion_requires_year(self):
        config = default_config()
        config.disease.scenario = "invasion"
        config.disease.invasion_year = None
        with pytest.raises(ValueError, match="invasion_year"):
            validate_config(config)

    def test_invasion_valid_with_year(self):
        config = default_config()
        config.disease.scenario = "invasion"
        config.disease.invasion_year = 2013
        config.disease.invasion_nodes = [0, 5]
        validate_config(config)  # should not raise

    def test_year_ordering(self):
        config = default_config()
        config.simulation.start_year = 2100
        config.simulation.end_year = 1910
        with pytest.raises(ValueError, match="start_year"):
            validate_config(config)

    def test_epidemic_before_start(self):
        config = default_config()
        config.simulation.epidemic_year = 1800
        with pytest.raises(ValueError, match="epidemic_year"):
            validate_config(config)

    def test_loci_consistency(self):
        config = default_config()
        config.genetics.n_loci = 100  # mismatch with n_additive=51
        with pytest.raises(ValueError, match="n_loci"):
            validate_config(config)

    def test_survival_array_wrong_length(self):
        config = default_config()
        config.population.annual_survival = [0.03, 0.90, 0.95]  # only 3
        with pytest.raises(ValueError, match="annual_survival"):
            validate_config(config)

    def test_ubiquitous_valid(self):
        """Ubiquitous scenario doesn't need invasion_year."""
        config = default_config()
        config.disease.scenario = "ubiquitous"
        validate_config(config)  # should not raise


# ── Section dataclass tests ───────────────────────────────────────────

class TestSections:
    def test_disease_section_defaults(self):
        ds = DiseaseSection()
        assert ds.scenario == "ubiquitous"
        assert ds.Ea_I2D == 2000.0  # ERRATA E1
        assert ds.mu_I2D_ref == 0.173
        assert ds.sigma_D == 150.0   # ERRATA E14
        assert ds.T_ref == 20.0

    def test_genetics_section_no_cost(self):
        """CE-1: No cost_resistance attribute."""
        gs = GeneticsSection()
        assert not hasattr(gs, 'cost_resistance')
        assert gs.n_loci == 52
        assert gs.s_het == 0.19

    def test_population_section_defaults(self):
        ps = PopulationSection()
        assert ps.alpha_srs == 1.35
        assert ps.gamma_fert == 4.5
        assert ps.L_inf == 1000.0
        assert len(ps.annual_survival) == 5
