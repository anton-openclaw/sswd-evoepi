"""Tests for sswd_evoepi.types — core data types, enums, and agent arrays."""

import numpy as np
import pytest

from sswd_evoepi.types import (
    AGENT_DTYPE,
    ANNUAL_SURVIVAL,
    N_ADDITIVE,
    N_LOCI,
    IDX_EF1A,
    STAGE_SIZE_THRESHOLDS,
    DiseaseState,
    LarvalCohort,
    NodeSnapshot,
    Origin,
    SettlerPacket,
    Stage,
    Tier,
    allocate_agents,
    allocate_genotypes,
)


# ── Enum tests ────────────────────────────────────────────────────────

class TestStageEnum:
    def test_values(self):
        assert Stage.EGG_LARVA == 0
        assert Stage.SETTLER == 1
        assert Stage.JUVENILE == 2
        assert Stage.SUBADULT == 3
        assert Stage.ADULT == 4

    def test_count(self):
        assert len(Stage) == 5

    def test_integer_compatible(self):
        """Stages can be used as numpy int8 array indices."""
        arr = np.zeros(5, dtype=np.float64)
        arr[Stage.ADULT] = 1.0
        assert arr[4] == 1.0


class TestDiseaseStateEnum:
    def test_values(self):
        assert DiseaseState.S == 0
        assert DiseaseState.E == 1
        assert DiseaseState.I1 == 2
        assert DiseaseState.I2 == 3
        assert DiseaseState.D == 4
        assert DiseaseState.R == 5

    def test_count(self):
        assert len(DiseaseState) == 6


class TestOriginEnum:
    def test_values(self):
        assert Origin.WILD == 0
        assert Origin.CAPTIVE_BRED == 1
        assert Origin.AGF_CROSS == 2
        assert Origin.WILD_SOURCE_REL == 3

    def test_count(self):
        assert len(Origin) == 4


class TestTierEnum:
    def test_values(self):
        assert Tier.TIER1 == 1
        assert Tier.TIER2 == 2


# ── Constants tests ───────────────────────────────────────────────────

class TestConstants:
    def test_loci_counts(self):
        assert N_LOCI == 52
        assert N_ADDITIVE == 51
        assert IDX_EF1A == 51
        assert N_LOCI == N_ADDITIVE + 1

    def test_stage_thresholds(self):
        assert STAGE_SIZE_THRESHOLDS[Stage.SETTLER] == 10.0
        assert STAGE_SIZE_THRESHOLDS[Stage.JUVENILE] == 150.0
        assert STAGE_SIZE_THRESHOLDS[Stage.SUBADULT] == 400.0
        # EGG_LARVA and ADULT have no size threshold (transitions are handled differently)
        assert Stage.EGG_LARVA not in STAGE_SIZE_THRESHOLDS
        assert Stage.ADULT not in STAGE_SIZE_THRESHOLDS

    def test_annual_survival(self):
        assert len(ANNUAL_SURVIVAL) == 5
        assert ANNUAL_SURVIVAL.dtype == np.float64
        # All survival values in [0, 1]
        assert np.all(ANNUAL_SURVIVAL >= 0.0)
        assert np.all(ANNUAL_SURVIVAL <= 1.0)
        # Specific values from spec
        assert ANNUAL_SURVIVAL[Stage.SETTLER] == pytest.approx(0.03)
        assert ANNUAL_SURVIVAL[Stage.ADULT] == pytest.approx(0.98)


# ── AGENT_DTYPE tests ─────────────────────────────────────────────────

class TestAgentDtype:
    def test_has_all_fields(self):
        """AGENT_DTYPE has every field from the integration spec §4.1 + spawning fields."""
        expected_fields = [
            'x', 'y', 'heading', 'speed',       # Spatial
            'size', 'age', 'stage', 'sex',       # Life history
            'disease_state', 'disease_timer',     # Disease
            'resistance', 'fecundity_mod',        # Genetics
            'spawning_ready', 'has_spawned', 'spawn_refractory',  # Spawning (Phase 1)
            'spawn_gravity_timer', 'immunosuppression_timer', 'last_spawn_day',  # Spawning (continued)
            'node_id', 'alive', 'origin', 'cause_of_death',  # Administrative
            'pathogen_virulence',  # Pathogen evolution
            'settlement_day',  # Juvenile immunity (Phase 11)
        ]
        actual_fields = [name for name in AGENT_DTYPE.names]
        assert actual_fields == expected_fields

    def test_field_dtypes(self):
        """Field types match the spec."""
        assert AGENT_DTYPE['x'] == np.float32
        assert AGENT_DTYPE['y'] == np.float32
        assert AGENT_DTYPE['size'] == np.float32
        assert AGENT_DTYPE['age'] == np.float32
        assert AGENT_DTYPE['stage'] == np.int8
        assert AGENT_DTYPE['sex'] == np.int8
        assert AGENT_DTYPE['disease_state'] == np.int8
        assert AGENT_DTYPE['disease_timer'] == np.int16
        assert AGENT_DTYPE['resistance'] == np.float32
        assert AGENT_DTYPE['fecundity_mod'] == np.float32
        assert AGENT_DTYPE['node_id'] == np.int16
        assert AGENT_DTYPE['alive'] == np.bool_
        assert AGENT_DTYPE['origin'] == np.int8

    def test_origin_field_present(self):
        """ERRATA E13: origin field must exist."""
        assert 'origin' in AGENT_DTYPE.names

    def test_disease_timer_is_int16(self):
        """ERRATA E4: disease_timer needs countdown range."""
        assert AGENT_DTYPE['disease_timer'] == np.int16
        # int16 range: -32768 to 32767 — enough for countdown days


# ── Allocation tests ──────────────────────────────────────────────────

class TestAllocateAgents:
    def test_shape_and_dtype(self):
        agents = allocate_agents(100)
        assert agents.shape == (100,)
        assert agents.dtype == AGENT_DTYPE

    def test_zeroed(self):
        agents = allocate_agents(50)
        assert np.all(agents['alive'] == False)
        assert np.all(agents['size'] == 0.0)
        assert np.all(agents['resistance'] == 0.0)

    def test_zero_size(self):
        agents = allocate_agents(0)
        assert agents.shape == (0,)

    def test_field_writeable(self):
        agents = allocate_agents(10)
        agents[0]['alive'] = True
        agents[0]['stage'] = Stage.ADULT
        agents[0]['disease_state'] = DiseaseState.S
        agents[0]['origin'] = Origin.CAPTIVE_BRED
        agents[0]['resistance'] = 0.42
        agents[0]['fecundity_mod'] = 1.0  # CE-1: always 1.0
        assert agents[0]['alive'] == True
        assert agents[0]['stage'] == Stage.ADULT
        assert agents[0]['origin'] == Origin.CAPTIVE_BRED
        assert agents[0]['resistance'] == pytest.approx(0.42, abs=1e-6)


class TestAllocateGenotypes:
    def test_shape_and_dtype(self):
        geno = allocate_genotypes(100)
        assert geno.shape == (100, N_LOCI, 2)
        assert geno.dtype == np.int8

    def test_zeroed(self):
        geno = allocate_genotypes(50)
        assert np.all(geno == 0)

    def test_diploid_encoding(self):
        """Can set alleles per individual per locus."""
        geno = allocate_genotypes(5)
        # Individual 0, locus 3, allele copy 1 → resistant
        geno[0, 3, 1] = 1
        assert geno[0, 3, 0] == 0  # other copy unchanged
        assert geno[0, 3, 1] == 1

    def test_ef1a_locus(self):
        """EF1A locus is at the expected index."""
        geno = allocate_genotypes(3)
        # Heterozygote at EF1A
        geno[0, IDX_EF1A, 0] = 0
        geno[0, IDX_EF1A, 1] = 1
        assert geno[0, IDX_EF1A].sum() == 1  # heterozygous
        # Homozygous resistant (lethal — flagged by genetics module, not here)
        geno[1, IDX_EF1A, :] = 1
        assert geno[1, IDX_EF1A].sum() == 2


# ── Data transfer object tests ────────────────────────────────────────

class TestLarvalCohort:
    def test_creation(self):
        cohort = LarvalCohort(
            source_node=5,
            n_competent=100,
            genotypes=np.zeros((100, N_LOCI, 2), dtype=np.int8),
            parent_pairs=np.zeros((100, 2), dtype=np.int32),
            pld_days=63.0,
        )
        assert cohort.source_node == 5
        assert cohort.n_competent == 100
        assert cohort.genotypes.shape == (100, N_LOCI, 2)
        assert cohort.pld_days == 63.0
        # Default values for new fields
        assert cohort.spawn_day == 0
        assert cohort.sst_at_spawn == 10.5

    def test_creation_with_spawn_day(self):
        cohort = LarvalCohort(
            source_node=3,
            n_competent=50,
            genotypes=np.zeros((50, N_LOCI, 2), dtype=np.int8),
            parent_pairs=np.zeros((50, 2), dtype=np.int32),
            pld_days=53.0,
            spawn_day=120,
            sst_at_spawn=14.0,
        )
        assert cohort.spawn_day == 120
        assert cohort.sst_at_spawn == 14.0
        assert cohort.pld_days == 53.0


class TestSettlerPacket:
    def test_creation(self):
        packet = SettlerPacket(
            dest_node=12,
            genotypes=np.zeros((30, N_LOCI, 2), dtype=np.int8),
            source_nodes=np.zeros(30, dtype=np.int16),
        )
        assert packet.dest_node == 12
        assert packet.genotypes.shape == (30, N_LOCI, 2)


class TestNodeSnapshot:
    def test_creation(self):
        snap = NodeSnapshot(
            node_id=0, tier=1, n_alive=450, n_adults=200,
            mean_resistance=0.15, var_resistance=0.008,
            vibrio_concentration=1200.0, sst=14.5,
            salinity=30.0, flushing_rate=0.1,
            R0_estimate=1.8, epidemic_active=True,
        )
        assert snap.node_id == 0
        assert snap.epidemic_active is True
        assert snap.mean_resistance == pytest.approx(0.15)
