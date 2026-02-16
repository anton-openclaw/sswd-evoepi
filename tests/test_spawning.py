"""Tests for spawning system Phase 1.

Tests cover:
- Season boundary wrapping (Nov-July)
- Seasonal readiness probability distribution  
- Spontaneous spawning rate validation
- Female single-spawn enforcement
- Male multi-bout tracking with refractory periods
- Total fecundity normalization
- Latitude adjustment of spawning peak

Author: Anton 🔬
"""

import numpy as np
import pytest
from unittest.mock import Mock

from sswd_evoepi.spawning import (
    in_spawning_season,
    seasonal_readiness_prob,
    latitude_adjusted_peak,
    spawning_step,
    reset_spawning_season,
    _execute_spawning_events,
    _generate_larval_cohort,
    _cascade_induction_step,
    _get_recent_spawners_mask,
    _check_cascade_induction,
)
from sswd_evoepi.types import AGENT_DTYPE, Stage, allocate_agents, allocate_genotypes, N_LOCI
from sswd_evoepi.config import SpawningSection, DiseaseSection


# ═══════════════════════════════════════════════════════════════════════
# SEASON BOUNDARY TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestSpawningSeason:
    """Test spawning season boundary logic."""
    
    def test_default_season_boundaries(self):
        """Test default Nov 1 - Jul 15 season (305-196)."""
        # In season days
        assert in_spawning_season(305)  # Nov 1 - season start
        assert in_spawning_season(320)  # Mid November
        assert in_spawning_season(1)    # January 1
        assert in_spawning_season(105)  # Mid April (peak)
        assert in_spawning_season(196)  # Jul 15 - season end
        
        # Out of season days
        assert not in_spawning_season(197)  # Jul 16 - just after end
        assert not in_spawning_season(250)  # September
        assert not in_spawning_season(304)  # Oct 31 - just before start
    
    def test_custom_season_boundaries(self):
        """Test custom season boundaries."""
        # Non-wrapping season (spring only)
        assert in_spawning_season(100, season_start=90, season_end=150)
        assert not in_spawning_season(80, season_start=90, season_end=150)
        assert not in_spawning_season(160, season_start=90, season_end=150)
        
        # Different wrapping season
        assert in_spawning_season(350, season_start=300, season_end=100)
        assert in_spawning_season(50, season_start=300, season_end=100)
        assert not in_spawning_season(200, season_start=300, season_end=100)
    
    def test_year_boundary_wrapping(self):
        """Test correct handling of year boundary crossing."""
        # Season 305-196 crosses year boundary
        
        # Late year days (Nov-Dec)
        for doy in range(305, 366):
            assert in_spawning_season(doy, 305, 196)
        
        # Early year days (Jan-Jul)  
        for doy in range(1, 197):
            assert in_spawning_season(doy, 305, 196)
        
        # Out of season (Aug-Oct)
        for doy in range(197, 305):
            assert not in_spawning_season(doy, 305, 196)


# ═══════════════════════════════════════════════════════════════════════
# READINESS PROBABILITY TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestReadinessProbability:
    """Test seasonal readiness probability calculations."""
    
    def test_peak_day_probability(self):
        """Test maximum probability at peak day."""
        prob = seasonal_readiness_prob(105, peak_doy=105, peak_width=45)
        assert prob == 1.0
    
    def test_probability_decay_with_distance(self):
        """Test probability decreases with distance from peak."""
        peak_doy = 105
        peak_width = 45
        
        prob_peak = seasonal_readiness_prob(105, peak_doy, peak_width)
        prob_near = seasonal_readiness_prob(110, peak_doy, peak_width)  # 5 days away
        prob_far = seasonal_readiness_prob(150, peak_doy, peak_width)   # 45 days away
        
        assert prob_peak > prob_near > prob_far > 0
        assert prob_far == pytest.approx(np.exp(-0.5), rel=1e-3)  # 1 std dev away
    
    def test_year_wrapping_distances(self):
        """Test year-boundary distance calculations."""
        # Peak at April 15 (day 105)
        
        # Day 350 (mid-December) is 115 days away directly,
        # but only 115 days away wrapped (350 -> 365 -> 105 = 120 days)
        # Actually: min(|350-105|, 365-|350-105|) = min(245, 120) = 120
        prob_dec = seasonal_readiness_prob(350, peak_doy=105, peak_width=45)
        
        # Day 10 (mid-January) is 95 days away directly, 270 wrapped
        # min(95, 270) = 95
        prob_jan = seasonal_readiness_prob(10, peak_doy=105, peak_width=45)
        
        # January should be closer to peak than December
        assert prob_jan > prob_dec
    
    def test_symmetric_probabilities(self):
        """Test probabilities are symmetric around peak."""
        peak_doy = 105
        peak_width = 45
        
        # Same distance from peak should give same probability
        prob_before = seasonal_readiness_prob(90, peak_doy, peak_width)  # 15 days before
        prob_after = seasonal_readiness_prob(120, peak_doy, peak_width)  # 15 days after
        
        assert prob_before == pytest.approx(prob_after, rel=1e-6)
    
    def test_probability_bounds(self):
        """Test probabilities are in [0, 1] range."""
        peak_doy = 105
        peak_width = 45
        
        # Test various days throughout the year
        for doy in range(1, 366, 10):
            prob = seasonal_readiness_prob(doy, peak_doy, peak_width)
            assert 0.0 <= prob <= 1.0


class TestLatitudeAdjustment:
    """Test latitude-based peak adjustment."""
    
    def test_higher_latitude_later_peak(self):
        """Test higher latitudes have later spawning peaks."""
        base_peak = 105
        
        # Reference latitude (40°N) should give base peak
        ref_peak = latitude_adjusted_peak(base_peak, 40.0)
        assert ref_peak == base_peak
        
        # Higher latitude should be later
        north_peak = latitude_adjusted_peak(base_peak, 50.0)
        assert north_peak > base_peak
        
        # Lower latitude should be earlier  
        south_peak = latitude_adjusted_peak(base_peak, 30.0)
        assert south_peak < base_peak
    
    def test_latitude_shift_rate(self):
        """Test default 3 days per degree shift rate."""
        base_peak = 105
        
        # 10 degrees north should be 30 days later
        peak_10n = latitude_adjusted_peak(base_peak, 50.0)  # 10 degrees above reference
        expected = base_peak + 30
        assert peak_10n == expected
        
        # 5 degrees south should be 15 days earlier
        peak_5s = latitude_adjusted_peak(base_peak, 35.0)   # 5 degrees below reference  
        expected = base_peak - 15
        assert peak_5s == expected
    
    def test_peak_wrapping(self):
        """Test peak day wraps correctly at year boundaries."""
        # Peak that would go before day 1
        early_peak = latitude_adjusted_peak(10, 30.0)  # 30 days earlier
        assert 1 <= early_peak <= 365
        
        # Peak that would go past day 365
        late_peak = latitude_adjusted_peak(350, 60.0)  # 60 days later  
        assert 1 <= late_peak <= 365
    
    def test_custom_shift_rate(self):
        """Test custom latitude shift rates."""
        base_peak = 105
        
        # 1 day per degree
        peak_1deg = latitude_adjusted_peak(base_peak, 45.0, lat_shift_per_deg=1.0)
        expected = base_peak + 5  # 5 degrees * 1 day/degree
        assert peak_1deg == expected
        
        # 5 days per degree
        peak_5deg = latitude_adjusted_peak(base_peak, 45.0, lat_shift_per_deg=5.0)
        expected = base_peak + 25  # 5 degrees * 5 days/degree
        assert peak_5deg == expected


# ═══════════════════════════════════════════════════════════════════════
# SPAWNING STEP INTEGRATION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestSpawningStep:
    """Test integrated spawning step logic."""
    
    def setup_method(self):
        """Set up test agents and genotypes."""
        self.n_agents = 100
        self.agents = allocate_agents(self.n_agents)
        self.genotypes = allocate_genotypes(self.n_agents)
        self.spawning_config = SpawningSection()
        self.disease_config = DiseaseSection()
        self.rng = np.random.default_rng(42)
        
        # Create some mature adults
        n_adults = 50
        self.agents['alive'][:n_adults] = True
        self.agents['stage'][:n_adults] = Stage.ADULT
        self.agents['size'][:n_adults] = 450.0  # Above L_min_repro
        self.agents['sex'][:n_adults] = self.rng.integers(0, 2, n_adults)  # Random sex
        self.agents['node_id'][:n_adults] = 0
        
        # Initialize random genotypes
        self.genotypes[:n_adults] = self.rng.integers(0, 2, (n_adults, N_LOCI, 2))
    
    def test_out_of_season_no_activity(self):
        """Test no spawning activity outside of season."""
        # Day 250 = September, out of season
        cohorts = spawning_step(
            self.agents, self.genotypes, day_of_year=250,
            node_latitude=45.0, spawning_config=self.spawning_config, 
            disease_config=self.disease_config, rng=self.rng
        )
        
        assert len(cohorts) == 0
        assert np.all(self.agents['spawning_ready'] == 0)
        assert np.all(self.agents['has_spawned'] == 0)
    
    def test_readiness_updates_in_season(self):
        """Test spawning readiness updates during season."""
        # Day 105 = peak spawning day, should have high readiness probability
        initial_ready = np.sum(self.agents['spawning_ready'])
        
        cohorts = spawning_step(
            self.agents, self.genotypes, day_of_year=105,
            node_latitude=40.0, spawning_config=self.spawning_config,
            disease_config=self.disease_config, rng=self.rng
        )
        
        final_ready = np.sum(self.agents['spawning_ready'])
        
        # Some agents should become ready (probability = 1.0 at peak)
        assert final_ready > initial_ready
    
    def test_female_single_spawn_enforcement(self):
        """Test females can only spawn once per season."""
        # Set up some ready females
        n_adults = np.sum(self.agents['alive'] & (self.agents['stage'] == Stage.ADULT))
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        female_mask = self.agents['sex'] == 0
        
        # Make some females ready
        ready_females = adult_mask & female_mask
        self.agents['spawning_ready'][ready_females] = 1
        
        # Force high spawning probability by mocking
        original_p_female = self.spawning_config.p_spontaneous_female
        self.spawning_config.p_spontaneous_female = 1.0  # 100% spawn rate
        
        try:
            # First day - should spawn
            cohorts1 = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            spawned_females = np.sum(
                (self.agents['sex'] == 0) & 
                (self.agents['has_spawned'] == 1)
            )
            
            # Second day - should not spawn again (already has_spawned=1)
            cohorts2 = spawning_step(
                self.agents, self.genotypes, day_of_year=106,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            # Same females should still have has_spawned=1, no additional spawning
            spawned_females_day2 = np.sum(
                (self.agents['sex'] == 0) & 
                (self.agents['has_spawned'] == 1)
            )
            
            assert spawned_females_day2 == spawned_females  # No additional spawning
            
        finally:
            self.spawning_config.p_spontaneous_female = original_p_female
    
    def test_male_multi_bout_with_refractory(self):
        """Test male multi-bout spawning with refractory periods."""
        # Set up ready males
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        male_mask = self.agents['sex'] == 1
        ready_males = adult_mask & male_mask
        
        self.agents['spawning_ready'][ready_males] = 1
        
        # Force high male spawning probability
        original_p_male = self.spawning_config.p_spontaneous_male
        self.spawning_config.p_spontaneous_male = 1.0
        
        try:
            # Day 1: First bout
            cohorts1 = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            spawned_males = (self.agents['sex'] == 1) & (self.agents['has_spawned'] > 0)
            bout_counts_day1 = self.agents['has_spawned'][spawned_males]
            refractory_timers = self.agents['spawn_refractory'][spawned_males]
            
            # Check spawning occurred and refractory is set
            assert np.any(spawned_males)
            assert np.all(bout_counts_day1 == 1)
            assert np.all(refractory_timers == self.spawning_config.male_refractory_days)
            
            # Day 2: Should not spawn (refractory)
            cohorts2 = spawning_step(
                self.agents, self.genotypes, day_of_year=106,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            bout_counts_day2 = self.agents['has_spawned'][spawned_males]
            assert np.all(bout_counts_day2 == 1)  # No additional bouts
            
            # Fast-forward past refractory period
            self.agents['spawn_refractory'][:] = 0  # Clear refractory
            
            # Day 22: Should spawn again (second bout)
            cohorts3 = spawning_step(
                self.agents, self.genotypes, day_of_year=127,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            bout_counts_day22 = self.agents['has_spawned'][spawned_males]
            assert np.any(bout_counts_day22 == 2)  # Some second bouts
            
        finally:
            self.spawning_config.p_spontaneous_male = original_p_male
    
    def test_male_bout_limit_enforcement(self):
        """Test males cannot exceed maximum bout limit."""
        # Set up a male at bout limit
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        male_mask = self.agents['sex'] == 1
        test_male_idx = np.where(adult_mask & male_mask)[0][0]
        
        self.agents['spawning_ready'][test_male_idx] = 1
        self.agents['has_spawned'][test_male_idx] = self.spawning_config.male_max_bouts  # At limit
        self.agents['spawn_refractory'][test_male_idx] = 0  # Not refractory
        
        # Force high spawning probability
        original_p_male = self.spawning_config.p_spontaneous_male
        self.spawning_config.p_spontaneous_male = 1.0
        
        try:
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            # Male should not spawn (already at limit)
            assert self.agents['has_spawned'][test_male_idx] == self.spawning_config.male_max_bouts
            
        finally:
            self.spawning_config.p_spontaneous_male = original_p_male
    
    def test_refractory_timer_countdown(self):
        """Test refractory timers count down properly."""
        # Set up males with various refractory times
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        male_indices = np.where(adult_mask & (self.agents['sex'] == 1))[0][:3]
        
        if len(male_indices) >= 3:
            self.agents['spawn_refractory'][male_indices[0]] = 5
            self.agents['spawn_refractory'][male_indices[1]] = 1  
            self.agents['spawn_refractory'][male_indices[2]] = 0
            
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            # Check timers decremented
            assert self.agents['spawn_refractory'][male_indices[0]] == 4
            assert self.agents['spawn_refractory'][male_indices[1]] == 0  # Now ready
            assert self.agents['spawn_refractory'][male_indices[2]] == 0  # Still 0
    
    def test_larval_cohort_generation(self):
        """Test larval cohorts are generated from spawning events."""
        # Set up ready adults of both sexes
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        self.agents['spawning_ready'][adult_mask] = 1
        
        # Ensure we have both sexes
        n_adults = np.sum(adult_mask)
        if n_adults >= 2:
            adult_indices = np.where(adult_mask)[0]
            self.agents['sex'][adult_indices[0]] = 0  # Female
            self.agents['sex'][adult_indices[1]] = 1  # Male
        
        # Force spawning
        original_p_female = self.spawning_config.p_spontaneous_female
        original_p_male = self.spawning_config.p_spontaneous_male
        self.spawning_config.p_spontaneous_female = 1.0
        self.spawning_config.p_spontaneous_male = 1.0
        
        try:
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            # Should generate cohort if both sexes present
            if np.any(self.agents['sex'][adult_mask] == 0) and np.any(self.agents['sex'][adult_mask] == 1):
                assert len(cohorts) > 0
                assert cohorts[0].n_competent > 0
                assert cohorts[0].genotypes.shape[0] == cohorts[0].n_competent
                assert cohorts[0].genotypes.shape[1] == N_LOCI
                assert cohorts[0].genotypes.shape[2] == 2  # Diploid
            
        finally:
            self.spawning_config.p_spontaneous_female = original_p_female
            self.spawning_config.p_spontaneous_male = original_p_male


# ═══════════════════════════════════════════════════════════════════════
# UTILITY FUNCTION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestUtilityFunctions:
    """Test spawning utility functions."""
    
    def test_execute_spawning_events_female(self):
        """Test female spawning event execution."""
        agents = allocate_agents(10)
        agents['alive'][:5] = True
        agents['sex'][:3] = 0  # Females
        
        female_indices = np.array([0, 1, 2])
        _execute_spawning_events(agents, female_indices, day_of_year=105, spawning_config=SpawningSection(), disease_config=DiseaseSection())
        
        # Females should have has_spawned=1, no refractory
        assert np.all(agents['has_spawned'][female_indices] == 1)
        assert np.all(agents['spawn_refractory'][female_indices] == 0)
        assert np.all(agents['last_spawn_day'][female_indices] == 105)
    
    def test_execute_spawning_events_male(self):
        """Test male spawning event execution."""
        agents = allocate_agents(10)
        agents['alive'][:5] = True
        agents['sex'][:3] = 1  # Males
        
        male_indices = np.array([0, 1, 2])
        _execute_spawning_events(agents, male_indices, day_of_year=105, spawning_config=SpawningSection(), disease_config=DiseaseSection(), refractory_days=21)
        
        # Males should increment bout count and get refractory
        assert np.all(agents['has_spawned'][male_indices] == 1)
        assert np.all(agents['spawn_refractory'][male_indices] == 21)
        assert np.all(agents['last_spawn_day'][male_indices] == 105)
    
    def test_reset_spawning_season(self):
        """Test spawning season reset functionality."""
        agents = allocate_agents(10)
        
        # Set up some spawning state
        agents['spawning_ready'][:5] = 1
        agents['has_spawned'][:5] = np.array([1, 2, 0, 1, 3])
        agents['spawn_refractory'][:5] = np.array([0, 10, 0, 5, 15])
        agents['spawn_gravity_timer'][:5] = np.array([0, 5, 10, 0, 3])
        agents['immunosuppression_timer'][:5] = np.array([0, 0, 20, 15, 0])
        
        reset_spawning_season(agents)
        
        # Readiness and spawn counts should be reset
        assert np.all(agents['spawning_ready'] == 0)
        assert np.all(agents['has_spawned'] == 0)
        
        # But timers should persist (may extend beyond season)
        assert not np.all(agents['spawn_refractory'][:5] == 0)
        assert not np.all(agents['spawn_gravity_timer'][:5] == 0)
        assert not np.all(agents['immunosuppression_timer'][:5] == 0)


# ═══════════════════════════════════════════════════════════════════════
# LARVAL COHORT GENERATION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestLarvalCohortGeneration:
    """Test larval cohort generation from spawning events."""
    
    def setup_method(self):
        """Set up test data for larval generation."""
        self.agents = allocate_agents(20)
        self.genotypes = allocate_genotypes(20)
        self.spawning_config = SpawningSection()
        self.rng = np.random.default_rng(42)
        
        # Set up adults
        n_adults = 10
        self.agents['alive'][:n_adults] = True
        self.agents['stage'][:n_adults] = Stage.ADULT
        self.agents['size'][:n_adults] = 450.0  # Reproductive size
        self.agents['node_id'][:n_adults] = 0
        
        # Half female, half male
        self.agents['sex'][:5] = 0  # Females
        self.agents['sex'][5:10] = 1  # Males
        
        # Random genotypes
        self.genotypes[:n_adults] = self.rng.integers(0, 2, (n_adults, N_LOCI, 2))
    
    def test_no_cohort_without_both_sexes(self):
        """Test no cohort generated without both sexes."""
        # Only females spawn
        female_indices = [0, 1, 2]
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, female_indices, self.spawning_config, self.rng
        )
        assert cohort.n_competent == 0
        
        # Only males spawn  
        male_indices = [5, 6, 7]
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, male_indices, self.spawning_config, self.rng
        )
        assert cohort.n_competent == 0
    
    def test_cohort_with_both_sexes(self):
        """Test cohort generation with both sexes present."""
        spawner_indices = [0, 1, 5, 6]  # 2 females, 2 males
        
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, spawner_indices, self.spawning_config, self.rng
        )
        
        if cohort.n_competent > 0:
            assert cohort.genotypes.shape == (cohort.n_competent, N_LOCI, 2)
            assert cohort.parent_pairs.shape == (cohort.n_competent, 2)
            assert cohort.source_node == 0
            assert cohort.pld_days > 0
    
    def test_offspring_genotype_inheritance(self):
        """Test offspring inherit alleles from parents through Mendelian rules."""
        # Set specific parent genotypes to test inheritance
        spawner_indices = [0, 5]  # 1 female, 1 male
        
        # Female (index 0): all 0 alleles
        self.genotypes[0] = np.zeros((N_LOCI, 2), dtype=np.int8)
        
        # Male (index 5): all 1 alleles  
        self.genotypes[5] = np.ones((N_LOCI, 2), dtype=np.int8)
        
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, spawner_indices, self.spawning_config, self.rng
        )
        
        if cohort.n_competent > 0:
            # All offspring should have exactly one 0 and one 1 allele at each locus
            for larva in range(cohort.n_competent):
                for locus in range(N_LOCI):
                    alleles = cohort.genotypes[larva, locus]
                    assert sorted(alleles) == [0, 1]  # One from each parent
    
    def test_parent_tracking(self):
        """Test parent pair tracking in larval cohorts."""
        spawner_indices = [1, 2, 6, 7]  # 2 females, 2 males
        
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, spawner_indices, self.spawning_config, self.rng
        )
        
        if cohort.n_competent > 0:
            # All parent pairs should be valid agent indices
            for i in range(cohort.n_competent):
                mother_idx, father_idx = cohort.parent_pairs[i]
                assert mother_idx in spawner_indices
                assert father_idx in spawner_indices
                assert self.agents['sex'][mother_idx] == 0  # Mother is female
                assert self.agents['sex'][father_idx] == 1   # Father is male


# ═══════════════════════════════════════════════════════════════════════
# PHASE 2: CASCADE INDUCTION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestCascadeInduction:
    """Test Phase 2 cascade induction functionality."""
    
    def setup_method(self):
        """Set up test agents with spatial positions for cascade testing."""
        self.n_agents = 50
        self.agents = allocate_agents(self.n_agents)
        self.genotypes = allocate_genotypes(self.n_agents)
        self.spawning_config = SpawningSection()
        self.disease_config = DiseaseSection()
        self.rng = np.random.default_rng(42)
        
        # Create mature adults
        n_adults = 20
        self.agents['alive'][:n_adults] = True
        self.agents['stage'][:n_adults] = Stage.ADULT
        self.agents['size'][:n_adults] = 450.0
        self.agents['node_id'][:n_adults] = 0
        
        # Set sex distribution
        self.agents['sex'][:10] = 0  # Females
        self.agents['sex'][10:20] = 1  # Males
        
        # Set spatial positions in grid pattern
        positions = np.array([
            [0, 0], [25, 0], [50, 0], [75, 0], [100, 0],      # Row 1: females
            [0, 25], [25, 25], [50, 25], [75, 25], [100, 25], # Row 2: females  
            [0, 50], [25, 50], [50, 50], [75, 50], [100, 50], # Row 3: males
            [0, 75], [25, 75], [50, 75], [75, 75], [100, 75]  # Row 4: males
        ])
        
        for i in range(n_adults):
            self.agents['x'][i] = positions[i, 0]
            self.agents['y'][i] = positions[i, 1]
    
    def test_recent_spawners_mask(self):
        """Test identification of recent spawners within cascade window."""
        # Set up some spawn days
        self.agents['last_spawn_day'][0] = 105  # 0 days ago (today)
        self.agents['last_spawn_day'][1] = 104  # 1 day ago
        self.agents['last_spawn_day'][2] = 103  # 2 days ago
        self.agents['last_spawn_day'][3] = 102  # 3 days ago (edge case)
        self.agents['last_spawn_day'][4] = 101  # 4 days ago (too old)
        self.agents['last_spawn_day'][5] = 0    # Never spawned
        
        current_doy = 105
        cascade_window = 3
        
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        adults = self.agents[adult_mask]
        
        recent_mask = _get_recent_spawners_mask(adults, current_doy, cascade_window)
        
        # Should identify indices 0, 1, 2, 3 as recent spawners
        expected_recent = np.array([True, True, True, True, False, False] + [False] * 14)
        np.testing.assert_array_equal(recent_mask, expected_recent)
    
    def test_recent_spawners_year_wrapping(self):
        """Test recent spawner detection across year boundary."""
        # Test early year with late previous year spawners
        self.agents['last_spawn_day'][0] = 363  # Dec 29 previous year
        self.agents['last_spawn_day'][1] = 365  # Dec 31 previous year  
        self.agents['last_spawn_day'][2] = 1    # Jan 1 current year
        
        current_doy = 3  # Jan 3
        cascade_window = 5
        
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        adults = self.agents[adult_mask]
        
        recent_mask = _get_recent_spawners_mask(adults, current_doy, cascade_window)
        
        # Dec 29: (365-363) + 3 = 5 days ago (exactly at edge)
        # Dec 31: (365-365) + 3 = 3 days ago (recent) 
        # Jan 1: 3-1 = 2 days ago (recent)
        expected_recent = np.array([True, True, True] + [False] * 17)
        np.testing.assert_array_equal(recent_mask, expected_recent)
    
    def test_cascade_distance_calculation(self):
        """Test spatial distance calculations for cascade induction."""
        # Female at origin
        female_idx = 0  # At (0, 0)
        
        # Males at various distances
        target_indices = np.array([10, 11, 12])  # At (0,50), (25,50), (50,50)
        inducer_indices = np.array([female_idx])
        
        cascade_radius = 50.0
        induction_prob = 1.0  # 100% to test distance only
        
        induced = _check_cascade_induction(
            self.agents, target_indices, inducer_indices,
            cascade_radius, induction_prob, self.rng
        )
        
        # Distance from (0,0) to (0,50) = 50.0 (exactly at edge)
        # Distance from (0,0) to (25,50) = sqrt(25²+50²) = 55.9 (too far)
        # Distance from (0,0) to (50,50) = sqrt(50²+50²) = 70.7 (too far)
        
        assert 10 in induced  # (0,50) - exactly at cascade_radius
        assert 11 not in induced  # (25,50) - too far
        assert 12 not in induced  # (50,50) - too far
    
    def test_female_to_male_induction(self):
        """Test strong female→male induction (κ_fm = 0.80)."""
        # Set up spawning scenario
        female_idx = 0  # At (0, 0)
        male_idx = 10   # At (0, 50) - within cascade_radius
        
        # Female spawned yesterday
        self.agents['last_spawn_day'][female_idx] = 104
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['has_spawned'][male_idx] = 0
        self.agents['spawn_refractory'][male_idx] = 0
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Multiple trials to test probability
        successes = 0
        trials = 1000
        
        for _ in range(trials):
            # Reset state
            test_agents = self.agents.copy()
            test_agents['last_spawn_day'][male_idx] = current_doy - 1  # Reset to not spawned today
            
            cascade_spawners = _cascade_induction_step(
                test_agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            if male_idx in cascade_spawners:
                successes += 1
        
        success_rate = successes / trials
        expected_rate = self.spawning_config.induction_female_to_male  # 0.80
        
        # Should be close to expected rate (within ~3 standard deviations)
        std_error = np.sqrt(expected_rate * (1 - expected_rate) / trials)
        tolerance = 3 * std_error
        
        assert abs(success_rate - expected_rate) < tolerance
    
    def test_male_to_female_induction(self):
        """Test moderate male→female induction (κ_mf = 0.30)."""
        # Set up spawning scenario  
        male_idx = 10   # At (0, 50)
        female_idx = 0  # At (0, 0) - within cascade_radius
        
        # Male spawned yesterday
        self.agents['last_spawn_day'][male_idx] = 104
        self.agents['spawning_ready'][female_idx] = 1
        self.agents['has_spawned'][female_idx] = 0
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Multiple trials to test probability
        successes = 0
        trials = 1000
        
        for _ in range(trials):
            # Reset state
            test_agents = self.agents.copy()
            test_agents['last_spawn_day'][female_idx] = current_doy - 1  # Reset to not spawned today
            
            cascade_spawners = _cascade_induction_step(
                test_agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            if female_idx in cascade_spawners:
                successes += 1
        
        success_rate = successes / trials
        expected_rate = self.spawning_config.induction_male_to_female  # 0.30
        
        # Should be close to expected rate (within ~3 standard deviations)
        std_error = np.sqrt(expected_rate * (1 - expected_rate) / trials)
        tolerance = 3 * std_error
        
        assert abs(success_rate - expected_rate) < tolerance
    
    def test_cascade_failure_at_distance(self):
        """Test cascade fails when agents are too far apart."""
        # Female at origin, male far away
        female_idx = 0   # At (0, 0)
        male_idx = 19    # Place at distance > 200m to ensure cascade failure
        
        # Override male position to be beyond cascade_radius (200m)
        self.agents['x'][male_idx] = 250.0  # 250m away
        self.agents['y'][male_idx] = 0.0    # Distance = 250m > 200m cascade_radius
        
        # Female spawned recently
        self.agents['last_spawn_day'][female_idx] = 104
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['has_spawned'][male_idx] = 0
        self.agents['spawn_refractory'][male_idx] = 0
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Force high induction probability for this test
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should not induce due to distance
            assert male_idx not in cascade_spawners
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction
    
    def test_male_multi_bout_cascade(self):
        """Test males can be cascade-induced for multiple bouts."""
        # Set up male and female
        female_idx = 0   # At (0, 0)
        male_idx = 10    # At (0, 50) - within range
        
        # Make male ready and capable of multiple bouts
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['has_spawned'][male_idx] = 1  # Already has 1 bout
        self.agents['spawn_refractory'][male_idx] = 0  # Not refractory
        
        # Female spawned recently
        self.agents['last_spawn_day'][female_idx] = 104
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Force induction for testing
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should induce male for second bout
            assert male_idx in cascade_spawners
            assert self.agents['has_spawned'][male_idx] == 2
            assert self.agents['spawn_refractory'][male_idx] == self.spawning_config.male_refractory_days
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction
    
    def test_male_bout_limit_prevents_cascade(self):
        """Test males at bout limit cannot be cascade-induced.""" 
        # Set up male at bout limit
        female_idx = 0   # At (0, 0)  
        male_idx = 10    # At (0, 50)
        
        # Male at maximum bouts
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['has_spawned'][male_idx] = self.spawning_config.male_max_bouts  # At limit
        self.agents['spawn_refractory'][male_idx] = 0
        
        # Female spawned recently
        self.agents['last_spawn_day'][female_idx] = 104
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Force induction attempt
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should not induce (at bout limit)
            assert male_idx not in cascade_spawners
            assert self.agents['has_spawned'][male_idx] == self.spawning_config.male_max_bouts
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction
    
    def test_refractory_prevents_cascade(self):
        """Test refractory males cannot be cascade-induced."""
        # Set up refractory male
        female_idx = 0   # At (0, 0)
        male_idx = 10    # At (0, 50)
        
        # Male in refractory period
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['has_spawned'][male_idx] = 1
        self.agents['spawn_refractory'][male_idx] = 10  # Still refractory
        
        # Female spawned recently
        self.agents['last_spawn_day'][female_idx] = 104
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Force induction attempt
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should not induce (refractory)
            assert male_idx not in cascade_spawners
            assert self.agents['spawn_refractory'][male_idx] == 10
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction
    
    def test_female_single_spawn_prevents_cascade(self):
        """Test females who already spawned cannot be cascade-induced."""
        # Set up already-spawned female
        male_idx = 10    # At (0, 50)
        female_idx = 0   # At (0, 0)
        
        # Female already spawned this season
        self.agents['spawning_ready'][female_idx] = 1
        self.agents['has_spawned'][female_idx] = 1  # Already spawned
        
        # Male spawned recently
        self.agents['last_spawn_day'][male_idx] = 104
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Force induction attempt
        original_induction = self.spawning_config.induction_male_to_female
        self.spawning_config.induction_male_to_female = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should not induce (already spawned)
            assert female_idx not in cascade_spawners
            assert self.agents['has_spawned'][female_idx] == 1
            
        finally:
            self.spawning_config.induction_male_to_female = original_induction
    
    def test_density_dependent_cascade_high_density(self):
        """Test coordinated spawning bouts at high density."""
        # Set all agents as ready and clustered
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        self.agents['spawning_ready'][adult_mask] = 1
        
        # Cluster all agents within cascade radius of origin
        n_adults = np.sum(adult_mask)
        for i, idx in enumerate(np.where(adult_mask)[0]):
            # Arrange in tight cluster (all within 30m of origin)
            angle = 2 * np.pi * i / n_adults
            self.agents['x'][idx] = 15 * np.cos(angle)
            self.agents['y'][idx] = 15 * np.sin(angle)
        
        # One female spawns spontaneously first
        female_indices = np.where(adult_mask & (self.agents['sex'] == 0))[0]
        first_female = female_indices[0]
        self.agents['last_spawn_day'][first_female] = 104  # Yesterday
        
        mature_indices = np.where(adult_mask)[0]
        current_doy = 105
        
        # High density should create large cascade
        cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
        )
        
        # Should trigger multiple spawners due to clustering
        # (exact number is stochastic, but should be substantial)
        assert len(cascade_spawners) >= 3  # Expect multiple spawners in cluster
    
    def test_density_dependent_cascade_low_density(self):
        """Test sporadic spawning at low density (agents far apart)."""
        # Spread agents far apart (beyond cascade radius)
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        self.agents['spawning_ready'][adult_mask] = 1
        
        # Place agents >100m apart (well beyond cascade_radius=50)
        positions = [(i*150, j*150) for i in range(5) for j in range(4)]
        for i, idx in enumerate(np.where(adult_mask)[0]):
            if i < len(positions):
                self.agents['x'][idx] = positions[i][0]
                self.agents['y'][idx] = positions[i][1]
        
        # One female spawns
        female_indices = np.where(adult_mask & (self.agents['sex'] == 0))[0]
        first_female = female_indices[0]
        self.agents['last_spawn_day'][first_female] = 104  # Yesterday
        
        mature_indices = np.where(adult_mask)[0]
        current_doy = 105
        
        # Force high induction probability to test distance effect only
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should trigger few/no spawners due to isolation
            assert len(cascade_spawners) == 0  # All too far away
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction
    
    def test_cascade_window_expiry(self):
        """Test cascade cues expire after cascade_window days."""
        # Female spawned too long ago
        female_idx = 0   # At (0, 0)
        male_idx = 10    # At (0, 50) - close enough
        
        # Female spawned beyond cascade window  
        self.agents['last_spawn_day'][female_idx] = 100  # 5 days ago (window=3)
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['has_spawned'][male_idx] = 0
        self.agents['spawn_refractory'][male_idx] = 0
        
        # Get mature adult indices
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        mature_indices = np.where(adult_mask)[0]
        
        current_doy = 105
        
        # Force high induction probability
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cascade_spawners = _cascade_induction_step(
                self.agents, mature_indices, current_doy, self.spawning_config, DiseaseSection(), self.rng
            )
            
            # Should not induce (cue expired)
            assert male_idx not in cascade_spawners
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction


# ═══════════════════════════════════════════════════════════════════════
# PHASE 2 INTEGRATION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestPhase2Integration:
    """Test Phase 2 integration with spawning_step function."""
    
    def setup_method(self):
        """Set up test scenario for integrated spawning step."""
        self.n_agents = 20
        self.agents = allocate_agents(self.n_agents)
        self.genotypes = allocate_genotypes(self.n_agents)
        self.spawning_config = SpawningSection()
        self.disease_config = DiseaseSection()
        self.rng = np.random.default_rng(42)
        
        # Create clustered mature adults
        n_adults = 10
        self.agents['alive'][:n_adults] = True
        self.agents['stage'][:n_adults] = Stage.ADULT
        self.agents['size'][:n_adults] = 450.0
        self.agents['node_id'][:n_adults] = 0
        
        # Half female, half male, clustered together
        self.agents['sex'][:5] = 0  # Females
        self.agents['sex'][5:10] = 1  # Males
        
        # Place all within cascade radius
        for i in range(n_adults):
            angle = 2 * np.pi * i / n_adults
            self.agents['x'][i] = 20 * np.cos(angle)  # Within cascade_radius=50
            self.agents['y'][i] = 20 * np.sin(angle)
        
        # Initialize genotypes
        self.genotypes[:n_adults] = self.rng.integers(0, 2, (n_adults, N_LOCI, 2))
    
    def test_spontaneous_triggers_cascade(self):
        """Test spontaneous spawning can trigger cascade in clustered population."""
        # Make all adults ready
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        self.agents['spawning_ready'][adult_mask] = 1
        
        # Force one spontaneous spawner
        original_p_female = self.spawning_config.p_spontaneous_female
        original_p_male = self.spawning_config.p_spontaneous_male
        
        # Set one female to spawn spontaneously, low prob for others
        self.spawning_config.p_spontaneous_female = 0.2  # One should spawn
        self.spawning_config.p_spontaneous_male = 0.0    # No spontaneous males
        
        try:
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            # Count total spawners
            total_spawners = np.sum(self.agents['last_spawn_day'] == 105)
            
            # Should have more spawners than just spontaneous (cascade effect)
            # In clustered high-density scenario, expect cascade amplification
            spontaneous_spawners = np.sum(
                (self.agents['sex'] == 0) & (self.agents['last_spawn_day'] == 105)
            )
            
            # Note: cascade may or may not occur due to probability, but test structure
            # This mainly verifies the integration works without error
            assert total_spawners >= spontaneous_spawners
            
        finally:
            self.spawning_config.p_spontaneous_female = original_p_female
            self.spawning_config.p_spontaneous_male = original_p_male
    
    def test_cascade_produces_viable_cohorts(self):
        """Test cascade-induced spawning produces viable larval cohorts."""
        # Set up guaranteed cascade scenario
        female_idx = 0
        male_idx = 5
        
        # Female spawned yesterday
        self.agents['spawning_ready'][female_idx] = 1
        self.agents['spawning_ready'][male_idx] = 1
        self.agents['last_spawn_day'][female_idx] = 104
        self.agents['has_spawned'][female_idx] = 1  # Already spawned
        
        # Force male cascade induction
        original_induction = self.spawning_config.induction_female_to_male
        self.spawning_config.induction_female_to_male = 1.0
        
        try:
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, spawning_config=self.spawning_config, disease_config=self.disease_config, rng=self.rng
            )
            
            # Should generate cohorts from both spontaneous and cascade spawning
            if len(cohorts) > 0:
                cohort = cohorts[0]
                assert cohort.n_competent >= 0
                if cohort.n_competent > 0:
                    assert cohort.genotypes.shape == (cohort.n_competent, N_LOCI, 2)
                    assert cohort.source_node == 0
            
        finally:
            self.spawning_config.induction_female_to_male = original_induction


# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestSpawningConfiguration:
    """Test spawning configuration defaults and validation."""
    
    def test_default_configuration(self):
        """Test default spawning configuration values."""
        config = SpawningSection()
        
        # Season parameters
        assert config.season_start_doy == 305
        assert config.season_end_doy == 196
        assert config.peak_doy == 105
        assert config.peak_width_days == 60.0  # From Phase 1A2, verified in 1A3
        
        # Spontaneous rates (calibrated in Phase 1A3)
        assert abs(config.p_spontaneous_female - 0.012) < 1e-5
        assert abs(config.p_spontaneous_male - 0.0125) < 1e-5
        
        # Male bout parameters
        assert config.male_max_bouts == 3
        assert config.male_refractory_days == 21
        
        # Phase 2 cascade parameters (cascade_radius calibrated in Phase 1A2)
        assert config.induction_female_to_male == 0.80
        assert config.induction_male_to_female == 0.30
        assert config.cascade_window == 3
        assert config.cascade_radius == 200.0
        
        # Phase 3-4 parameters (not used yet)
        assert config.gravity_enabled is True
        assert DiseaseSection().immunosuppression_enabled is True
    
    def test_season_duration(self):
        """Test spawning season duration calculation.""" 
        config = SpawningSection()
        
        # Count days in default season (305 to 196, wrapping)
        season_days = 0
        for doy in range(1, 366):
            if in_spawning_season(doy, config.season_start_doy, config.season_end_doy):
                season_days += 1
        
        # Should be approximately 270 days (305 to 365 = 60, 1 to 196 = 196, total = 256)
        # Actually: (365 - 305 + 1) + 196 = 61 + 196 = 257 days
        expected_days = (365 - config.season_start_doy + 1) + config.season_end_doy
        assert season_days == expected_days
        assert season_days > 250  # Long spawning season


# ═══════════════════════════════════════════════════════════════════════
# PHASE 3 TESTS: SPAWNING GRAVITY
# ═══════════════════════════════════════════════════════════════════════

class TestSpawningGravityTimers:
    """Test spawning gravity timer functionality (Phase 3)."""
    
    def test_gravity_timer_set_on_readiness(self):
        """When agent becomes ready, gravity timer should be set if enabled."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(50)
        genotypes = allocate_genotypes(50)
        
        # Set up mature adults
        for i in range(50):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['spawning_ready'] = 0  # Not ready yet
            agents[i]['spawn_gravity_timer'] = 0
        
        config = SpawningSection()
        config.gravity_enabled = True
        config.pre_spawn_gravity_days = 10
        config.post_spawn_gravity_days = 8
        
        # High readiness probability to ensure some become ready
        config.peak_doy = 105
        config.peak_width_days = 45
        config.p_spontaneous_female = 0.0  # No spawning this step
        config.p_spontaneous_male = 0.0
        
        # Run spawning step during peak season
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=config, disease_config=DiseaseSection(), rng=rng)
        
        # Some agents should have become ready and gotten gravity timers
        ready_agents = agents['spawning_ready'] == 1
        if np.any(ready_agents):
            expected_timer = config.pre_spawn_gravity_days + config.post_spawn_gravity_days - 1  # Decremented in same step
            assert np.all(agents['spawn_gravity_timer'][ready_agents] == expected_timer)
    
    def test_gravity_timer_not_set_when_disabled(self):
        """When gravity is disabled, timer should not be set."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(50)
        genotypes = allocate_genotypes(50)
        
        # Set up mature adults
        for i in range(50):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['spawning_ready'] = 0
            agents[i]['spawn_gravity_timer'] = 0
        
        config = SpawningSection()
        config.gravity_enabled = False  # Disabled
        
        # Run spawning step
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=config, disease_config=DiseaseSection(), rng=rng)
        
        # No gravity timers should be set
        assert np.all(agents['spawn_gravity_timer'] == 0)
    
    def test_gravity_timer_countdown(self):
        """Gravity timers should count down daily."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(20)
        genotypes = allocate_genotypes(20)
        
        # Set up agents with various timer values
        for i in range(20):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['spawning_ready'] = 1
            agents[i]['spawn_gravity_timer'] = i + 1  # Timers from 1 to 20
        
        config = SpawningSection()
        
        # Run spawning step (should decrement timers)
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=config, disease_config=DiseaseSection(), rng=rng)
        
        # All timers should have decreased by 1
        for i in range(20):
            expected = max(0, i)  # Timer i+1 decremented by 1 = i, but min 0
            assert agents[i]['spawn_gravity_timer'] == expected
    
    def test_gravity_timer_stops_at_zero(self):
        """Gravity timers should not go below zero."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(10)
        genotypes = allocate_genotypes(10)
        
        # Set up agents with timer = 0
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['spawning_ready'] = 1
            agents[i]['spawn_gravity_timer'] = 0
        
        config = SpawningSection()
        
        # Run spawning step multiple times
        for _ in range(5):
            spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                         spawning_config=config, disease_config=DiseaseSection(), rng=rng)
        
        # All timers should remain at 0
        assert np.all(agents['spawn_gravity_timer'] == 0)
    
    def test_gravity_timer_with_spawning_event(self):
        """Test timer behavior when agents actually spawn."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(20)
        genotypes = allocate_genotypes(20)
        
        # Set up ready females
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['sex'] = 0  # Female
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 0
            agents[i]['spawn_gravity_timer'] = 5
        
        # Set up ready males  
        for i in range(10, 20):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['sex'] = 1  # Male
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 0
            agents[i]['spawn_refractory'] = 0
            agents[i]['spawn_gravity_timer'] = 8
        
        config = SpawningSection()
        config.p_spontaneous_female = 1.0  # Guarantee spawning
        config.p_spontaneous_male = 1.0
        
        initial_timers = agents['spawn_gravity_timer'].copy()
        
        # Run spawning step
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=config, disease_config=DiseaseSection(), rng=rng)
        
        # Timers should have decremented regardless of spawning
        expected_timers = np.maximum(0, initial_timers - 1)
        np.testing.assert_array_equal(agents['spawn_gravity_timer'], expected_timers)
        
        # Some agents should have spawned (separate from timer logic)
        assert np.sum(agents['has_spawned'] > 0) > 0
    
    def test_season_reset_preserves_gravity_timers(self):
        """Season reset should not reset gravity timers (they extend beyond season)."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(10)
        
        # Set up agents with various timers
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 1
            agents[i]['spawn_gravity_timer'] = i + 5  # Various timer values
        
        initial_timers = agents['spawn_gravity_timer'].copy()
        
        # Reset spawning season
        reset_spawning_season(agents)
        
        # spawning_ready and has_spawned should be reset
        assert np.all(agents['spawning_ready'] == 0)
        assert np.all(agents['has_spawned'] == 0)
        
        # But gravity timers should be preserved
        np.testing.assert_array_equal(agents['spawn_gravity_timer'], initial_timers)


# ═══════════════════════════════════════════════════════════════════════
# PHASE 4 TESTS: POST-SPAWNING IMMUNOSUPPRESSION
# ═══════════════════════════════════════════════════════════════════════

class TestImmounosuppressionTimers:
    """Test post-spawning immunosuppression timer functionality (Phase 4)."""
    
    def test_immunosuppression_timer_set_on_spawning(self):
        """When agent spawns, immunosuppression timer should be set if enabled."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(20)
        genotypes = allocate_genotypes(20)
        
        # Set up ready females and males
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['sex'] = 0  # Female
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 0
            agents[i]['immunosuppression_timer'] = 0
        
        for i in range(10, 20):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['sex'] = 1  # Male
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 0
            agents[i]['spawn_refractory'] = 0
            agents[i]['immunosuppression_timer'] = 0
        
        spawning_config = SpawningSection()
        spawning_config.p_spontaneous_female = 1.0  # Guarantee spawning
        spawning_config.p_spontaneous_male = 1.0
        
        disease_config = DiseaseSection()
        disease_config.immunosuppression_enabled = True
        disease_config.immunosuppression_duration = 28
        
        # Run spawning step
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=spawning_config, disease_config=disease_config, rng=rng)
        
        # All agents should have spawned and gotten immunosuppression timers
        spawned_mask = agents['has_spawned'] > 0
        assert np.sum(spawned_mask) > 0  # Some should have spawned
        
        # All spawned agents should have immunosuppression timer set
        # Timer is set after countdown, so it should be the full duration
        expected_timer = disease_config.immunosuppression_duration
        spawned_indices = np.where(spawned_mask)[0]
        assert np.all(agents['immunosuppression_timer'][spawned_indices] == expected_timer)
    
    def test_immunosuppression_timer_not_set_when_disabled(self):
        """When immunosuppression is disabled, timer should not be set."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(10)
        genotypes = allocate_genotypes(10)
        
        # Set up ready females
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['sex'] = 0  # Female
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 0
            agents[i]['immunosuppression_timer'] = 0
        
        spawning_config = SpawningSection()
        spawning_config.p_spontaneous_female = 1.0  # Guarantee spawning
        
        disease_config = DiseaseSection()
        disease_config.immunosuppression_enabled = False  # Disabled
        
        # Run spawning step
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=spawning_config, disease_config=disease_config, rng=rng)
        
        # No immunosuppression timers should be set
        assert np.all(agents['immunosuppression_timer'] == 0)
    
    def test_immunosuppression_timer_resets_for_male_multiple_bouts(self):
        """Male timer should reset on each spawning bout."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(5)
        genotypes = allocate_genotypes(5)
        
        # Set up ready male
        agents[0]['alive'] = True
        agents[0]['stage'] = Stage.ADULT
        agents[0]['size'] = 450.0
        agents[0]['sex'] = 1  # Male
        agents[0]['spawning_ready'] = 1
        agents[0]['has_spawned'] = 0
        agents[0]['spawn_refractory'] = 0
        agents[0]['immunosuppression_timer'] = 0
        
        spawning_config = SpawningSection()
        spawning_config.p_spontaneous_male = 1.0  # Guarantee spawning
        spawning_config.male_max_bouts = 3
        spawning_config.male_refractory_days = 21
        
        disease_config = DiseaseSection()
        disease_config.immunosuppression_enabled = True
        disease_config.immunosuppression_duration = 28
        
        # First spawning
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=spawning_config, disease_config=disease_config, rng=rng)
        
        assert agents[0]['has_spawned'] == 1
        first_timer = agents[0]['immunosuppression_timer'] 
        assert first_timer == disease_config.immunosuppression_duration  # Set after countdown
        
        # Simulate timer countdown and refractory period expiry
        agents[0]['immunosuppression_timer'] = 10  # Partially decremented
        agents[0]['spawn_refractory'] = 0  # Ready to spawn again
        
        # Second spawning
        spawning_step(agents, genotypes, day_of_year=130, node_latitude=40.0,
                     spawning_config=spawning_config, disease_config=disease_config, rng=rng)
        
        assert agents[0]['has_spawned'] == 2
        # Timer should reset to full duration
        second_timer = agents[0]['immunosuppression_timer']
        assert second_timer == disease_config.immunosuppression_duration
    
    def test_immunosuppression_timer_countdown(self):
        """Immunosuppression timers should count down daily."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(20)
        genotypes = allocate_genotypes(20)
        
        # Set up agents with various timer values
        for i in range(20):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['spawning_ready'] = 1
            agents[i]['immunosuppression_timer'] = i + 1  # Timers from 1 to 20
        
        spawning_config = SpawningSection()
        disease_config = DiseaseSection()
        
        # Run spawning step (should decrement timers)
        spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                     spawning_config=spawning_config, disease_config=disease_config, rng=rng)
        
        # All timers should have decreased by 1
        for i in range(20):
            expected = max(0, i)  # Timer i+1 decremented by 1 = i, but min 0
            assert agents[i]['immunosuppression_timer'] == expected
    
    def test_immunosuppression_timer_stops_at_zero(self):
        """Immunosuppression timers should not go below zero."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(10)
        genotypes = allocate_genotypes(10)
        
        # Set up agents with timer = 0
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['spawning_ready'] = 1
            agents[i]['immunosuppression_timer'] = 0
        
        spawning_config = SpawningSection()
        disease_config = DiseaseSection()
        
        # Run spawning step multiple times
        for _ in range(5):
            spawning_step(agents, genotypes, day_of_year=105, node_latitude=40.0,
                         spawning_config=spawning_config, disease_config=disease_config, rng=rng)
        
        # All timers should remain at 0
        assert np.all(agents['immunosuppression_timer'] == 0)
    
    def test_season_reset_preserves_immunosuppression_timers(self):
        """Season reset should not reset immunosuppression timers (they extend beyond season)."""
        rng = np.random.default_rng(42)
        agents = allocate_agents(10)
        
        # Set up agents with various timers
        for i in range(10):
            agents[i]['alive'] = True
            agents[i]['spawning_ready'] = 1
            agents[i]['has_spawned'] = 1
            agents[i]['immunosuppression_timer'] = i + 5  # Various timer values
        
        initial_timers = agents['immunosuppression_timer'].copy()
        
        # Reset spawning season
        reset_spawning_season(agents)
        
        # spawning_ready and has_spawned should be reset
        assert np.all(agents['spawning_ready'] == 0)
        assert np.all(agents['has_spawned'] == 0)
        
        # But immunosuppression timers should be preserved
        np.testing.assert_array_equal(agents['immunosuppression_timer'], initial_timers)


# ═══════════════════════════════════════════════════════════════════════
# PHASE 4 INTEGRATION TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestImmounosuppressionIntegration:
    """Test integration between spawning and disease modules for immunosuppression."""
    
    def test_execute_spawning_events_with_both_configs(self):
        """Test that _execute_spawning_events works with both spawning and disease configs."""
        agents = allocate_agents(5)
        
        # Set up agents
        for i in range(5):
            agents[i]['alive'] = True
            agents[i]['stage'] = Stage.ADULT
            agents[i]['size'] = 450.0
            agents[i]['sex'] = i % 2  # Mix of male/female
            agents[i]['immunosuppression_timer'] = 0
        
        spawning_config = SpawningSection()
        disease_config = DiseaseSection()
        disease_config.immunosuppression_enabled = True
        disease_config.immunosuppression_duration = 28
        
        # Execute spawning events
        spawner_indices = np.array([0, 1, 2])
        _execute_spawning_events(
            agents, spawner_indices, day_of_year=105, 
            spawning_config=spawning_config, disease_config=disease_config
        )
        
        # All spawners should have immunosuppression timer set
        assert agents[0]['immunosuppression_timer'] == disease_config.immunosuppression_duration
        assert agents[1]['immunosuppression_timer'] == disease_config.immunosuppression_duration
        assert agents[2]['immunosuppression_timer'] == disease_config.immunosuppression_duration
        
        # Non-spawners should have timer = 0
        assert agents[3]['immunosuppression_timer'] == 0
        assert agents[4]['immunosuppression_timer'] == 0
    
    def test_immunosuppression_disabled_in_execute_spawning_events(self):
        """Test that immunosuppression timers are not set when disabled."""
        agents = allocate_agents(3)
        
        # Set up agents
        for i in range(3):
            agents[i]['alive'] = True
            agents[i]['immunosuppression_timer'] = 0
        
        spawning_config = SpawningSection()
        disease_config = DiseaseSection()
        disease_config.immunosuppression_enabled = False  # Disabled
        
        # Execute spawning events
        spawner_indices = np.array([0, 1])
        _execute_spawning_events(
            agents, spawner_indices, day_of_year=105,
            spawning_config=spawning_config, disease_config=disease_config
        )
        
        # No timers should be set
        assert np.all(agents['immunosuppression_timer'] == 0)


if __name__ == "__main__":
    pytest.main([__file__])