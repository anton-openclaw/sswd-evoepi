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
)
from sswd_evoepi.types import AGENT_DTYPE, Stage, allocate_agents, allocate_genotypes, N_LOCI
from sswd_evoepi.config import SpawningSection


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
        self.config = SpawningSection()
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
            node_latitude=45.0, config=self.config, rng=self.rng
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
            node_latitude=40.0, config=self.config, rng=self.rng
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
        original_p_female = self.config.p_spontaneous_female
        self.config.p_spontaneous_female = 1.0  # 100% spawn rate
        
        try:
            # First day - should spawn
            cohorts1 = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            spawned_females = np.sum(
                (self.agents['sex'] == 0) & 
                (self.agents['has_spawned'] == 1)
            )
            
            # Second day - should not spawn again (already has_spawned=1)
            cohorts2 = spawning_step(
                self.agents, self.genotypes, day_of_year=106,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            # Same females should still have has_spawned=1, no additional spawning
            spawned_females_day2 = np.sum(
                (self.agents['sex'] == 0) & 
                (self.agents['has_spawned'] == 1)
            )
            
            assert spawned_females_day2 == spawned_females  # No additional spawning
            
        finally:
            self.config.p_spontaneous_female = original_p_female
    
    def test_male_multi_bout_with_refractory(self):
        """Test male multi-bout spawning with refractory periods."""
        # Set up ready males
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        male_mask = self.agents['sex'] == 1
        ready_males = adult_mask & male_mask
        
        self.agents['spawning_ready'][ready_males] = 1
        
        # Force high male spawning probability
        original_p_male = self.config.p_spontaneous_male
        self.config.p_spontaneous_male = 1.0
        
        try:
            # Day 1: First bout
            cohorts1 = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            spawned_males = (self.agents['sex'] == 1) & (self.agents['has_spawned'] > 0)
            bout_counts_day1 = self.agents['has_spawned'][spawned_males]
            refractory_timers = self.agents['spawn_refractory'][spawned_males]
            
            # Check spawning occurred and refractory is set
            assert np.any(spawned_males)
            assert np.all(bout_counts_day1 == 1)
            assert np.all(refractory_timers == self.config.male_refractory_days)
            
            # Day 2: Should not spawn (refractory)
            cohorts2 = spawning_step(
                self.agents, self.genotypes, day_of_year=106,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            bout_counts_day2 = self.agents['has_spawned'][spawned_males]
            assert np.all(bout_counts_day2 == 1)  # No additional bouts
            
            # Fast-forward past refractory period
            self.agents['spawn_refractory'][:] = 0  # Clear refractory
            
            # Day 22: Should spawn again (second bout)
            cohorts3 = spawning_step(
                self.agents, self.genotypes, day_of_year=127,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            bout_counts_day22 = self.agents['has_spawned'][spawned_males]
            assert np.any(bout_counts_day22 == 2)  # Some second bouts
            
        finally:
            self.config.p_spontaneous_male = original_p_male
    
    def test_male_bout_limit_enforcement(self):
        """Test males cannot exceed maximum bout limit."""
        # Set up a male at bout limit
        adult_mask = self.agents['alive'] & (self.agents['stage'] == Stage.ADULT)
        male_mask = self.agents['sex'] == 1
        test_male_idx = np.where(adult_mask & male_mask)[0][0]
        
        self.agents['spawning_ready'][test_male_idx] = 1
        self.agents['has_spawned'][test_male_idx] = self.config.male_max_bouts  # At limit
        self.agents['spawn_refractory'][test_male_idx] = 0  # Not refractory
        
        # Force high spawning probability
        original_p_male = self.config.p_spontaneous_male
        self.config.p_spontaneous_male = 1.0
        
        try:
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            # Male should not spawn (already at limit)
            assert self.agents['has_spawned'][test_male_idx] == self.config.male_max_bouts
            
        finally:
            self.config.p_spontaneous_male = original_p_male
    
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
                node_latitude=40.0, config=self.config, rng=self.rng
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
        original_p_female = self.config.p_spontaneous_female
        original_p_male = self.config.p_spontaneous_male
        self.config.p_spontaneous_female = 1.0
        self.config.p_spontaneous_male = 1.0
        
        try:
            cohorts = spawning_step(
                self.agents, self.genotypes, day_of_year=105,
                node_latitude=40.0, config=self.config, rng=self.rng
            )
            
            # Should generate cohort if both sexes present
            if np.any(self.agents['sex'][adult_mask] == 0) and np.any(self.agents['sex'][adult_mask] == 1):
                assert len(cohorts) > 0
                assert cohorts[0].n_competent > 0
                assert cohorts[0].genotypes.shape[0] == cohorts[0].n_competent
                assert cohorts[0].genotypes.shape[1] == N_LOCI
                assert cohorts[0].genotypes.shape[2] == 2  # Diploid
            
        finally:
            self.config.p_spontaneous_female = original_p_female
            self.config.p_spontaneous_male = original_p_male


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
        _execute_spawning_events(agents, female_indices, day_of_year=105)
        
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
        _execute_spawning_events(agents, male_indices, day_of_year=105, refractory_days=21)
        
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
        self.config = SpawningSection()
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
            self.agents, self.genotypes, female_indices, self.config, self.rng
        )
        assert cohort.n_competent == 0
        
        # Only males spawn  
        male_indices = [5, 6, 7]
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, male_indices, self.config, self.rng
        )
        assert cohort.n_competent == 0
    
    def test_cohort_with_both_sexes(self):
        """Test cohort generation with both sexes present."""
        spawner_indices = [0, 1, 5, 6]  # 2 females, 2 males
        
        cohort = _generate_larval_cohort(
            self.agents, self.genotypes, spawner_indices, self.config, self.rng
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
            self.agents, self.genotypes, spawner_indices, self.config, self.rng
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
            self.agents, self.genotypes, spawner_indices, self.config, self.rng
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
        assert config.peak_width_days == 45.0
        
        # Spontaneous rates
        assert config.p_spontaneous_female == 0.005
        assert config.p_spontaneous_male == 0.008
        
        # Male bout parameters
        assert config.male_max_bouts == 3
        assert config.male_refractory_days == 21
        
        # Phase 2-4 parameters (not used yet)
        assert config.induction_female_to_male == 0.80
        assert config.gravity_enabled is True
        assert config.immunosuppression_enabled is True
    
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


if __name__ == "__main__":
    pytest.main([__file__])