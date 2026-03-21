#!/usr/bin/env python3
"""
Verify calibrated spawning parameters produce correct biological behavior.

Uses simplified rate-based simulation to test spawning targets.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from pathlib import Path

from sswd_evoepi import config


def calculate_seasonal_modifier(doy, season_start_doy, season_end_doy, peak_doy, peak_width):
    """Calculate Gaussian seasonal modifier for a given day of year."""
    # Check if in season (handle year wrap)
    in_season = False
    if season_start_doy > season_end_doy:  # Wraps year (Nov-Jul)
        in_season = (doy >= season_start_doy) or (doy <= season_end_doy)
    else:
        in_season = (season_start_doy <= doy <= season_end_doy)
    
    if not in_season:
        return 0.0
    
    # Calculate distance from peak (handle year wrap)
    if peak_doy < season_start_doy:  # Peak in next year (Apr < Nov)
        if doy >= season_start_doy:  # Currently in first part (Nov-Dec)
            days_from_peak = (365 - doy) + peak_doy
        else:  # Currently in second part (Jan-Jul)
            days_from_peak = abs(doy - peak_doy)
    else:
        days_from_peak = abs(doy - peak_doy)
    
    # Gaussian seasonal modifier
    return np.exp(-(days_from_peak**2) / (2 * peak_width**2))


def simulate_spawning_season(cfg, n_males=250, n_females=250, seed=42):
    """Simulate spawning behavior over a full year."""
    np.random.seed(seed)
    
    # Get parameters
    season_start_doy = cfg.spawning.season_start_doy  # 305 (Nov 1)
    season_end_doy = cfg.spawning.season_end_doy      # 196 (Jul 15)
    peak_doy = cfg.spawning.peak_doy                  # 105 (Apr 15)
    peak_width = cfg.spawning.peak_width_days         # 60
    p_female = cfg.spawning.p_spontaneous_female      # 0.004205
    p_male = cfg.spawning.p_spontaneous_male          # 0.006308
    
    # Track results
    male_bout_counts = [0] * n_males
    females_spawned = set()
    daily_spawning = []
    total_spawning_events = 0
    season_spawning = 0
    major_days = 0
    
    print(f"Parameters: p_female={p_female:.6f}, p_male={p_male:.6f}, peak_width={peak_width}")
    print(f"Season: {season_start_doy}-{season_end_doy} (wraps), peak at {peak_doy}")
    
    for day in range(365):
        doy = day + 1
        
        # Calculate seasonal modifier
        seasonal_modifier = calculate_seasonal_modifier(doy, season_start_doy, season_end_doy, peak_doy, peak_width)
        
        if seasonal_modifier == 0:
            daily_spawning.append(0)
            continue
        
        # Apply spawning probabilities
        p_female_today = p_female * seasonal_modifier
        p_male_today = p_male * seasonal_modifier
        
        # Sample spawning events
        male_spawns = np.random.binomial(n_males, p_male_today)
        female_spawns = np.random.binomial(n_females, p_female_today)
        
        # Update counters
        if male_spawns > 0:
            spawning_males = np.random.choice(n_males, male_spawns, replace=False)
            for male_id in spawning_males:
                male_bout_counts[male_id] += 1
        
        if female_spawns > 0:
            spawning_females = np.random.choice(n_females, female_spawns, replace=False)
            for female_id in spawning_females:
                females_spawned.add(female_id)
        
        daily_total = male_spawns + female_spawns
        daily_spawning.append(daily_total)
        total_spawning_events += daily_total
        season_spawning += daily_total
        
        if daily_total > 10:
            major_days += 1
        
        # Progress for peak days
        if seasonal_modifier > 0.5:  # Near peak
            print(f"  Day {day} (DOY {doy}): modifier={seasonal_modifier:.3f}, spawning={daily_total}")
    
    return {
        'male_bout_counts': male_bout_counts,
        'females_spawned': females_spawned,
        'daily_spawning': daily_spawning,
        'total_spawning': total_spawning_events,
        'major_days': major_days
    }


def verify_targets(sim_results, cfg):
    """Verify all spawning targets against simulation results."""
    results = {}
    
    # a) Male mean bouts: 1.8-2.6 (target 2.2)
    male_bouts = [count for count in sim_results['male_bout_counts'] if count > 0]  # Only males that spawned
    mean_bouts = np.mean(sim_results['male_bout_counts'])  # All males (including 0s)
    results['male_mean_bouts'] = {
        'value': mean_bouts,
        'target': (1.8, 2.6),
        'pass': 1.8 <= mean_bouts <= 2.6,
        'spawning_males': len(male_bouts),
        'mean_spawning_males': np.mean(male_bouts) if male_bouts else 0
    }
    
    # b) Major spawning days: 2-5 per season
    major_days = sim_results['major_days']
    results['major_spawning_days'] = {
        'value': major_days,
        'target': (2, 5),
        'pass': 2 <= major_days <= 5
    }
    
    # d) Female participation: >80%
    n_females = 250  # Known from simulation
    participation = len(sim_results['females_spawned']) / n_females
    results['female_participation'] = {
        'value': participation,
        'target': (0.8, 1.0),
        'pass': participation >= 0.8,
        'spawning_females': len(sim_results['females_spawned'])
    }
    
    # Additional diagnostics
    results['total_spawning_events'] = sim_results['total_spawning']
    results['peak_daily_spawning'] = max(sim_results['daily_spawning'])
    
    return results


def create_diagnostic_plot(sim_results, cfg):
    """Create diagnostic plot of spawning patterns."""
    plt.style.use('dark_background')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # Top-left: Daily spawning events
    ax = axes[0, 0]
    ax.plot(range(365), sim_results['daily_spawning'], color='cyan', linewidth=2)
    ax.axhline(10, color='red', linestyle='--', alpha=0.7, label='Major Day Threshold')
    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Spawning Events')
    ax.set_title('Daily Spawning Events')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Top-right: Male bout distribution
    ax = axes[0, 1]
    bout_counts = sim_results['male_bout_counts']
    max_bouts = max(bout_counts) if bout_counts else 0
    if max_bouts > 0:
        bins = range(0, max_bouts + 2)
        ax.hist(bout_counts, bins=bins, alpha=0.7, color='orange', edgecolor='white')
        mean_bouts = np.mean(bout_counts)
        ax.axvline(mean_bouts, color='red', linestyle='--', label=f'Mean: {mean_bouts:.1f}')
        ax.axvspan(1.8, 2.6, alpha=0.2, color='green', label='Target Range')
    ax.set_xlabel('Bouts per Male')
    ax.set_ylabel('Count')
    ax.set_title('Male Spawning Bout Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Bottom-left: Seasonal pattern
    ax = axes[1, 0]
    # Show seasonal modifier curve
    season_start = cfg.spawning.season_start_doy
    season_end = cfg.spawning.season_end_doy
    peak_doy = cfg.spawning.peak_doy
    peak_width = cfg.spawning.peak_width_days
    
    modifiers = [calculate_seasonal_modifier(doy, season_start, season_end, peak_doy, peak_width) 
                for doy in range(1, 366)]
    ax.plot(range(1, 366), modifiers, color='orange', label='Seasonal Modifier')
    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Seasonal Modifier')
    ax.set_title('Spawning Season Pattern')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Bottom-right: Participation summary
    ax = axes[1, 1]
    n_males, n_females = 250, 250
    spawning_males = sum(1 for count in bout_counts if count > 0)
    spawning_females = len(sim_results['females_spawned'])
    
    categories = ['Males\nSpawning', 'Males\nNot Spawning', 'Females\nSpawning', 'Females\nNot Spawning']
    values = [spawning_males, n_males - spawning_males, spawning_females, n_females - spawning_females]
    colors = ['lightblue', 'darkblue', 'lightcoral', 'darkred']
    
    ax.bar(categories, values, color=colors, alpha=0.7)
    ax.set_ylabel('Count')
    ax.set_title('Spawning Participation')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add percentage labels
    for i, (cat, val) in enumerate(zip(categories, values)):
        if 'Males' in cat:
            pct = val / n_males * 100
        else:
            pct = val / n_females * 100
        ax.text(i, val + 5, f'{pct:.1f}%', ha='center', va='bottom')
    
    plt.tight_layout()
    return fig


def main():
    """Main verification workflow."""
    print("SSWD-EvoEpi Spawning Calibration Verification")
    print("=" * 50)
    
    # Create results directory
    Path("results/calibration").mkdir(parents=True, exist_ok=True)
    
    # Load config
    cfg = config.default_config()
    
    # Run simulation
    print("Running spawning simulation...")
    sim_results = simulate_spawning_season(cfg, seed=42)
    
    # Verify targets
    verification = verify_targets(sim_results, cfg)
    
    # Print results
    print("\n" + "="*60)
    print("CALIBRATION VERIFICATION RESULTS")
    print("="*60)
    
    all_pass = True
    for target_name, result in verification.items():
        if isinstance(result, dict) and 'pass' in result:
            status = "PASS" if result['pass'] else "FAIL"
            if not result['pass']:
                all_pass = False
            
            print(f"{target_name:20s}: {result['value']:8.3f} "
                  f"(target: {result['target'][0]:.1f}-{result['target'][1]:.1f}) [{status}]")
            
            # Extra info for male bouts
            if target_name == 'male_mean_bouts':
                print(f"{'':21s}  {result['spawning_males']} males spawned, "
                      f"mean among spawners: {result['mean_spawning_males']:.2f}")
            elif target_name == 'female_participation':
                print(f"{'':21s}  {result['spawning_females']}/250 females spawned")
    
    print(f"{'total_events':20s}: {verification['total_spawning_events']:8d}")
    print(f"{'peak_daily':20s}: {verification['peak_daily_spawning']:8d}")
    
    print("="*60)
    print(f"OVERALL STATUS: {'PASS' if all_pass else 'FAIL'}")
    print("="*60)
    
    # Create diagnostic plot
    fig = create_diagnostic_plot(sim_results, cfg)
    plot_path = "results/calibration/calibration_summary.png"
    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved: {plot_path}")
    
    # Generate report
    report_path = "results/calibration/CALIBRATION_REPORT.md"
    with open(report_path, 'w') as f:
        f.write("# SSWD-EvoEpi Spawning Calibration Verification\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Final Parameter Values\n\n")
        f.write(f"- `p_spontaneous_female`: {cfg.spawning.p_spontaneous_female:.6f}\n")
        f.write(f"- `p_spontaneous_male`: {cfg.spawning.p_spontaneous_male:.6f}\n")
        f.write(f"- `peak_width_days`: {cfg.spawning.peak_width_days}\n")
        f.write(f"- `peak_doy`: {cfg.spawning.peak_doy} (April 15)\n")
        f.write(f"- `season_start_doy`: {cfg.spawning.season_start_doy} (Nov 1)\n")
        f.write(f"- `season_end_doy`: {cfg.spawning.season_end_doy} (Jul 15)\n\n")
        
        f.write("## Target Verification Results\n\n")
        for target_name, result in verification.items():
            if isinstance(result, dict) and 'pass' in result:
                status_emoji = "✅" if result['pass'] else "❌"
                f.write(f"- **{target_name}**: {result['value']:.3f} "
                       f"(target: {result['target'][0]:.1f}-{result['target'][1]:.1f}) {status_emoji}\n")
        
        f.write(f"\n## Simulation Summary\n\n")
        f.write(f"- **Total spawning events**: {verification['total_spawning_events']}\n")
        f.write(f"- **Peak daily spawning**: {verification['peak_daily_spawning']}\n")
        f.write(f"- **Major spawning days**: {verification['major_spawning_days']['value']}\n")
        f.write(f"- **Males that spawned**: {verification['male_mean_bouts']['spawning_males']}/250\n")
        f.write(f"- **Females that spawned**: {verification['female_participation']['spawning_females']}/250\n\n")
        
        f.write(f"## Overall Status\n\n")
        overall_emoji = "✅" if all_pass else "❌"
        f.write(f"{overall_emoji} **{'PASS' if all_pass else 'FAIL'}** - {'All' if all_pass else 'Some'} targets met\n\n")
        
        f.write("## Diagnostic Plot\n\n")
        f.write("![Calibration Summary](calibration_summary.png)\n\n")
        
        f.write("## Methods\n\n")
        f.write("- **Approach**: Direct parameter simulation (simplified model)\n")
        f.write("- **Population**: 250 males + 250 females\n")
        f.write("- **Duration**: 365 days\n")
        f.write("- **Seed**: 42\n")
        f.write("- **Seasonal model**: Gaussian around peak with year-wrap handling\n")
    
    print(f"Report saved: {report_path}")
    
    return all_pass


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)