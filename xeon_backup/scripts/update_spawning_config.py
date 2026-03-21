#!/usr/bin/env python3
"""
Update spawning configuration based on calibration results.

Based on Phase 1A2 calibration sweep results, the spawning probabilities 
are too high by ~4-16x. This script calculates and applies corrected values.

Author: Anton 🔬
Date: 2026-02-15
"""

import json
from pathlib import Path
from sswd_evoepi.config import SimulationConfig

def update_spawning_parameters():
    """Update spawning config based on calibration results."""
    
    print("=" * 60)
    print("SPAWNING CONFIG UPDATE - Phase 1A2")
    print("=" * 60)
    
    # Load calibration results
    results_dir = Path("results/spawning_calibration")
    best_params_file = sorted(results_dir.glob("best_parameters_*.json"))[-1]
    
    with open(best_params_file, 'r') as f:
        calibration_results = json.load(f)
    
    print(f"Loaded calibration results from: {best_params_file.name}")
    print(f"Best calibrated values (90-day simulation, n=200):")
    for key, value in calibration_results.items():
        print(f"  {key}: {value}")
    
    # Calculate scaling factors
    # Target: 2.2 male bouts per full season (270 days)
    # Achieved: 3.0 male bouts in 90-day simulation
    # Scale factor for male spawning: (0.8 target / 3.0 achieved) = 0.267
    
    male_bout_scale = 0.8 / calibration_results['male_bouts_achieved']  # 0.8 was 90-day target
    
    # For major bouts: achieved 13, target was 2 → scale by 2/13 = 0.154
    major_bout_scale = 2.0 / calibration_results['major_bouts_achieved']
    
    # For larvae: achieved 79974, target was 5000 → scale by 5000/79974 = 0.063
    larvae_scale = 5000 / calibration_results['larvae_achieved']
    
    print(f"\nScaling factors:")
    print(f"  Male bout scale: {male_bout_scale:.3f}")
    print(f"  Major bout scale: {major_bout_scale:.3f}")
    print(f"  Larvae scale: {larvae_scale:.3f}")
    
    # Use conservative scaling (average of male bout and major bout scales)
    overall_scale = (male_bout_scale + major_bout_scale) / 2
    print(f"  Overall scale (conservative): {overall_scale:.3f}")
    
    # Apply scaling to calibrated parameters
    new_p_female = calibration_results['p_spontaneous_female'] * overall_scale
    new_p_male = calibration_results['p_spontaneous_male'] * overall_scale
    
    # Keep other parameters as they performed reasonably
    new_cascade_radius = calibration_results['cascade_radius']
    new_peak_width = calibration_results['peak_width_days']
    
    print(f"\nNew parameter values:")
    print(f"  p_spontaneous_female: {new_p_female:.5f} (was {calibration_results['p_spontaneous_female']:.3f})")
    print(f"  p_spontaneous_male: {new_p_male:.5f} (was {calibration_results['p_spontaneous_male']:.3f})")
    print(f"  cascade_radius: {new_cascade_radius:.1f} (unchanged)")
    print(f"  peak_width_days: {new_peak_width:.1f} (unchanged)")
    
    # Load current config to see changes
    config = SimulationConfig()
    current_spawning = config.spawning
    
    print(f"\nCurrent config values:")
    print(f"  p_spontaneous_female: {current_spawning.p_spontaneous_female}")
    print(f"  p_spontaneous_male: {current_spawning.p_spontaneous_male}")
    print(f"  cascade_radius: {current_spawning.cascade_radius}")
    print(f"  peak_width_days: {current_spawning.peak_width_days}")
    
    # Calculate change ratios
    female_ratio = new_p_female / current_spawning.p_spontaneous_female
    male_ratio = new_p_male / current_spawning.p_spontaneous_male
    
    print(f"\nChange ratios (new/old):")
    print(f"  p_spontaneous_female: {female_ratio:.3f}x")
    print(f"  p_spontaneous_male: {male_ratio:.3f}x")
    
    # Update the config file directly
    config_file = Path("sswd_evoepi/config.py")
    
    print(f"\nUpdating {config_file}...")
    
    # Read current config file
    with open(config_file, 'r') as f:
        config_content = f.read()
    
    # Find and replace spawning parameters
    import re
    
    # Replace p_spontaneous_female
    pattern = r'p_spontaneous_female:\s*float\s*=\s*[\d\.]+\s*'
    replacement = f'p_spontaneous_female: float = {new_p_female:.6f}  '
    config_content = re.sub(pattern, replacement, config_content)
    
    # Replace p_spontaneous_male
    pattern = r'p_spontaneous_male:\s*float\s*=\s*[\d\.]+\s*'
    replacement = f'p_spontaneous_male: float = {new_p_male:.6f}  '
    config_content = re.sub(pattern, replacement, config_content)
    
    # Replace cascade_radius
    pattern = r'cascade_radius:\s*float\s*=\s*[\d\.]+\s*'
    replacement = f'cascade_radius: float = {new_cascade_radius:.1f}  '
    config_content = re.sub(pattern, replacement, config_content)
    
    # Replace peak_width_days
    pattern = r'peak_width_days:\s*int\s*=\s*\d+\s*'
    replacement = f'peak_width_days: int = {int(new_peak_width)}  '
    config_content = re.sub(pattern, replacement, config_content)
    
    # Write updated config
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    print(f"✅ Config file updated successfully!")
    
    # Save the new parameters
    new_params = {
        'p_spontaneous_female': new_p_female,
        'p_spontaneous_male': new_p_male,
        'cascade_radius': new_cascade_radius,
        'peak_width_days': new_peak_width,
        'scaling_applied': overall_scale,
        'calibration_source': str(best_params_file.name),
        'male_bout_scale_factor': male_bout_scale,
        'major_bout_scale_factor': major_bout_scale,
        'larvae_scale_factor': larvae_scale,
    }
    
    output_file = results_dir / "final_parameters_1A2.json"
    with open(output_file, 'w') as f:
        json.dump(new_params, f, indent=2)
    
    print(f"Final parameters saved to: {output_file}")
    
    return new_params

if __name__ == "__main__":
    params = update_spawning_parameters()
    print(f"\nSpawning config update complete!")
    print(f"New parameters: {params}")