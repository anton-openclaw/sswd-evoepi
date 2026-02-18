"""SSWD-EvoEpi visualization library.

Modules:
  - style: Dark theme colours and helpers
  - population: Population dynamics & demographics (10 plots)
  - disease: Disease & epidemiology (12 plots)
  - genetics: Host genetics & evolution (12 plots)
  - coevolution: Pathogen co-evolution & arms race (8 plots)
  - spatial: Spatial metapopulation & geography (8 plots)
"""

from sswd_evoepi.viz.style import (  # noqa: F401
    ACCENT_COLORS,
    DARK_BG,
    DARK_PANEL,
    DEATH_COLORS,
    GRID_COLOR,
    NODE_COLORS,
    STAGE_COLORS,
    TEXT_COLOR,
    apply_dark_theme,
    dark_figure,
    save_figure,
)

from sswd_evoepi.viz.population import (  # noqa: F401
    plot_population_trajectory,
    plot_stage_composition,
    plot_cause_of_death_breakdown,
    plot_age_size_pyramid,
    plot_population_heatmap,
    plot_recruitment_timeseries,
    plot_survival_curves,
    plot_sex_ratio_over_time,
    plot_density_dependence,
    plot_node_comparison_bars,
)

from sswd_evoepi.viz.disease import (  # noqa: F401
    COMPARTMENT_COLORS,
    plot_epidemic_curve,
    plot_vibrio_concentration,
    plot_force_of_infection_distribution,
    plot_R0_over_time,
    plot_disease_mortality_by_node,
    plot_epidemic_wave_timing,
    plot_compartment_flow_sankey,
    plot_shedding_timeseries,
    plot_disease_state_heatmap,
    plot_immunosuppression_overlap,
    plot_recovery_vs_resistance,
    plot_cfr_over_time,
)

from sswd_evoepi.viz.genetics import (  # noqa: F401
    RESISTANCE_LOW,
    RESISTANCE_MID,
    RESISTANCE_HIGH,
    plot_resistance_trajectory,
    plot_resistance_distribution,
    plot_allele_frequency_spaghetti,
    plot_additive_variance_over_time,
    plot_ef1a_dynamics,
    plot_selection_differential,
    plot_heritability_over_time,
    plot_genotype_phenotype_map,
    plot_locus_effect_size_distribution,
    plot_resistance_by_node_violin,
    plot_genetic_drift_null,
    plot_beta_init_visualization,
)

from sswd_evoepi.viz.coevolution import (  # noqa: F401
    VIRULENCE_LOW,
    VIRULENCE_HIGH,
    plot_virulence_trajectory,
    plot_coevolution_phase_portrait,
    plot_virulence_distribution_over_time,
    plot_tradeoff_curve,
    plot_R0_by_virulence,
    plot_virulence_vs_host_density,
    plot_strain_competition,
    plot_coevolution_multi_seed,
)

from sswd_evoepi.viz.spatial import (  # noqa: F401
    plot_network_map,
    plot_connectivity_heatmap,
    plot_north_south_gradient,
    plot_fjord_vs_open,
    plot_metapopulation_timeseries,
    plot_larval_flow_diagram,
    plot_spatial_epidemic_timeline,
    plot_node_fate_matrix,
)
