"""SSWD-EvoEpi visualization library.

Modules:
  - style: Dark theme colours and helpers
  - population: Population dynamics & demographics (10 plots)
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
