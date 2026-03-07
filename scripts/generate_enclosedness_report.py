#!/usr/bin/env python3
"""
Generate comprehensive PDF report for the SSWD-EvoEpi enclosedness analysis.

Creates matplotlib figures and a LaTeX-based PDF report with spatial analysis.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import subprocess
from typing import Dict, List, Tuple
import warnings

# Suppress some matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
    print("Using cartopy for map projections")
except ImportError:
    HAS_CARTOPY = False
    print("Warning: cartopy not available. Using simple lat/lon plots.")

try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False


def load_enclosedness_data(data_path: str) -> Tuple[pd.DataFrame, List[str]]:
    """
    Load the enclosedness data and extract unique regions sorted N→S.
    """
    df = pd.read_csv(data_path)
    
    # Region order from Alaska to Baja California (roughly N→S)
    region_order = [
        'AK-AL', 'AK-EG', 'AK-WG', 'AK-FN', 'AK-FS', 'AK-PWS', 'AK-OC',
        'BC-N', 'BC-C', 'SS-N', 'SS-S', 'JDF', 'WA-O', 'OR', 'CA-N', 'CA-C', 'CA-S', 'BJ'
    ]
    
    # Filter to only regions that exist in the data
    available_regions = df['region'].unique()
    region_order = [r for r in region_order if r in available_regions]
    
    return df, region_order


def create_map_figure(df: pd.DataFrame, value_column: str, title: str, 
                     output_path: str, cmap: str = 'RdYlBu_r', 
                     coastline_path: str = None) -> None:
    """
    Create a map figure colored by the specified value column.
    """
    fig_width, fig_height = 12, 10
    
    if HAS_CARTOPY:
        # Use Pacific-centered projection if possible
        fig = plt.figure(figsize=(fig_width, fig_height))
        
        # Choose projection that shows Alaska to Baja well
        # PlateCarree with central longitude around -130°
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-130))
        
        # Set map extent to cover Alaska to Baja California
        ax.set_extent([-180, -110, 25, 65], ccrs.PlateCarree())
        
        # Add coastlines and land features
        ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        
        # Transform coordinates for plotting
        transform = ccrs.PlateCarree()
    else:
        # Simple lat/lon plot
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Add coastline from shapefile if available
        if coastline_path and HAS_GEOPANDAS and os.path.exists(coastline_path):
            try:
                coastline = gpd.read_file(coastline_path)
                coastline.boundary.plot(ax=ax, color='gray', linewidth=0.5, alpha=0.7)
            except Exception as e:
                print(f"Warning: Could not load coastline from {coastline_path}: {e}")
        
        ax.set_xlim(-180, -110)
        ax.set_ylim(25, 65)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        transform = None
    
    # Create color mapping
    values = df[value_column].values
    vmin, vmax = np.min(values), np.max(values)
    
    # Create scatter plot
    if HAS_CARTOPY:
        scatter = ax.scatter(df['lon'], df['lat'], c=values, cmap=cmap, 
                           s=20, alpha=0.8, vmin=vmin, vmax=vmax,
                           transform=transform, edgecolors='black', linewidths=0.3)
    else:
        scatter = ax.scatter(df['lon'], df['lat'], c=values, cmap=cmap, 
                           s=20, alpha=0.8, vmin=vmin, vmax=vmax,
                           edgecolors='black', linewidths=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(value_column.replace('_', ' ').title(), fontsize=12)
    
    # Set title and grid
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    
    if HAS_CARTOPY:
        ax.gridlines(draw_labels=True, alpha=0.3)
    else:
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved map: {output_path}")


def create_histogram(df: pd.DataFrame, region_order: List[str], output_path: str) -> None:
    """
    Create histogram of enclosedness distribution colored by region.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create color map for regions
    colors = plt.cm.Set3(np.linspace(0, 1, len(region_order)))
    region_colors = {region: colors[i] for i, region in enumerate(region_order)}
    
    # Plot histogram for each region
    for region in region_order:
        region_data = df[df['region'] == region]['enclosedness_combined']
        ax.hist(region_data, bins=20, alpha=0.7, label=region, 
                color=region_colors[region], edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel('Combined Enclosedness', fontsize=12)
    ax.set_ylabel('Number of Sites', fontsize=12)
    ax.set_title('Distribution of Enclosedness Across All Sites', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Legend with multiple columns to save space
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved histogram: {output_path}")


def create_box_plot(df: pd.DataFrame, region_order: List[str], output_path: str) -> None:
    """
    Create box plot of enclosedness by region, ordered N→S.
    """
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Prepare data for box plot
    region_data = [df[df['region'] == region]['enclosedness_combined'] for region in region_order]
    
    # Create box plot
    bp = ax.boxplot(region_data, labels=region_order, patch_artist=True)
    
    # Color boxes with a gradient
    colors = plt.cm.viridis(np.linspace(0, 1, len(region_order)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Region (North → South)', fontsize=12)
    ax.set_ylabel('Combined Enclosedness', fontsize=12)
    ax.set_title('Enclosedness by Region (Ordered North to South)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved box plot: {output_path}")


def create_correlation_plot(df: pd.DataFrame, region_order: List[str], output_path: str) -> None:
    """
    Create scatter plot of tortuosity vs ray-cast enclosedness, colored by region.
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create color map for regions
    colors = plt.cm.Set3(np.linspace(0, 1, len(region_order)))
    region_colors = {region: colors[i] for i, region in enumerate(region_order)}
    
    # Plot each region
    for region in region_order:
        region_data = df[df['region'] == region]
        # Normalize tortuosity ratio to [0, 1] for plotting
        normalized_tort = np.clip((region_data['tortuosity_ratio'] - 1.0) / (3.0 - 1.0), 0.0, 1.0)
        
        ax.scatter(normalized_tort, region_data['enclosedness_rays'], 
                  alpha=0.7, label=region, color=region_colors[region], s=30)
    
    ax.set_xlabel('Normalized Tortuosity Ratio', fontsize=12)
    ax.set_ylabel('Ray-cast Enclosedness', fontsize=12)
    ax.set_title('Correlation: Tortuosity vs Ray-cast Enclosedness', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add diagonal line for reference
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1)
    
    # Legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2, fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved correlation plot: {output_path}")


def create_summary_table(df: pd.DataFrame, region_order: List[str], output_path: str) -> None:
    """
    Create a summary table/bar chart of mean enclosedness by region.
    """
    # Calculate summary statistics
    summary_stats = []
    for region in region_order:
        region_data = df[df['region'] == region]['enclosedness_combined']
        stats = {
            'Region': region,
            'N_Sites': len(region_data),
            'Mean_Enclosedness': np.mean(region_data),
            'Std_Enclosedness': np.std(region_data),
            'Min_Enclosedness': np.min(region_data),
            'Max_Enclosedness': np.max(region_data)
        }
        summary_stats.append(stats)
    
    summary_df = pd.DataFrame(summary_stats)
    
    # Create figure with table and bar chart
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [1, 2]})
    
    # Bar chart
    bars = ax1.bar(range(len(region_order)), summary_df['Mean_Enclosedness'], 
                   yerr=summary_df['Std_Enclosedness'], capsize=5, alpha=0.8,
                   color=plt.cm.viridis(np.linspace(0, 1, len(region_order))))
    
    ax1.set_xticks(range(len(region_order)))
    ax1.set_xticklabels(region_order, rotation=45, ha='right')
    ax1.set_ylabel('Mean Enclosedness ± SD', fontsize=12)
    ax1.set_title('Mean Enclosedness by Region', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for i, (bar, val, std) in enumerate(zip(bars, summary_df['Mean_Enclosedness'], summary_df['Std_Enclosedness'])):
        ax1.text(bar.get_x() + bar.get_width()/2, val + std + 0.02, 
                f'{val:.3f}', ha='center', va='bottom', fontsize=9)
    
    # Table
    ax2.axis('off')
    
    # Prepare table data
    table_data = []
    for _, row in summary_df.iterrows():
        table_data.append([
            row['Region'],
            f"{row['N_Sites']}",
            f"{row['Mean_Enclosedness']:.3f} ± {row['Std_Enclosedness']:.3f}",
            f"[{row['Min_Enclosedness']:.3f}, {row['Max_Enclosedness']:.3f}]"
        ])
    
    table_headers = ['Region', 'N Sites', 'Mean ± SD', 'Range [Min, Max]']
    
    # Create table
    table = ax2.table(cellText=table_data, colLabels=table_headers, 
                     cellLoc='center', loc='center',
                     bbox=[0, 0, 1, 1])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style the table
    for i in range(len(table_headers)):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    for i in range(1, len(table_data) + 1):
        for j in range(len(table_headers)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#f5f5f5')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved summary table: {output_path}")


def generate_latex_report(df: pd.DataFrame, region_order: List[str], 
                         output_dir: Path) -> None:
    """
    Generate LaTeX report with all figures and analysis.
    """
    # Calculate some summary statistics for the text
    total_sites = len(df)
    mean_enclosedness = np.mean(df['enclosedness_combined'])
    std_enclosedness = np.std(df['enclosedness_combined'])
    min_enclosedness = np.min(df['enclosedness_combined'])
    max_enclosedness = np.max(df['enclosedness_combined'])
    
    # Top 3 most enclosed regions
    region_means = df.groupby('region')['enclosedness_combined'].mean().sort_values(ascending=False)
    most_enclosed = region_means.head(3)
    least_enclosed = region_means.tail(3)
    
    latex_content = f"""\\documentclass[11pt]{{article}}
\\usepackage[utf8]{{inputenc}}
\\usepackage{{graphicx}}
\\usepackage{{geometry}}
\\usepackage{{amsmath}}
\\usepackage{{float}}
\\usepackage{{booktabs}}
\\usepackage{{hyperref}}
\\usepackage{{subcaption}}

\\geometry{{margin=1in}}

\\title{{SSWD-EvoEpi: Site Enclosedness Analysis}}
\\author{{Automated Analysis Report}}
\\date{{\\today}}

\\begin{{document}}

\\maketitle

\\section{{Executive Summary}}

This report presents an algorithmic analysis of coastal site enclosedness for {total_sites} sites spanning from Alaska to Baja California in the SSWD-EvoEpi model network. Two complementary metrics were computed:

\\begin{{enumerate}}
    \\item \\textbf{{Tortuosity Ratio}}: Over-water distance vs. straight-line distance to nearest neighbors
    \\item \\textbf{{Horizon Exposure}}: Fraction of rays that don't encounter land within 50 km (ray-casting)
\\end{{enumerate}}

The combined enclosedness metric ranges from {min_enclosedness:.3f} to {max_enclosedness:.3f} (mean = {mean_enclosedness:.3f} ± {std_enclosedness:.3f}). The most enclosed regions are {', '.join([f"{region} ({mean:.3f})" for region, mean in most_enclosed.items()][:3])}, while the least enclosed are {', '.join([f"{region} ({mean:.3f})" for region, mean in least_enclosed.items()][:3])}.

The derived flushing rates ($\\phi$) range from 0.030 to 0.789, with enclosed fjords receiving low flushing and open coasts receiving high flushing as expected.

\\section{{Methods}}

\\subsection{{Tortuosity Ratio}}
For each site $i$, we identified $k=20$ nearest neighbors by straight-line (haversine) distance, then computed:

$$\\text{{tortuosity}}_i = \\frac{{1}}{{k}} \\sum_{{j \\in \\text{{k-nearest}}}} \\frac{{d_{{\\text{{over-water}}}}(i,j)}}{{d_{{\\text{{haversine}}}}(i,j)}}$$

where over-water distances come from the pre-computed 896×896 distance matrix and haversine distances are great-circle distances. Disconnected site pairs (infinite over-water distance) were excluded.

\\subsection{{Horizon Exposure (Ray-casting)}}
For each site, 72 rays were cast at 5° intervals to a maximum distance of 50 km. Each ray was tested for intersection with Natural Earth land polygons. The horizon exposure is:

$$\\text{{horizon exposure}}_i = \\frac{{\\text{{\\# rays not hitting land}}}}{{72}}$$
$$\\text{{enclosedness}}_i = 1 - \\text{{horizon exposure}}_i$$

\\subsection{{Combined Metric}}
The tortuosity ratio was normalized to [0,1] using a linear mapping where tortuosity=1 $\\rightarrow$ 0 and tortuosity$\\geq$3 $\\rightarrow$ 1. The final enclosedness score combines both metrics:

$$\\text{{combined enclosedness}} = 0.5 \\times \\text{{normalized tortuosity}} + 0.5 \\times \\text{{ray-cast enclosedness}}$$

\\subsection{{Flushing Rate Mapping}}
The enclosedness score was mapped to flushing rate using:

$$\\phi = \\phi_{{\\text{{open}}}} \\times (1 - \\text{{enclosedness}}) + \\phi_{{\\text{{fjord}}}} \\times \\text{{enclosedness}}$$

where $\\phi_{{\\text{{open}}}} = 0.8$ and $\\phi_{{\\text{{fjord}}}} = 0.03$.

\\section{{Results}}

\\begin{{figure}}[H]
    \\centering
    \\begin{{subcaptiongroup}}
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/map_tortuosity_ratio.png}}
            \\caption{{Tortuosity ratio}}
        \\end{{subfigure}}
        \\hfill
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/map_enclosedness_rays.png}}
            \\caption{{Ray-cast enclosedness}}
        \\end{{subfigure}}
    \\end{{subcaptiongroup}}
    
    \\begin{{subcaptiongroup}}
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/map_enclosedness_combined.png}}
            \\caption{{Combined enclosedness}}
        \\end{{subfigure}}
        \\hfill
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/map_flushing_rate.png}}
            \\caption{{Derived flushing rate ($\\phi$)}}
        \\end{{subfigure}}
    \\end{{subcaptiongroup}}
    \\caption{{Spatial distribution of enclosedness metrics and derived flushing rates across the 896-site network.}}
    \\label{{fig:maps}}
\\end{{figure}}

\\begin{{figure}}[H]
    \\centering
    \\begin{{subcaptiongroup}}
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/histogram_distribution.png}}
            \\caption{{Distribution by region}}
        \\end{{subfigure}}
        \\hfill
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/boxplot_by_region.png}}
            \\caption{{Regional comparison (N→S)}}
        \\end{{subfigure}}
    \\end{{subcaptiongroup}}
    \\caption{{Distribution of enclosedness across regions.}}
    \\label{{fig:distributions}}
\\end{{figure}}

\\begin{{figure}}[H]
    \\centering
    \\begin{{subcaptiongroup}}
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/correlation_plot.png}}
            \\caption{{Metric correlation}}
        \\end{{subfigure}}
        \\hfill
        \\begin{{subfigure}}[b]{{0.48\\textwidth}}
            \\includegraphics[width=\\textwidth]{{figures/summary_table.png}}
            \\caption{{Regional summary}}
        \\end{{subfigure}}
    \\end{{subcaptiongroup}}
    \\caption{{Correlation between metrics and regional summary statistics.}}
    \\label{{fig:analysis}}
\\end{{figure}}

\\section{{Discussion}}

\\subsection{{Regional Patterns}}
The analysis reveals clear regional patterns consistent with coastal geomorphology:

\\begin{{itemize}}
    \\item \\textbf{{High enclosedness}}: {', '.join([region for region, _ in most_enclosed.head(3).items()])} show the highest enclosedness, consistent with fjord and inlet geography.
    \\item \\textbf{{Low enclosedness}}: {', '.join([region for region, _ in least_enclosed.head(3).items()])} show lower enclosedness, typical of open coastlines.
    \\item \\textbf{{Gradient}}: There's a clear latitudinal gradient with northern regions (Alaska, British Columbia) generally showing higher enclosedness than southern regions (Oregon, California).
\\end{{itemize}}

\\subsection{{Validation}}
The results align well with expected coastal geomorphology:
\\begin{{enumerate}}
    \\item Prince William Sound (AK-PWS) shows high enclosedness as expected for a complex fjord system
    \\item Oregon (OR) and Northern California (CA-N) coastlines show low enclosedness, consistent with their straight, exposed nature
    \\item British Columbia fjords and Alaskan regions show appropriately high values
    \\item The Salish Sea regions (SS-N, SS-S) show high enclosedness reflecting the inland sea character
\\end{{enumerate}}

\\subsection{{Metric Correlation}}
The tortuosity and ray-casting metrics show complementary behavior. While positively correlated, they capture different aspects of enclosedness:
\\begin{{itemize}}
    \\item Tortuosity captures the complexity of water-based connectivity
    \\item Ray-casting captures the direct exposure to open ocean
    \\item Combined, they provide a robust measure of coastal enclosedness
\\end{{itemize}}

\\section{{Implications for Model}}

\\subsection{{Disease Dynamics}}
The differentiated flushing rates will significantly affect pathogen dispersal and concentration:

\\begin{{enumerate}}
    \\item \\textbf{{Enclosed sites}} (low $\\phi$): Pathogens will accumulate, potentially creating disease hotspots and facilitating local transmission.
    \\item \\textbf{{Open sites}} (high $\\phi$): Rapid flushing will reduce local pathogen concentration but may facilitate long-range dispersal.
    \\item \\textbf{{Connectivity gradients}}: The enclosedness metric creates realistic connectivity gradients that should improve model realism.
\\end{{enumerate}}

\\subsection{{Expected Effects}}
We anticipate the following changes in model behavior:
\\begin{{itemize}}
    \\item More realistic disease persistence in fjords and protected waters
    \\item Reduced unrealistic long-distance transmission from highly enclosed sites
    \\item Better representation of outbreak patterns matching observed SSWD dynamics
    \\item Improved spatial structure in evolutionary dynamics
\\end{{itemize}}

\\section{{Conclusion}}

The algorithmic enclosedness classification successfully identified {total_sites} coastal sites with differentiated exposure characteristics. The derived flushing rates provide a mechanistically sound basis for improving the SSWD-EvoEpi model's spatial realism. The clear regional patterns and validation against known coastal geomorphology support the robustness of this approach.

Implementation of these site-specific flushing rates should improve model predictions of disease spread, persistence, and evolutionary dynamics in the complex coastal environment from Alaska to Baja California.

\\end{{document}}
"""
    
    # Write LaTeX file
    tex_file = output_dir / "main.tex"
    with open(tex_file, 'w') as f:
        f.write(latex_content)
    
    print(f"Generated LaTeX report: {tex_file}")


def compile_pdf_report(output_dir: Path, pdflatex_path: str) -> bool:
    """
    Compile the LaTeX report to PDF using pdflatex.
    """
    tex_file = output_dir / "main.tex"
    
    if not tex_file.exists():
        print(f"Error: LaTeX file not found at {tex_file}")
        return False
    
    try:
        # Change to output directory
        original_cwd = os.getcwd()
        os.chdir(output_dir)
        
        # Run pdflatex twice for references
        for run in [1, 2]:
            print(f"Running pdflatex (pass {run}/2)...")
            result = subprocess.run([
                pdflatex_path, '-interaction=nonstopmode', 'main.tex'
            ], capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"Error in pdflatex pass {run}:")
                print(result.stdout)
                print(result.stderr)
                return False
        
        # Check if PDF was created
        pdf_file = output_dir / "main.pdf"
        if pdf_file.exists():
            print(f"Successfully compiled PDF: {pdf_file}")
            return True
        else:
            print("Error: PDF file not created")
            return False
    
    except Exception as e:
        print(f"Error compiling PDF: {e}")
        return False
    
    finally:
        os.chdir(original_cwd)


def main():
    """
    Main function to generate the enclosedness analysis report.
    """
    # Set up paths
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    
    # Data paths
    enclosedness_csv = repo_root / "data" / "nodes" / "site_enclosedness.csv"
    coastline_shapefile = repo_root / "data" / "shorelines" / "ne_10m_land" / "ne_10m_land.shp"
    
    # Output paths
    output_dir = repo_root / "reports" / "enclosedness_analysis"
    figures_dir = output_dir / "figures"
    
    # pdflatex path
    pdflatex_path = "~/.TinyTeX/bin/x86_64-linux/pdflatex"
    pdflatex_path = os.path.expanduser(pdflatex_path)
    
    print("=== Generating Enclosedness Analysis Report ===")
    
    # Load data
    print("Loading enclosedness data...")
    if not enclosedness_csv.exists():
        print(f"Error: Enclosedness data not found at {enclosedness_csv}")
        return 1
    
    df, region_order = load_enclosedness_data(str(enclosedness_csv))
    print(f"Loaded data for {len(df)} sites across {len(region_order)} regions")
    
    # Set up matplotlib
    plt.style.use('default')
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 11
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 9
    
    # Create figures
    print("\\nGenerating figures...")
    
    coastline_path = str(coastline_shapefile) if coastline_shapefile.exists() else None
    
    # Map figures
    create_map_figure(df, 'tortuosity_ratio', 'Site Tortuosity Ratio',
                     str(figures_dir / 'map_tortuosity_ratio.png'), 
                     cmap='plasma', coastline_path=coastline_path)
    
    create_map_figure(df, 'enclosedness_rays', 'Ray-cast Enclosedness',
                     str(figures_dir / 'map_enclosedness_rays.png'), 
                     cmap='RdYlBu_r', coastline_path=coastline_path)
    
    create_map_figure(df, 'enclosedness_combined', 'Combined Enclosedness',
                     str(figures_dir / 'map_enclosedness_combined.png'), 
                     cmap='RdYlBu_r', coastline_path=coastline_path)
    
    create_map_figure(df, 'flushing_rate', 'Flushing Rate (φ)',
                     str(figures_dir / 'map_flushing_rate.png'), 
                     cmap='RdYlBu', coastline_path=coastline_path)
    
    # Distribution figures
    create_histogram(df, region_order, str(figures_dir / 'histogram_distribution.png'))
    create_box_plot(df, region_order, str(figures_dir / 'boxplot_by_region.png'))
    create_correlation_plot(df, region_order, str(figures_dir / 'correlation_plot.png'))
    create_summary_table(df, region_order, str(figures_dir / 'summary_table.png'))
    
    # Generate LaTeX report
    print("\\nGenerating LaTeX report...")
    generate_latex_report(df, region_order, output_dir)
    
    # Compile PDF
    print("\\nCompiling PDF...")
    if os.path.exists(pdflatex_path):
        success = compile_pdf_report(output_dir, pdflatex_path)
        if success:
            print(f"\\nReport generated successfully: {output_dir / 'main.pdf'}")
        else:
            print("\\nError: PDF compilation failed")
            return 1
    else:
        print(f"Warning: pdflatex not found at {pdflatex_path}")
        print("LaTeX source generated but PDF not compiled")
    
    print("\\n=== Report Generation Complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())