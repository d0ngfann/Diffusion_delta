"""
Generate discrete delta-star (δ*) heatmap plots for parameter sweep results with p1=0.01

This script processes summary CSV files from output_data/20260108_p1_0.01/ and generates
heatmap plots showing delta_star_mean across 2D parameter spaces.

Features:
- Automatic file detection and parsing
- 3 experiment types:
  * Exp1: p2 vs thrshld (beta fixed) - 16 plots
  * Exp2: beta vs thrshld (p2 fixed) - 12 plots
  * Exp3: p2 vs beta (thrshld fixed) - 12 plots
- Discrete heatmaps (no interpolation or smoothing - scientifically honest!)
- Outputs to: output_plots/20260108_p1_0.01/

Usage:
    cd scripts
    python plot_deltastar_heatmaps_p1_0.01.py
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from glob import glob
import warnings
warnings.filterwarnings('ignore')

# Get the script's directory and compute paths relative to it
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# Configuration - Modified for p1=0.01 data
DATA_DIR = os.path.join(PROJECT_ROOT, 'output_data', '20260108_p1_0.01')
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'output_plots', '20260108_p1_0.01')

# Key contour levels
CONTOUR_LEVELS = [0.3, 0.4, 0.5, 0.6, 0.7]

# High quality settings
GRID_RESOLUTION_X = 600
GRID_RESOLUTION_Y = 300
CONTOUR_LEVELS_FILLED = 101
GAUSSIAN_SIGMA = 1.0


def parse_filename(filename):
    """
    Parse filename to extract experiment information.

    Parameters:
        filename: CSV filename (e.g., '20260108_171915_sweep_summary_p2_q4.0_p2_vs_thrshld_beta1.0.csv')

    Returns:
        dict with keys: p, q, exp_type, x_param, y_param, fixed_param_name, fixed_param_value
        Returns None if parsing fails
    """
    # Extract p and q
    p_match = re.search(r'_p(\d+)_q', filename)
    q_match = re.search(r'_q([\d.]+)_', filename)

    if not p_match or not q_match:
        return None

    p = int(p_match.group(1))
    q = float(q_match.group(1))

    # Determine experiment type and parameters
    if 'p2_vs_thrshld' in filename:
        # Exp1: p2 vs thrshld (beta fixed)
        exp_type = 'p2_vs_thrshld'
        x_param = 'p2'
        y_param = 'thrshld'
        fixed_param_name = 'beta'

        beta_match = re.search(r'beta([\d.]+?)\.csv', filename)
        if not beta_match:
            return None
        fixed_param_value = float(beta_match.group(1))

    elif 'beta_vs_thrshld' in filename:
        # Exp2: beta vs thrshld (p2 fixed)
        exp_type = 'beta_vs_thrshld'
        x_param = 'beta'
        y_param = 'thrshld'
        fixed_param_name = 'p2'

        p2_match = re.search(r'p2_([\d.]+?)\.csv', filename)
        if not p2_match:
            return None
        fixed_param_value = float(p2_match.group(1))

    elif 'p2_vs_beta' in filename:
        # Exp3: p2 vs beta (thrshld fixed)
        exp_type = 'p2_vs_beta'
        x_param = 'p2'
        y_param = 'beta'
        fixed_param_name = 'thrshld'

        thrshld_match = re.search(r'thrshld([\d.]+?)\.csv', filename)
        if not thrshld_match:
            return None
        fixed_param_value = float(thrshld_match.group(1))

    else:
        return None

    return {
        'p': p,
        'q': q,
        'exp_type': exp_type,
        'x_param': x_param,
        'y_param': y_param,
        'fixed_param_name': fixed_param_name,
        'fixed_param_value': fixed_param_value
    }


def plot_case_distribution_single_file(csv_path, file_info, output_path):
    """
    Create and save a case distribution plot (no contours, just background colors).
    Used when delta_star_mean is not available.

    Parameters:
        csv_path: Path to summary CSV file
        file_info: Dictionary from parse_filename()
        output_path: Path to save the figure
    """
    # Read data
    df = pd.read_csv(csv_path)

    # Extract data
    x_param = file_info['x_param']
    y_param = file_info['y_param']

    x_data = df[x_param].values
    y_data = df[y_param].values

    # Get case percentages (compute from counts if needed)
    if 'case_A_pct' in df.columns:
        case_A_pct = df['case_A_pct'].values
        case_B_pct = df['case_B_pct'].values
        case_C_pct = df['case_C_pct'].values
    elif 'case_A_count' in df.columns:
        # Calculate percentages from counts (assuming 500 total pairs)
        total_pairs = 500
        case_A_pct = (df['case_A_count'].values / total_pairs) * 100
        case_B_pct = (df['case_B_count'].values / total_pairs) * 100
        case_C_pct = (df['case_C_count'].values / total_pairs) * 100
    else:
        case_A_pct = np.zeros(len(df))
        case_B_pct = np.zeros(len(df))
        case_C_pct = np.zeros(len(df))

    # Get unique values for axis
    x_unique = sorted(df[x_param].unique())
    y_unique = sorted(df[y_param].unique())

    # Create figure
    width = 12
    height = max(3.5, len(y_unique) * 0.4 + 2.0)
    fig, ax = plt.subplots(figsize=(width, height))

    # Plot background colored rectangles based on case distribution
    for i, row in df.iterrows():
        x_val = row[x_param]
        y_val = row[y_param]

        # Find cell boundaries
        x_idx = x_unique.index(x_val)
        y_idx = y_unique.index(y_val)

        # Calculate cell boundaries
        if x_idx < len(x_unique) - 1:
            x_left = x_val
            x_right = x_unique[x_idx + 1]
        else:
            x_left = x_unique[x_idx - 1] if x_idx > 0 else x_val - 0.01
            x_right = x_val + (x_val - x_left)

        if y_idx < len(y_unique) - 1:
            y_bottom = y_val
            y_top = y_unique[y_idx + 1]
        else:
            y_bottom = y_unique[y_idx - 1] if y_idx > 0 else y_val - 0.1
            y_top = y_val + (y_val - y_bottom)

        # Determine color based on case distribution
        # Red for Case B (Random dominates), Blue for Case A (Cpq dominates), Green for Case C (intersection)
        r = case_B_pct[i] / 100.0
        g = case_C_pct[i] / 100.0
        b = case_A_pct[i] / 100.0

        # Draw rectangle
        rect = Rectangle((x_left, y_bottom), x_right - x_left, y_top - y_bottom,
                         facecolor=(r, g, b, 0.8), edgecolor='none', linewidth=0)
        ax.add_patch(rect)

    # Set axis labels based on parameter
    if x_param == 'p2':
        ax.set_xlabel('Reinforced Adoption Probability (p2)', fontsize=12)
    elif x_param == 'beta':
        ax.set_xlabel('Time of Influence (beta)', fontsize=12)

    if y_param == 'thrshld':
        ax.set_ylabel('Threshold (thrshld)', fontsize=12)
    elif y_param == 'beta':
        ax.set_ylabel('Time of Influence (beta)', fontsize=12)

    # Create title
    p, q = file_info['p'], file_info['q']
    fixed_name = file_info['fixed_param_name']
    fixed_val = file_info['fixed_param_value']

    if fixed_name == 'beta':
        fixed_str = f'β={fixed_val:.1f}'
    elif fixed_name == 'p2':
        fixed_str = f'p2={fixed_val:.1f}'
    elif fixed_name == 'thrshld':
        fixed_str = f'thrshld={fixed_val:.1f}'

    title = f'Case Distribution (p1=0.01): {x_param} vs {y_param} ({fixed_str}) | C_{{{p},{q:.0f}/n}}'
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Set axis limits
    x_min, x_max = min(x_unique), max(x_unique)
    y_min, y_max = min(y_unique), max(y_unique)

    ax.set_xlim(x_min - 0.01, x_max + 0.01)

    if y_param in ['thrshld', 'beta']:
        y_range = y_max - y_min
        padding = max(0.3, y_range * 0.05)
        ax.set_ylim(y_min - padding, y_max + padding)
        ax.set_yticks(y_unique)

    # Add grid
    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

    # Add legend
    textstr = 'Color: Red=Random dominates, Blue=Cpq dominates, Green=Intersection exists\n(No δ* data available for contours)'
    props = dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray')
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved (case distribution only): {os.path.basename(output_path)}")


def plot_contour_single_file(csv_path, file_info, output_path):
    """
    Create and save a discrete heatmap for one CSV file.

    Parameters:
        csv_path: Path to summary CSV file
        file_info: Dictionary from parse_filename()
        output_path: Path to save the figure
    """
    # Read data
    df = pd.read_csv(csv_path)

    # Check if we have delta_star_mean data
    df_with_delta = df[df['delta_star_mean'].notna()]

    if len(df_with_delta) == 0:
        print(f"    ⚠ No delta_star_mean data, plotting case distribution instead...")
        plot_case_distribution_single_file(csv_path, file_info, output_path)
        return

    # Use data with delta_star_mean
    df = df_with_delta

    # Extract data
    x_param = file_info['x_param']
    y_param = file_info['y_param']

    # Get unique sorted values for both axes
    x_unique = np.array(sorted(df[x_param].unique()))
    y_unique = np.array(sorted(df[y_param].unique()))

    # Create 2D grid for heatmap
    z_grid = np.full((len(y_unique), len(x_unique)), np.nan)

    # Fill grid with actual data
    for _, row in df.iterrows():
        x_idx = np.where(x_unique == row[x_param])[0][0]
        y_idx = np.where(y_unique == row[y_param])[0][0]
        z_grid[y_idx, x_idx] = row['delta_star_mean']

    # Create figure
    width = 12
    height = max(3.5, len(y_unique) * 0.4 + 2.0)
    fig, ax = plt.subplots(figsize=(width, height))

    # Create discrete heatmap using pcolormesh
    # Use edges for proper alignment
    x_edges = np.concatenate([
        [x_unique[0] - (x_unique[1] - x_unique[0])/2],
        (x_unique[:-1] + x_unique[1:]) / 2,
        [x_unique[-1] + (x_unique[-1] - x_unique[-2])/2]
    ])
    y_edges = np.concatenate([
        [y_unique[0] - (y_unique[1] - y_unique[0])/2],
        (y_unique[:-1] + y_unique[1:]) / 2,
        [y_unique[-1] + (y_unique[-1] - y_unique[-2])/2]
    ])

    # Plot heatmap
    mesh = ax.pcolormesh(x_edges, y_edges, z_grid,
                         cmap='coolwarm',
                         vmin=0, vmax=1,
                         shading='flat',
                         edgecolors='face',
                         linewidth=0)

    # Add colorbar
    cbar = plt.colorbar(mesh, ax=ax, label='δ* (delta_star_mean)')
    cbar.ax.tick_params(labelsize=10)

    # Add reference lines on colorbar at key thresholds
    for level in CONTOUR_LEVELS:
        cbar.ax.axhline(level, color='black', linewidth=1.5, alpha=0.6, linestyle='--')

    # Set axis labels based on parameter
    if x_param == 'p2':
        ax.set_xlabel('Reinforced Adoption Probability (p2)', fontsize=12)
    elif x_param == 'beta':
        ax.set_xlabel('Time of Influence (beta)', fontsize=12)

    if y_param == 'thrshld':
        ax.set_ylabel('Threshold (thrshld)', fontsize=12)
    elif y_param == 'beta':
        ax.set_ylabel('Time of Influence (beta)', fontsize=12)

    # Create title
    p, q = file_info['p'], file_info['q']
    fixed_name = file_info['fixed_param_name']
    fixed_val = file_info['fixed_param_value']

    if fixed_name == 'beta':
        fixed_str = f'β={fixed_val:.1f}'
    elif fixed_name == 'p2':
        fixed_str = f'p2={fixed_val:.1f}'
    elif fixed_name == 'thrshld':
        fixed_str = f'thrshld={fixed_val:.1f}'

    title = f'δ* Heatmap (p1=0.01): {x_param} vs {y_param} ({fixed_str}) | C_{{{p},{q:.0f}/n}}'
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Set axis limits
    x_min, x_max = x_unique.min(), x_unique.max()
    y_min, y_max = y_unique.min(), y_unique.max()

    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])

    # Set y-axis ticks for discrete parameters
    if y_param in ['thrshld', 'beta']:
        # Show subset of ticks to avoid overcrowding
        if len(y_unique) > 20:
            tick_indices = np.linspace(0, len(y_unique)-1, 11, dtype=int)
            ax.set_yticks(y_unique[tick_indices])
        else:
            ax.set_yticks(y_unique)

    # Add grid
    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

    # Add annotation
    textstr = f'p1=0.01 | Discrete Heatmap: {len(x_unique)}×{len(y_unique)} parameter combinations\nDashed lines on colorbar: δ* = 0.3, 0.4, 0.5, 0.6, 0.7'
    props = dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray')
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: {os.path.basename(output_path)}")


def main():
    """
    Main function to process all CSV files and generate discrete heatmap plots.
    """
    print("="*70)
    print("DELTA-STAR HEATMAP GENERATOR (p1=0.01)")
    print("Settings:")
    print(f"  - Base adoption probability: p1=0.01 (NEW)")
    print(f"  - Plot type: Discrete heatmap (no interpolation)")
    print(f"  - Colorbar reference levels: {CONTOUR_LEVELS}")
    print(f"  - Output DPI: 200")
    print("="*70)

    # Find all summary CSV files
    pattern = os.path.join(DATA_DIR, '*sweep_summary*.csv')
    csv_files = sorted(glob(pattern))

    if not csv_files:
        print(f"\n✗ No summary CSV files found in: {DATA_DIR}")
        return

    print(f"\nFound {len(csv_files)} summary CSV files")
    print(f"Data directory: {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Track statistics
    exp_counts = {'p2_vs_thrshld': 0, 'beta_vs_thrshld': 0, 'p2_vs_beta': 0}
    skipped = 0

    # Process each file
    for csv_path in csv_files:
        filename = os.path.basename(csv_path)

        # Parse filename
        file_info = parse_filename(filename)

        if file_info is None:
            print(f"✗ Could not parse: {filename}")
            skipped += 1
            continue

        # Determine output directory based on experiment type
        exp_type = file_info['exp_type']
        exp_folder = f"exp_{exp_type}"
        output_subdir = os.path.join(OUTPUT_DIR, exp_folder)
        os.makedirs(output_subdir, exist_ok=True)

        # Create output filename
        p = file_info['p']
        q = int(file_info['q'])
        fixed_name = file_info['fixed_param_name']
        fixed_val = file_info['fixed_param_value']

        output_filename = f"p{p}q{q}_{fixed_name}{fixed_val:.1f}.png"
        output_path = os.path.join(output_subdir, output_filename)

        # Generate plot
        print(f"\nProcessing: {filename}")
        print(f"  Network: C_{{{p},{file_info['q']:.0f}/n}}")
        print(f"  Type: {exp_type} ({fixed_name}={fixed_val:.1f})")

        plot_contour_single_file(csv_path, file_info, output_path)

        exp_counts[exp_type] += 1

    # Print summary
    print("\n" + "="*70)
    print("DELTA-STAR HEATMAP GENERATION COMPLETE (p1=0.01)!")
    print("="*70)
    print(f"\nGenerated plots:")
    print(f"  - Exp1 (p2 vs thrshld): {exp_counts['p2_vs_thrshld']} plots")
    print(f"  - Exp2 (beta vs thrshld): {exp_counts['beta_vs_thrshld']} plots")
    print(f"  - Exp3 (p2 vs beta): {exp_counts['p2_vs_beta']} plots")
    print(f"  - Total: {sum(exp_counts.values())} plots")
    if skipped > 0:
        print(f"  - Skipped: {skipped} files")

    print(f"\nOutput saved to: {OUTPUT_DIR}/")
    print("  - exp_p2_vs_thrshld/")
    print("  - exp_beta_vs_thrshld/")
    print("  - exp_p2_vs_beta/")
    print()


if __name__ == "__main__":
    main()
