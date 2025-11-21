"""
Parameter Sweep Script for Delta-Star Analysis

This script performs a systematic parameter sweep to find delta star across
multiple combinations of network and diffusion parameters.

The sweep explores:
- Network degree (k)
- Social reinforcement threshold (thrshld)
- Base and reinforced adoption probabilities (p1, p2)
- Time of influence (beta)

Results are saved to CSV with summary statistics and visualizations.
"""

import os
import sys
import time
import itertools
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

# IMPORTANT: Set matplotlib to non-interactive mode to prevent plot windows from appearing
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    print("Warning: seaborn not installed - some visualizations may be simplified")

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    # Create a simple progress indicator if tqdm is not available
    class tqdm:
        def __init__(self, iterable, desc="", unit=""):
            self.iterable = iterable
            self.desc = desc
        def __iter__(self):
            return iter(self.iterable)
        @staticmethod
        def write(msg):
            print(msg)

# Import the main analysis function
from DH_delta_star_analysis import main as run_delta_star_analysis


def generate_parameter_combinations():
    """
    Generate all valid parameter combinations for the sweep.

    Returns:
        list of dicts: Each dict contains one parameter combination
    """
    # Define parameter grids
    k_values = [4, 8, 12]
    thrshld_values = [2, 3, 4]
    p1_values = [0.05, 0.1, 0.15]
    p2_values = list(np.linspace(0.3, 1.0, 30)) 
    beta_values = [1, 3, 5, 7]

    # Generate all combinations
    all_combos = list(itertools.product(
        k_values, thrshld_values, p1_values, p2_values, beta_values
    ))

    # Filter: only keep combinations where p1 <= p2
    valid_combos = []
    combination_id = 1
    for k, thrshld, p1, p2, beta in all_combos:
        if p1 <= p2:
            valid_combos.append({
                'combination_id': combination_id,
                'k': k,
                'thrshld': thrshld,
                'p1': p1,
                'p2': p2,
                'beta': beta,
                'seeds': thrshld  # seeds always equals thrshld
            })
            combination_id += 1

    return valid_combos


def estimate_runtime(n_combinations, seconds_per_combo=5):
    """
    Estimate total runtime for the sweep.

    Parameters:
        n_combinations: int
            Number of parameter combinations to test
        seconds_per_combo: float
            Estimated seconds per combination (default: 5)

    Returns:
        str: Formatted time estimate
    """
    total_seconds = n_combinations * seconds_per_combo
    return str(timedelta(seconds=int(total_seconds)))


def run_single_combination(params, verbose=False):
    """
    Run delta-star analysis for a single parameter combination.

    Parameters:
        params: dict
            Parameter dictionary with keys: k, thrshld, p1, p2, beta, seeds
        verbose: bool
            If True, print detailed output from the analysis

    Returns:
        dict: Results dictionary with delta_star, num_intersections, etc.
    """
    try:
        # Run the analysis (with output suppressed if not verbose)
        delta_star, results_summary, all_intersections = run_delta_star_analysis(
            k=params['k'],
            n=None,  # Use default: k * 250
            thrshld=params['thrshld'],
            p1=params['p1'],
            p2=params['p2'],
            beta=params['beta'],
            trials=50,  # Standard number of trials
            seeds=params['seeds'],
            delta_values=None,  
            verbose=verbose
        )

        # Extract dominant network info from results_summary (take first row since values are repeated)
        dominant_info = {
            'delta_star_count': results_summary['delta_star_count'].iloc[0],
            'delta_star_values': results_summary['delta_star_values'].iloc[0],
            'dominant_below': results_summary['dominant_below'].iloc[0],
            'dominant_above': results_summary['dominant_above'].iloc[0],
            'always_dominant': results_summary['always_dominant'].iloc[0]
        }

        # Extract results
        result = {
            'primary_delta_star': delta_star,
            'num_intersections': len(all_intersections) if all_intersections else 0,
            'status': 'not_found' if delta_star is None else ('multiple' if len(all_intersections) > 1 else 'found'),
            # Add dominant network info
            'delta_star_count': dominant_info['delta_star_count'],
            'delta_star_values': dominant_info['delta_star_values'],
            'dominant_below': dominant_info['dominant_below'],
            'dominant_above': dominant_info['dominant_above'],
            'always_dominant': dominant_info['always_dominant']
        }

        # Add individual intersection values (up to 5)
        for i in range(5):
            if all_intersections and i < len(all_intersections):
                result[f'intersection_{i+1}'] = all_intersections[i]['delta_star']
            else:
                result[f'intersection_{i+1}'] = None

        return result

    except Exception as e:
        # Log error and return error status
        print(f"\n✗ Error in combination {params['combination_id']}: {str(e)}")
        return {
            'primary_delta_star': None,
            'num_intersections': 0,
            'status': 'error',
            'delta_star_count': 0,
            'delta_star_values': None,
            'dominant_below': None,
            'dominant_above': None,
            'always_dominant': None,
            'intersection_1': None,
            'intersection_2': None,
            'intersection_3': None,
            'intersection_4': None,
            'intersection_5': None,
            'error_message': str(e)
        }


def run_parameter_sweep(combinations, output_dir='../output_data', n_workers=8):
    """
    Run the complete parameter sweep using parallel processing.

    Parameters:
        combinations: list of dicts
            Parameter combinations to test
        output_dir: str
            Directory to save results
        n_workers: int
            Number of parallel workers (CPU cores) to use (default: 8)

    Returns:
        pd.DataFrame: Results dataframe
    """
    n_combos = len(combinations)
    results = []

    print(f"\nStarting parameter sweep with {n_workers} parallel workers...")
    print(f"Total combinations to test: {n_combos}")
    print(f"Estimated sequential time: {estimate_runtime(n_combos)}")
    print(f"Expected parallel speedup: ~{n_workers}x faster")
    print()

    start_time = time.time()

    # Use ProcessPoolExecutor for parallel execution
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all jobs
        future_to_params = {
            executor.submit(run_single_combination, params, False): params
            for params in combinations
        }

        # Process results as they complete
        completed = 0
        if HAS_TQDM:
            progress_bar = tqdm(total=n_combos, desc="Progress", unit="combo")

        for future in as_completed(future_to_params):
            params = future_to_params[future]
            completed += 1

            try:
                result = future.result()

                # Merge parameters and results
                full_result = {**params, **result}
                results.append(full_result)

                # Calculate timing
                elapsed = time.time() - start_time
                remaining = (elapsed / completed) * (n_combos - completed)

                # Print status for this combination
                status_symbol = "✓" if result['status'] == 'found' else ("✗" if result['status'] == 'error' else "○")
                delta_str = f"{result['primary_delta_star']:.4f}" if result['primary_delta_star'] is not None else "None"

                msg = (f"{status_symbol} [{completed}/{n_combos}] k={params['k']}, i={params['thrshld']}, "
                       f"p1={params['p1']}, p2={params['p2']}, β={params['beta']} → δ*={delta_str} | "
                       f"Elapsed: {timedelta(seconds=int(elapsed))} | "
                       f"Remaining: {timedelta(seconds=int(remaining))}")

                if HAS_TQDM:
                    tqdm.write(msg)
                    progress_bar.update(1)
                else:
                    print(msg)

            except Exception as e:
                print(f"\n✗ Unexpected error in combination {params['combination_id']}: {str(e)}")
                # Add error result
                error_result = {
                    **params,
                    'primary_delta_star': None,
                    'num_intersections': 0,
                    'status': 'error',
                    'delta_star_count': 0,
                    'delta_star_values': None,
                    'dominant_below': None,
                    'dominant_above': None,
                    'always_dominant': None,
                    'error_message': str(e)
                }
                results.append(error_result)

                if HAS_TQDM:
                    progress_bar.update(1)

        if HAS_TQDM:
            progress_bar.close()

    total_time = time.time() - start_time
    print(f"\n✓ Sweep complete! Total time: {timedelta(seconds=int(total_time))}")
    print(f"  Average time per combination: {total_time/n_combos:.2f} seconds")
    print(f"  Speedup achieved: ~{estimate_runtime(n_combos, 5).split(':')[0]}h → {total_time/3600:.1f}h")

    # Convert to DataFrame
    df = pd.DataFrame(results)

    return df


def save_results(df, timestamp, output_dir='../output_data'):
    """
    Save results to CSV file.

    Parameters:
        df: pd.DataFrame
            Results dataframe
        timestamp: str
            Timestamp string for filename (YYYYMMDD_HHMMSS)
        output_dir: str
            Directory to save results
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{timestamp}_sweep_results.csv')
    df.to_csv(output_file, index=False, encoding='utf-8')
    print(f"\n✓ Results saved to: {output_file}")
    return output_file


def generate_summary_statistics(df, timestamp, output_dir='../output_data'):
    """
    Generate summary statistics and save to text file.

    Parameters:
        df: pd.DataFrame
            Results dataframe
        timestamp: str
            Timestamp string for filename (YYYYMMDD_HHMMSS)
        output_dir: str
            Directory to save summary
    """
    os.makedirs(output_dir, exist_ok=True)
    summary_file = os.path.join(output_dir, f'{timestamp}_sweep_summary.txt')

    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("DELTA-STAR PARAMETER SWEEP SUMMARY\n")
        f.write("="*70 + "\n\n")

        # Overall statistics
        f.write("OVERALL STATISTICS\n")
        f.write("-"*70 + "\n")
        f.write(f"Total combinations tested: {len(df)}\n")
        f.write(f"  - Delta-star found (single): {len(df[df['status'] == 'found'])}\n")
        f.write(f"  - Multiple intersections: {len(df[df['status'] == 'multiple'])}\n")
        f.write(f"  - Not found: {len(df[df['status'] == 'not_found'])}\n")
        if 'error_message' in df.columns:
            f.write(f"  - Errors: {len(df[df['status'] == 'error'])}\n")
        f.write("\n")

        # Delta-star statistics (when found)
        df_found = df[df['primary_delta_star'].notna()]
        if len(df_found) > 0:
            f.write("DELTA-STAR STATISTICS (when found)\n")
            f.write("-"*70 + "\n")
            f.write(f"Average δ*: {df_found['primary_delta_star'].mean():.4f}\n")
            f.write(f"Median δ*: {df_found['primary_delta_star'].median():.4f}\n")
            f.write(f"Std Dev: {df_found['primary_delta_star'].std():.4f}\n")
            f.write(f"Min δ*: {df_found['primary_delta_star'].min():.4f}\n")
            f.write(f"Max δ*: {df_found['primary_delta_star'].max():.4f}\n")
            f.write("\n")

        # Statistics by parameter
        for param in ['k', 'thrshld', 'p1', 'p2', 'beta']:
            f.write(f"\nSTATISTICS BY {param.upper()}\n")
            f.write("-"*70 + "\n")
            grouped = df_found.groupby(param)['primary_delta_star'].agg(['count', 'mean', 'std'])
            f.write(grouped.to_string())
            f.write("\n")

        # Success rate by parameter
        f.write("\nSUCCESS RATE BY PARAMETER\n")
        f.write("-"*70 + "\n")
        for param in ['k', 'thrshld', 'beta']:
            f.write(f"\n{param}:\n")
            grouped = df.groupby(param)['status'].apply(
                lambda x: (x == 'found').sum() / len(x) * 100
            )
            for val, rate in grouped.items():
                f.write(f"  {param}={val}: {rate:.1f}% success\n")

        # Dominant network statistics
        f.write("\n\nDOMINANT NETWORK STATISTICS\n")
        f.write("-"*70 + "\n")

        # Count cases with always-dominant network
        always_dominant_cases = df[df['always_dominant'].notna()]
        if len(always_dominant_cases) > 0:
            f.write(f"\nCases with always-dominant network: {len(always_dominant_cases)}\n")
            dominant_counts = always_dominant_cases['always_dominant'].value_counts()
            for network_type, count in dominant_counts.items():
                pct = count / len(always_dominant_cases) * 100
                f.write(f"  - Always {network_type}: {count} ({pct:.1f}%)\n")

        # Count cases with single delta* (clear transition)
        single_delta_cases = df[df['delta_star_count'] == 1]
        if len(single_delta_cases) > 0:
            f.write(f"\nCases with single δ* (clear transition): {len(single_delta_cases)}\n")
            # Count transition types
            clustered_to_random = len(single_delta_cases[single_delta_cases['dominant_below'] == 'Clustered'])
            random_to_clustered = len(single_delta_cases[single_delta_cases['dominant_below'] == 'Random'])
            f.write(f"  - Clustered → Random transitions: {clustered_to_random}\n")
            f.write(f"  - Random → Clustered transitions: {random_to_clustered}\n")

        # Count cases with multiple delta* (complex behavior)
        multiple_delta_cases = df[df['delta_star_count'] > 1]
        if len(multiple_delta_cases) > 0:
            f.write(f"\nCases with multiple δ* (complex behavior): {len(multiple_delta_cases)}\n")
            delta_count_dist = multiple_delta_cases['delta_star_count'].value_counts().sort_index()
            for count, freq in delta_count_dist.items():
                f.write(f"  - {int(count)} intersections: {freq} cases\n")

        f.write("\n" + "="*70 + "\n")

    print(f"✓ Summary saved to: {summary_file}")
    return summary_file


def save_parameter_info(combinations, timestamp, n_workers, output_dir='../output_data'):
    """
    Save parameter grid information to text file.

    Parameters:
        combinations: list of dicts
            Parameter combinations used in the sweep
        timestamp: str
            Timestamp string for filename (YYYYMMDD_HHMMSS)
        n_workers: int
            Number of parallel workers used
        output_dir: str
            Directory to save parameter info
    """
    os.makedirs(output_dir, exist_ok=True)
    param_file = os.path.join(output_dir, f'{timestamp}_parameters.txt')

    # Extract unique values for each parameter
    k_values = sorted(set(c['k'] for c in combinations))
    thrshld_values = sorted(set(c['thrshld'] for c in combinations))
    p1_values = sorted(set(c['p1'] for c in combinations))
    p2_values = sorted(set(c['p2'] for c in combinations))
    beta_values = sorted(set(c['beta'] for c in combinations))

    # Get trials from first combination (all have same trials value)
    trials = 50  # Default from run_single_combination

    with open(param_file, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("PARAMETER SWEEP CONFIGURATION\n")
        f.write("="*70 + "\n\n")

        # Execution time
        exec_time = datetime.strptime(timestamp, '%Y%m%d_%H%M%S')
        f.write(f"Execution Time: {exec_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Timestamp: {timestamp}\n\n")

        # Parameter ranges
        f.write("PARAMETER RANGES\n")
        f.write("-"*70 + "\n")
        f.write(f"Network Degree (k): {k_values}\n")
        f.write(f"  → {len(k_values)} values\n\n")

        f.write(f"Threshold (thrshld/i): {thrshld_values}\n")
        f.write(f"  → {len(thrshld_values)} values\n\n")

        f.write(f"Base Probability (p1): {p1_values}\n")
        f.write(f"  → {len(p1_values)} values\n\n")

        f.write(f"Reinforced Probability (p2): {p2_values}\n")
        f.write(f"  → {len(p2_values)} values\n\n")

        f.write(f"Time of Influence (beta): {beta_values}\n")
        f.write(f"  → {len(beta_values)} values\n\n")

        # Execution parameters
        f.write("EXECUTION PARAMETERS\n")
        f.write("-"*70 + "\n")
        f.write(f"Trials per combination: {trials}\n")
        f.write(f"Parallel workers: {n_workers}\n")
        f.write(f"Network size (n): k × 250 nodes\n\n")

        # Combination statistics
        f.write("COMBINATION STATISTICS\n")
        f.write("-"*70 + "\n")
        f.write(f"Total parameter combinations: {len(combinations)}\n")
        f.write(f"  (Only combinations with p1 ≤ p2 are included)\n\n")

        f.write(f"Theoretical maximum: {len(k_values) * len(thrshld_values) * len(p1_values) * len(p2_values) * len(beta_values)}\n")
        f.write(f"Valid combinations: {len(combinations)}\n")
        f.write(f"Filtered out: {len(k_values) * len(thrshld_values) * len(p1_values) * len(p2_values) * len(beta_values) - len(combinations)} (p1 > p2)\n\n")

        # Additional notes
        f.write("NOTES\n")
        f.write("-"*70 + "\n")
        f.write("• seeds = threshold (i) for all combinations\n")
        f.write("• Network type: Watts-Strogatz (WS)\n")
        f.write("• Seeding strategy: one_nbrs for clustered, random for random\n")
        f.write("• Delta values tested: np.linspace(0, 1, 1001)\n")
        f.write("• Minimum delta threshold: 0.0001\n")
        f.write("• Full spread threshold: 60% of population\n\n")

        f.write("="*70 + "\n")

    print(f"✓ Parameter info saved to: {param_file}")
    return param_file


def create_visualizations(df, timestamp, output_dir='../output_plots/sweep_plots'):
    """
    Create visualization plots for the sweep results.

    Parameters:
        df: pd.DataFrame
            Results dataframe
        timestamp: str
            Timestamp string for filename (YYYYMMDD_HHMMSS)
        output_dir: str
            Directory to save plots
    """
    os.makedirs(output_dir, exist_ok=True)

    # Filter to only rows where delta_star was found
    df_found = df[df['primary_delta_star'].notna()].copy()

    if len(df_found) == 0:
        print("⚠ No delta-star values found - skipping visualizations")
        return

    # Set style
    if HAS_SEABORN:
        sns.set_style("whitegrid")

    # 1. Heatmap: k vs thrshld
    print("Creating heatmap (k vs thrshld)...")
    fig, ax = plt.subplots(figsize=(10, 6))
    pivot = df_found.pivot_table(
        values='primary_delta_star',
        index='k',
        columns='thrshld',
        aggfunc='mean'
    )
    if HAS_SEABORN:
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='viridis', ax=ax,
                    cbar_kws={'label': 'Average δ*'})
    else:
        # Fallback to matplotlib's imshow
        im = ax.imshow(pivot, cmap='viridis', aspect='auto')
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_yticks(range(len(pivot.index)))
        ax.set_xticklabels(pivot.columns)
        ax.set_yticklabels(pivot.index)
        plt.colorbar(im, ax=ax, label='Average δ*')
        # Add value annotations
        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                text = ax.text(j, i, f'{pivot.iloc[i, j]:.3f}',
                             ha="center", va="center", color="w")
    ax.set_title('Average Delta-Star by Network Degree and Threshold', fontweight='bold')
    ax.set_xlabel('Threshold (i)', fontweight='bold')
    ax.set_ylabel('Network Degree (k)', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{timestamp}_heatmap_k_vs_thrshld.png'), dpi=300)
    plt.close()

    # 2. Scatter plot: p1 vs p2
    print("Creating scatter plot (p1 vs p2)...")
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(
        df_found['p1'],
        df_found['p2'],
        c=df_found['primary_delta_star'],
        s=df_found['beta'] * 50,  # Size by beta
        cmap='coolwarm',
        alpha=0.6,
        edgecolors='black',
        linewidth=0.5
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('δ*', fontweight='bold')
    ax.set_xlabel('Base Adoption Probability (p₁)', fontweight='bold')
    ax.set_ylabel('Reinforced Adoption Probability (p₂)', fontweight='bold')
    ax.set_title('Delta-Star by Adoption Probabilities\n(point size = β)', fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Add legend for beta sizes
    for beta_val in sorted(df_found['beta'].unique()):
        ax.scatter([], [], s=beta_val*50, c='gray', alpha=0.6,
                  edgecolors='black', linewidth=0.5,
                  label=f'β={beta_val}')
    ax.legend(scatterpoints=1, frameon=True, labelspacing=1, title='Time of Influence')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{timestamp}_scatter_p1_vs_p2.png'), dpi=300)
    plt.close()

    # 3. Box plot: beta vs delta_star
    print("Creating box plot (beta distribution)...")
    fig, ax = plt.subplots(figsize=(8, 6))
    df_found.boxplot(column='primary_delta_star', by='beta', ax=ax)
    ax.set_xlabel('Time of Influence (β)', fontweight='bold')
    ax.set_ylabel('Delta-Star (δ*)', fontweight='bold')
    ax.set_title('Delta-Star Distribution by Time of Influence', fontweight='bold')
    plt.suptitle('')  # Remove default title
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{timestamp}_boxplot_beta.png'), dpi=300)
    plt.close()

    # 4. Pie chart: status distribution
    print("Creating pie chart (status distribution)...")
    fig, ax = plt.subplots(figsize=(8, 8))
    status_counts = df['status'].value_counts()
    colors = {'found': '#2ecc71', 'multiple': '#f39c12', 'not_found': '#e74c3c', 'error': '#95a5a6'}
    pie_colors = [colors.get(status, '#95a5a6') for status in status_counts.index]

    wedges, texts, autotexts = ax.pie(
        status_counts.values,
        labels=status_counts.index,
        autopct='%1.1f%%',
        startangle=90,
        colors=pie_colors,
        textprops={'fontweight': 'bold'}
    )
    ax.set_title('Distribution of Sweep Results', fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{timestamp}_pie_chart_status.png'), dpi=300)
    plt.close()

    print(f"\n✓ Visualizations saved to: {output_dir}/")


def main():
    """
    Main function to run the complete parameter sweep.
    """
    # Generate timestamp for this run
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print("="*70)
    print("DELTA-STAR PARAMETER SWEEP")
    print("Systematic exploration of parameter space")
    print("="*70)
    print(f"\nRun Timestamp: {timestamp}")
    print(f"  ({datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")

    # Generate parameter combinations
    print("\nGenerating parameter combinations...")
    combinations = generate_parameter_combinations()
    print(f"✓ Generated {len(combinations)} valid combinations (p₁ ≤ p₂)")

    # Get number of workers
    n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))

    # Save parameter information
    print("\nSaving parameter configuration...")
    save_parameter_info(combinations, timestamp, n_workers)

    # Run the sweep
    results_df = run_parameter_sweep(combinations, n_workers=n_workers)

    # Save results
    save_results(results_df, timestamp)

    # Generate summary statistics
    generate_summary_statistics(results_df, timestamp)

    # Create visualizations
    create_visualizations(results_df, timestamp)

    print("\n" + "="*70)
    print("SWEEP COMPLETE!")
    print("="*70)
    print(f"\nFiles created (timestamp: {timestamp}):")
    print(f"  - ../output_data/{timestamp}_parameters.txt")
    print(f"  - ../output_data/{timestamp}_sweep_results.csv")
    print(f"  - ../output_data/{timestamp}_sweep_summary.txt")
    print(f"  - ../output_plots/sweep_plots/{timestamp}_*.png")
    print()


if __name__ == "__main__":
    main()
