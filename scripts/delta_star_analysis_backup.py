#MEMO
# 초기 시작점 확인 

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


# δ* (delta-star) 임계값 탐색 실험: 시간적 할인 요인을 통한 네트워크 성능 비교
"""
This script implements an experiment to find the critical temporal discount factor δ*
where clustered and random networks perform equally in complex contagion simulations.

Key Concepts:
- δ (delta): Temporal discount factor ranging from 0 to 1
  - When δ is close to 0: Only immediate spread matters
  - When δ is close to 1: All time steps are valued equally

- V: Discounted cumulative activation size
  - Formula: V = c₁×δ + c₂×δ² + c₃×δ³ + ... + cₜ×δᵗ
  - where cₜ is the number of newly activated nodes at time step t

- δ*: Critical threshold where V_clustered = V_random
  - Below δ*: Clustered networks outperform (reinforcement matters more)
  - Above δ*: Random networks outperform (reach matters more)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from helper_functions import (
    sirsic_wmem, run_simulation,
    seed_strat_one_nbrs, seed_strat_random, seed_strat_adj
)

###################
# Core Functions for δ* Analysis

def cumulative_to_new_infections(cumulative_series):
    """
    Convert cumulative adoption counts to per-timestep new infections.

    Parameters:
        cumulative_series: list of integers
            Cumulative number of adopters at each time step
            Example: [5, 8, 12, 15, 15, 15] means 5 at t=0, 8 at t=1, etc.

    Returns:
        new_infections: list of integers
            Number of newly infected nodes at each time step
            Example: [5, 3, 4, 3, 0, 0] for the above input
    """
    # Initialize with the first value (number of initial seeds)
    # This represents how many nodes adopted at timestep 0
    new_infections = [cumulative_series[0]]

    # Loop through remaining timesteps to calculate new adoptions
    for t in range(1, len(cumulative_series)):
        # Calculate new infections = cumulative at t - cumulative at t-1
        # This gives us how many NEW nodes adopted at this specific timestep
        new_at_t = cumulative_series[t] - cumulative_series[t-1]
        new_infections.append(new_at_t)

    return new_infections


def calculate_V(new_infections, delta):
    """
    Calculate the discounted cumulative activation size V.

    This metric applies temporal discounting to value early adoptions more heavily
    than later ones (when delta < 1).

    Formula: V = Σ(cₜ × δᵗ) for t = 1, 2, 3, ...

    Parameters:
        new_infections: list of integers
            Number of newly activated nodes at each time step (cₜ values)
        delta: float between 0 and 1
            Temporal discount factor
            - delta = 0: Only first timestep matters
            - delta = 1: All timesteps weighted equally (equivalent to total spread)

    Returns:
        V: float
            Discounted cumulative activation size

    Example:
        If new_infections = [5, 3, 4, 2] and delta = 0.9:
        V = 5×0.9¹ + 3×0.9² + 4×0.9³ + 2×0.9⁴
          = 5×0.9 + 3×0.81 + 4×0.729 + 2×0.6561
          = 4.5 + 2.43 + 2.916 + 1.3122
          = 11.1582
    """
    # Initialize the discounted value V to zero
    V = 0.0

    # Loop through each timestep and apply exponential discounting
    for t in range(len(new_infections)):
        # Apply discount factor: delta^(t+1)
        # - t+1 is used because we start counting from t=0 but discount from time 1
        # - Lower delta values make later adoptions worth less
        # - When delta=0.9, adoptions at t=0 are worth 0.9, at t=1 worth 0.81, etc.
        V += new_infections[t] * (delta ** (t + 1))

    return V


def run_delta_simulations(k=8, n=None, thrshld=2, p1=0.3, p2=0.6,
                          beta=1, trials=20, seeds=5):
    """
    Run simulations for both clustered and random networks with full time series.

    Parameters:
        k: int
            Network degree (default: 8)
        n: int or None
            Number of nodes. If None, uses k*250 (default: None)
        thrshld: int
            Social reinforcement threshold (default: 2)
        p1: float
            Base adoption probability (default: 0.3)
        p2: float
            Reinforced adoption probability (default: 0.6)
        beta: int
            Time of influence (default: 1)
        trials: int
            Number of simulation trials per network type (default: 20)
        seeds: int
            Number of initial seeds (default: 5)

    Returns:
        results_clustered: DataFrame
            Simulation results for clustered network (perc_rewire=0)
        results_random: DataFrame
            Simulation results for random network (perc_rewire=1)
    """
    # Set default network size if not specified (k*250 is a common convention)
    if n is None:
        n = k * 250

    # Display simulation parameters for user verification
    print(f"Setting up simulations with parameters:")
    print(f"  Network: k={k}, n={n}")
    print(f"  Diffusion: i={thrshld}, p1={p1}, p2={p2}, beta={beta}")
    print(f"  Trials: {trials} per network type\n")

    # Create base network using Watts-Strogatz model
    # - n: number of nodes
    # - k: each node connects to k nearest neighbors
    # - 0: rewiring probability (0 = no rewiring, creates regular ring lattice)
    # - seed=42: fixed seed for reproducibility
    print("Creating base network (WS ring lattice)...")
    G = nx.connected_watts_strogatz_graph(n, k, 0, seed=42)

    # Build the adoption probability trajectory
    # - prob_list[i] = probability of adoption when exposed to i influential neighbors
    # - [0]: Can't adopt with 0 exposures
    # - [p1] * (thrshld-1): Low probability (p1) below threshold
    # - [p2] * n: High probability (p2) at or above threshold
    # Example with thrshld=2: [0, p1, p2, p2, p2, ...]
    prob_list = [0] + [p1] * (thrshld - 1) + [p2] * n

    # Run simulations for CLUSTERED network (perc_rewire = 0)
    # Clustered networks have high local clustering (friends-of-friends know each other)
    print("Running simulations on CLUSTERED network (perc_rewire=0)...")
    _, cf_clustered = run_simulation(
        G=G,                              # Base network structure
        G_name='WS',                      # Watts-Strogatz network type
        model=sirsic_wmem,                # Diffusion model with memory
        trials=trials,                    # Number of independent simulation runs
        r_start=0,                        # Starting trial number
        beta=beta,                        # Time of influence (memory parameter)
        prob_list=prob_list,              # Adoption probabilities by exposure count
        perc=0,                           # 0 = no rewiring = CLUSTERED network
        rand_path=None,                   # No pre-saved random network file
        seeds=seeds,                      # Number of initial adopters
        seed_strat0=seed_strat_one_nbrs,  # First seeding: select neighbor pairs
        seed_strat=seed_strat_random,     # Subsequent seeding: random selection
        sig_thresh_sd=0,                  # 0 = no heterogeneity in thresholds
        full_series=True,                 # CRITICAL: Save full time series for V calculation
        rand_network='ES'                 # 'ES' = Edge Swap randomization method
    )

    # Run simulations for RANDOM network (perc_rewire = 1)
    # Random networks have low clustering but shorter path lengths (more reach)
    print("Running simulations on RANDOM network (perc_rewire=1)...")
    _, cf_random = run_simulation(
        G=G,                              # Same base network structure
        G_name='WS',
        model=sirsic_wmem,
        trials=trials,
        r_start=0,
        beta=beta,
        prob_list=prob_list,
        perc=1,                           # 1 = full rewiring = RANDOM network
        rand_path=None,
        seeds=seeds,
        seed_strat0=seed_strat_one_nbrs,
        seed_strat=seed_strat_random,
        sig_thresh_sd=0,
        full_series=True,                 # CRITICAL: Save full time series for V calculation
        rand_network='ES'
    )

    print(f"\nSimulations complete!")
    print(f"  Clustered trials: {len(cf_clustered)}")
    print(f"  Random trials: {len(cf_random)}")

    return cf_clustered, cf_random


def compute_V_across_deltas(results_df, delta_values):
    """
    Compute V (discounted cumulative activation) for each trial across different δ values.

    Parameters:
        results_df: DataFrame
            Simulation results containing 'full_timeseries' column
        delta_values: list of floats
            Delta values to test (e.g., [0.5, 0.7, 0.9, 0.95, 0.99])

    Returns:
        V_matrix: numpy array of shape (n_trials, n_deltas)
            V values for each trial (rows) and delta value (columns)
    """
    # Get dimensions for the output matrix
    n_trials = len(results_df)
    n_deltas = len(delta_values)

    # Initialize matrix to store V values
    # Rows = individual trials, Columns = different delta values
    V_matrix = np.zeros((n_trials, n_deltas))

    # Loop through each simulation trial
    for trial_idx in range(n_trials):
        # Extract the cumulative adoption time series for this trial
        # This is a list like [5, 8, 12, 15, 15, ...] showing total adopters over time
        cumulative_series = results_df.iloc[trial_idx]['full_timeseries']

        # Convert cumulative counts to new infections per timestep
        # Transforms [5, 8, 12, 15] → [5, 3, 4, 3]
        new_infections = cumulative_to_new_infections(cumulative_series)

        # Calculate V for each delta value using the same trial data
        # This shows how the same diffusion trajectory is valued under different time preferences
        for delta_idx, delta in enumerate(delta_values):
            V = calculate_V(new_infections, delta)
            V_matrix[trial_idx, delta_idx] = V

    return V_matrix


def find_delta_star(delta_values, V_clustered_mean, V_random_mean, min_delta=0.05):
    """
    Find ALL δ* values where clustered and random networks have equal performance.

    Uses linear interpolation to find intersection points, filtering out trivial
    intersections near delta=0.

    Parameters:
        delta_values: array-like
            Delta values tested
        V_clustered_mean: array-like
            Mean V values for clustered network at each delta
        V_random_mean: array-like
            Mean V values for random network at each delta
        min_delta: float
            Minimum delta threshold to exclude trivial intersections (default: 0.05)

    Returns:
        tuple: (primary_delta_star, all_intersections)
            - primary_delta_star: float or None
                The first meaningful delta* value (None if no intersection found)
            - all_intersections: list of dicts
                List of all meaningful intersections, each containing:
                    'delta_star': float - the intersection value
                    'interval': tuple - (lower_delta, upper_delta) bracketing the intersection
                    'transition': str - direction of transition ('clustered_to_random' or 'random_to_clustered')
    """
    # Calculate the difference between random and clustered networks at each delta
    # Positive values = random network performing better
    # Negative values = clustered network performing better
    diff = np.array(V_random_mean) - np.array(V_clustered_mean)

    # Find where the difference changes sign (crosses zero)
    # np.diff(np.sign(diff)) will be non-zero where sign changes
    # This identifies the interval where the curves intersect
    sign_changes = np.where(np.diff(np.sign(diff)))[0]

    # List to store all meaningful intersections
    all_intersections = []

    # Process each sign change to find exact intersection points
    for idx in sign_changes:
        # Linear interpolation to find the exact delta* value
        # We know the intersection is between delta_values[idx] and delta_values[idx+1]
        x1, x2 = delta_values[idx], delta_values[idx + 1]  # Delta values bracketing the intersection
        y1, y2 = diff[idx], diff[idx + 1]                   # Difference values at those deltas

        # Use linear interpolation formula to find where the line crosses y=0
        # Formula: x* = x1 - y1 * (x2 - x1) / (y2 - y1)
        # This solves: y = y1 + (y2-y1)/(x2-x1) * (x - x1) = 0 for x
        delta_star_candidate = x1 - y1 * (x2 - x1) / (y2 - y1)

        # Filter out trivial intersections near delta=0
        if delta_star_candidate >= min_delta:
            # Determine transition direction
            # If diff goes from negative to positive: clustered was better, now random is better
            # If diff goes from positive to negative: random was better, now clustered is better
            if y1 < 0 and y2 > 0:
                transition = 'clustered_to_random'
            elif y1 > 0 and y2 < 0:
                transition = 'random_to_clustered'
            else:
                # Handle edge case where one of the values is exactly zero
                transition = 'clustered_to_random' if y2 > y1 else 'random_to_clustered'

            # Store this intersection
            all_intersections.append({
                'delta_star': delta_star_candidate,
                'interval': (x1, x2),
                'transition': transition
            })

    # Print detailed analysis of intersections
    print("\n" + "="*70)
    print("INTERSECTION ANALYSIS")
    print("="*70)

    if len(all_intersections) == 0:
        print(f"Found 0 meaningful intersections (threshold: δ ≥ {min_delta})")
        if len(sign_changes) > 0:
            print(f"\nNote: {len(sign_changes)} intersection(s) were found but all were")
            print(f"      below the minimum threshold of {min_delta} (trivial intersections).")
        else:
            print("\nNo intersections found in the given delta range.")
        print("="*70 + "\n")
        return None, []

    print(f"Found {len(all_intersections)} meaningful intersection(s):\n")

    # Print details for each intersection
    for i, intersection in enumerate(all_intersections, 1):
        print(f"Intersection {i}:")
        print(f"  δ* = {intersection['delta_star']:.4f}")
        print(f"  Interval: [{intersection['interval'][0]:.2f}, {intersection['interval'][1]:.2f}]")
        # Format transition direction for readability
        if intersection['transition'] == 'clustered_to_random':
            print(f"  Transition: Clustered → Random (clustered advantage ends)")
        else:
            print(f"  Transition: Random → Clustered (random advantage ends)")
        print()

    # Check for multiple intersections and warn about non-monotonic behavior
    if len(all_intersections) > 1:
        print("⚠ WARNING: MULTIPLE intersections detected - complex behavior!")
        print("  The relative performance of networks is NON-MONOTONIC with δ.")
        print("  This indicates that the optimal network structure changes multiple")
        print("  times as the temporal preference shifts.")
        print(f"  Recommending PRIMARY δ* = {all_intersections[0]['delta_star']:.4f} (first meaningful transition)")

    print("="*70 + "\n")

    # Return the first meaningful intersection as primary (for backward compatibility)
    # along with the complete list of all intersections
    primary_delta_star = all_intersections[0]['delta_star']
    return primary_delta_star, all_intersections


def visualize_delta_analysis(delta_values, V_clustered_mean, V_clustered_se,
                             V_random_mean, V_random_se, delta_star=None, all_intersections=None):
    """
    Create visualization showing how V changes with δ for both network types.

    Parameters:
        delta_values: array-like
            Delta values tested
        V_clustered_mean: array-like
            Mean V for clustered networks
        V_clustered_se: array-like
            Standard error for clustered networks
        V_random_mean: array-like
            Mean V for random networks
        V_random_se: array-like
            Standard error for random networks
        delta_star: float or None
            Primary critical delta value to mark on plot (for backward compatibility)
        all_intersections: list of dicts or None
            List of all intersection dictionaries with keys 'delta_star', 'interval', 'transition'
    """
    # Create figure and axis objects
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot clustered network performance curve
    # 'o-' = circles connected by lines
    ax.plot(delta_values, V_clustered_mean, 'o-', color='blue',
            linewidth=2, markersize=8, label='Clustered Network')
    # Add shaded region for standard error (uncertainty band)
    ax.fill_between(delta_values,
                    V_clustered_mean - V_clustered_se,  # Lower bound
                    V_clustered_mean + V_clustered_se,  # Upper bound
                    color='blue', alpha=0.2)            # Semi-transparent blue

    # Plot random network performance curve
    # 's-' = squares connected by lines
    ax.plot(delta_values, V_random_mean, 's-', color='red',
            linewidth=2, markersize=8, label='Random Network')
    # Add shaded region for standard error
    ax.fill_between(delta_values,
                    V_random_mean - V_random_se,
                    V_random_mean + V_random_se,
                    color='red', alpha=0.2)

    # Mark the critical threshold(s) δ* if found
    if delta_star is not None and all_intersections is not None:
        # Draw all intersections with appropriate styling
        for i, intersection in enumerate(all_intersections):
            delta_val = intersection['delta_star']

            if i == 0:
                # PRIMARY intersection: thick black dashed line with large label
                ax.axvline(delta_val, color='black', linestyle='--', linewidth=2.5,
                          label=f'δ* = {delta_val:.3f}', zorder=10)
                # Add prominent text label at the top
                ax.text(delta_val, ax.get_ylim()[1] * 0.95, f'  δ* = {delta_val:.3f}',
                       rotation=90, verticalalignment='top', fontsize=12,
                       fontweight='bold', color='black', zorder=11)
            else:
                # SECONDARY intersections: thinner gray dotted lines with smaller labels
                ax.axvline(delta_val, color='gray', linestyle=':', linewidth=1.5,
                          alpha=0.7, label=f'δ*₍{i+1}₎ = {delta_val:.3f}', zorder=9)
                # Add smaller text label
                ax.text(delta_val, ax.get_ylim()[1] * 0.88, f'  δ*₍{i+1}₎ = {delta_val:.3f}',
                       rotation=90, verticalalignment='top', fontsize=10,
                       color='gray', alpha=0.8, zorder=9)

    elif delta_star is not None:
        # Backward compatibility: single delta_star provided without all_intersections
        # Draw vertical dashed line at δ*
        ax.axvline(delta_star, color='black', linestyle='--', linewidth=2,
                  label=f'δ* = {delta_star:.3f}')
        # Add text label at the top of the line
        ax.text(delta_star, ax.get_ylim()[1] * 0.95, f'  δ* = {delta_star:.3f}',
               rotation=90, verticalalignment='top', fontsize=12, fontweight='bold')
    else:
        # δ* was not found - add a prominent warning annotation
        # Position: center-top of the plot for maximum visibility
        # Get the middle x position (in data coordinates)
        x_mid = (delta_values[0] + delta_values[-1]) / 2
        # Get the upper y position (95% of the way up)
        y_top = ax.get_ylim()[1] * 0.95

        # Create the warning text with line breaks for readability
        warning_text = "⚠ δ* NOT FOUND ⚠\nNo intersection in tested range\n(Consider expanding delta values)"

        # Add the text annotation with a highly visible yellow box and red text
        ax.text(x_mid, y_top, warning_text,
               horizontalalignment='center',      # Center the text horizontally
               verticalalignment='top',           # Anchor at top
               fontsize=13,                       # Large, readable font
               fontweight='bold',                 # Bold for emphasis
               color='darkred',                   # Dark red text for warning
               bbox=dict(
                   boxstyle='round,pad=0.8',      # Rounded corners with padding
                   facecolor='yellow',            # Yellow background (warning color)
                   edgecolor='red',               # Red border for extra emphasis
                   linewidth=3,                   # Thick border
                   alpha=0.9                      # Slightly transparent
               ),
               zorder=100)                        # Draw on top of everything

    # Set axis labels with formatting
    ax.set_xlabel('Temporal Discount Factor (δ)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Discounted Cumulative Activation (V)', fontsize=14, fontweight='bold')
    ax.set_title('Network Performance vs. Temporal Discounting\n' +
                'Finding δ*: The Critical Threshold', fontsize=16, fontweight='bold')

    # Add legend to identify the curves
    ax.legend(fontsize=12, loc='best')

    # Add grid for easier reading
    ax.grid(True, alpha=0.3)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Create output directory if it doesn't exist
    os.makedirs('../output_plots', exist_ok=True)

    # Save the figure as high-resolution PNG
    plt.savefig('../output_plots/delta_star_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nVisualization saved to: output_plots/delta_star_analysis.png")

    # Display the plot
    plt.show()

    return fig, ax


###################
# Main Execution

def main(k=4, n=None, thrshld=4, p1=0.75, p2=0.75, beta=1, trials=50, seeds=None,
         delta_values=None, verbose=True):
    """
    Main function to run the complete δ* analysis.

    Parameters:
        k: int
            Network degree (default: 4)
        n: int or None
            Number of nodes. If None, uses k*250 (default: None)
        thrshld: int
            Social reinforcement threshold (default: 4)
        p1: float
            Base adoption probability (default: 0.75)
        p2: float
            Reinforced adoption probability (default: 0.75)
        beta: int
            Time of influence (default: 1)
        trials: int
            Number of simulation trials per network type (default: 50)
        seeds: int or None
            Number of initial seeds. If None, uses thrshld (default: None)
        delta_values: list or None
            Delta values to test. If None, uses np.linspace(0, 1, 51) (default: None)
        verbose: bool
            If True, print detailed output. If False, suppress most output (default: True)

    Returns:
        tuple: (delta_star, results_summary, all_intersections)
            - delta_star: float or None (primary intersection)
            - results_summary: DataFrame with detailed results
            - all_intersections: list of all intersection dicts

    This function orchestrates the entire workflow:
    1. Define simulation parameters
    2. Run simulations on both network types
    3. Calculate V for different δ values
    4. Find the critical threshold δ*
    5. Visualize and save results
    """
    # Set defaults for None parameters
    if n is None:
        n = k * 250
    if seeds is None:
        seeds = thrshld
    if delta_values is None:
        delta_values = np.linspace(0, 1, 51).tolist()

    # Suppress output if not verbose
    import sys
    import io
    old_stdout = None
    if not verbose:
        # Redirect stdout to suppress print statements
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()

    try:
        # Print header for the analysis
        print("="*70)
        print("δ* (DELTA-STAR) ANALYSIS")
        print("Finding the critical temporal discount threshold")
        print("="*70 + "\n")

        # ========== STEP 1: Define parameters ==========
        print("STEP 1: Setting parameters")
        print("-" * 50)
        print(f"Delta values to test: {delta_values}\n")

        # ========== STEP 2: Run simulations ==========
        print("STEP 2: Running simulations")
        print("-" * 50)

        # Run simulations on both clustered (perc=0) and random (perc=1) networks
        # This will give us time series data for each trial
        cf_clustered, cf_random = run_delta_simulations(
            k=k, n=n, thrshld=thrshld, p1=p1, p2=p2,
            beta=beta, trials=trials, seeds=seeds
        )

        # ========== STEP 3: Calculate V for different delta values ==========
        print("\nSTEP 3: Calculating V for different δ values")
        print("-" * 50)

        # For clustered network: compute V at each delta value for all trials
        # This creates a matrix: rows=trials, columns=delta values
        print("Computing V for clustered network...")
        V_clustered = compute_V_across_deltas(cf_clustered, delta_values)

        # Same for random network
        print("Computing V for random network...")
        V_random = compute_V_across_deltas(cf_random, delta_values)

        # Calculate summary statistics across trials
        # Mean: average performance at each delta
        V_clustered_mean = np.mean(V_clustered, axis=0)
        # Standard error: uncertainty in the mean estimate
        V_clustered_se = np.std(V_clustered, axis=0) / np.sqrt(trials)

        V_random_mean = np.mean(V_random, axis=0)
        V_random_se = np.std(V_random, axis=0) / np.sqrt(trials)

        # Print results table showing performance comparison at each delta
        print("\nResults Summary:")
        print("=" * 70)
        print(f"{'δ':>8} | {'V_clustered':>15} | {'V_random':>15} | {'Difference':>15}")
        print("-" * 70)
        for i, delta in enumerate(delta_values):
            # Calculate difference: positive means random is better, negative means clustered is better
            diff = V_random_mean[i] - V_clustered_mean[i]
            winner = "Random" if diff > 0 else "Clustered"
            print(f"{delta:>8.2f} | {V_clustered_mean[i]:>15.2f} | "
            f"{V_random_mean[i]:>15.2f} | {diff:>+15.2f} ({winner})")
        print("=" * 70)

        # ========== STEP 4: Find delta* ==========
        print("\nSTEP 4: Finding δ* (critical threshold)")
        print("-" * 50)

        # Find ALL delta values where the two curves intersect
        # This is the critical threshold (or thresholds) that separate the two regimes
        delta_star, all_intersections = find_delta_star(delta_values, V_clustered_mean, V_random_mean)

        # Interpret the result
        if delta_star is not None:
            print(f"\n*** PRIMARY δ* found: {delta_star:.4f} ***")

            # Check for suspiciously low delta_star values
        if delta_star < 0.1:
            print(f"\n⚠ WARNING: δ* = {delta_star:.4f} is suspiciously LOW!")
            print(f"  This indicates an extreme early-adoption regime.")
            print(f"  Potential issues to check:")
            print(f"    - Is the minimum threshold (0.05) too permissive?")
            print(f"    - Are the simulation parameters creating trivial dynamics?")
            print(f"    - Consider validating with different parameter settings.")
            print()

            # Check for multiple intersections and provide extended interpretation
        if len(all_intersections) > 1:
            print(f"\n⚠ COMPLEX BEHAVIOR DETECTED:")
            print(f"  Multiple intersections indicate NON-MONOTONIC relative performance.")
            print(f"  The optimal network structure CHANGES {len(all_intersections)} times as δ varies.")
            print(f"\nDetailed regime breakdown:")

            # Add beginning regime
            first_transition = all_intersections[0]['transition']
            if first_transition == 'clustered_to_random':
                print(f"  1. δ < {all_intersections[0]['delta_star']:.4f}: Clustered networks dominate")
                print(f"     (High clustering advantage in early-adoption regime)")
            else:
                print(f"  1. δ < {all_intersections[0]['delta_star']:.4f}: Random networks dominate")
                print(f"     (High reach advantage in early-adoption regime)")

            # Add intermediate regimes
            for i in range(len(all_intersections)):
                intersection = all_intersections[i]
                if i < len(all_intersections) - 1:
                    next_intersection = all_intersections[i + 1]
                    # Determine which network wins in this interval
                    if intersection['transition'] == 'clustered_to_random':
                        winner = "Random"
                        reason = "Reach advantage dominates"
                    else:
                        winner = "Clustered"
                        reason = "Reinforcement advantage dominates"

                    print(f"  {i+2}. {intersection['delta_star']:.4f} < δ < {next_intersection['delta_star']:.4f}: {winner} networks dominate")
                    print(f"     ({reason})")
                else:
                    # Final regime (after last intersection)
                    if intersection['transition'] == 'clustered_to_random':
                        print(f"  {i+2}. δ > {intersection['delta_star']:.4f}: Random networks dominate")
                        print(f"     (Reach advantage in late-adoption regime)")
                    else:
                        print(f"  {i+2}. δ > {intersection['delta_star']:.4f}: Clustered networks dominate")
                        print(f"     (Reinforcement advantage in late-adoption regime)")

            print(f"\n  💡 Interpretation: This non-monotonic pattern suggests complex")
            print(f"     interaction between timing, reinforcement, and reach effects.")
            print(f"     The 'best' network depends critically on temporal preferences.")

        else:
            # Single intersection - standard interpretation
            print(f"\nInterpretation:")
            print(f"  • For δ < {delta_star:.4f}: Clustered networks outperform")
            print(f"    → Social reinforcement is more valuable")
            print(f"    → Early, reinforced adoptions matter most")
            print(f"  • For δ > {delta_star:.4f}: Random networks outperform")
            print(f"    → Network reach is more valuable")
            print(f"    → Spreading widely over time matters most")
        else:
            # δ* was not found in the tested range - provide detailed analysis
            print("\n" + "!"*70)
            print("*** δ* NOT FOUND in the tested range ***")
            print("!"*70)

            # Analyze which network type dominates across the tested range
            differences = V_random_mean - V_clustered_mean
            random_wins = np.sum(differences > 0)
            clustered_wins = np.sum(differences < 0)
            total_tests = len(differences)

            # Calculate average magnitude of differences to assess closeness
            avg_abs_diff = np.mean(np.abs(differences))
            # Normalize by mean V values to get relative difference
            avg_V = (np.mean(V_clustered_mean) + np.mean(V_random_mean)) / 2
            relative_diff = avg_abs_diff / avg_V if avg_V > 0 else 0

            print(f"\nNetwork Performance Summary:")
            print(f"  • Random network wins: {random_wins}/{total_tests} delta values")
            print(f"  • Clustered network wins: {clustered_wins}/{total_tests} delta values")
            print(f"  • Average relative difference: {relative_diff:.2%}")

            # Determine dominant pattern and provide interpretation
            if random_wins == total_tests:
            # Random network always wins
            print(f"\n⚠ Analysis: Random network ALWAYS outperforms in tested range")
            print(f"\nInterpretation:")
            print(f"  • δ* is likely BELOW the tested range (δ < {delta_values[0]})")
            print(f"  • OR δ* may not exist for these parameters")
            print(f"  • Network 'reach' dominates over 'reinforcement' across all tested δ values")
            print(f"  • Even when valuing early adoptions heavily, random networks spread better")
            print(f"\n💡 Recommendation:")
            print(f"  → Test LOWER delta values: try [0.01, 0.02, 0.05, 0.1, ...]")
            print(f"  → If δ* exists, it's in the extreme early-adoption regime")

            elif clustered_wins == total_tests:
            # Clustered network always wins
            print(f"\n⚠ Analysis: Clustered network ALWAYS outperforms in tested range")
            print(f"\nInterpretation:")
            print(f"  • δ* is likely ABOVE the tested range (δ > {delta_values[-1]})")
            print(f"  • Social 'reinforcement' dominates over 'reach' across all tested δ values")
            print(f"  • Even when valuing late adoptions equally, clustered networks perform better")
            print(f"\n💡 Recommendation:")
            print(f"  → Test HIGHER delta values: try [..., 0.995, 0.999, 0.9999]")
            print(f"  → If δ* exists, it's in the extreme late-adoption regime")

            elif relative_diff < 0.05:
            # Networks perform very similarly
            print(f"\n⚠ Analysis: Networks perform VERY SIMILARLY across tested range")
            print(f"\nInterpretation:")
            print(f"  • Average difference is only {relative_diff:.2%} - negligible")
            print(f"  • δ* may not be a meaningful concept for these parameters")
            print(f"  • Network structure has minimal impact on diffusion performance")
            print(f"\n💡 Recommendation:")
            print(f"  → Consider testing different parameter combinations:")
            print(f"     - Different thresholds (current: i={thrshld})")
            print(f"     - Different p1/p2 ratios (current: p1={p1}, p2={p2})")
            print(f"     - Different beta values (current: β={beta})")

            else:
            # Mixed results - intersection might be near boundaries
            print(f"\n⚠ Analysis: Mixed results - performance varies across range")
            print(f"\nInterpretation:")
            print(f"  • Both networks win in different regimes, but no clear intersection found")
            print(f"  • δ* might be very close to one of the tested values")
            print(f"  • OR the transition is non-monotonic (curves cross multiple times)")
            print(f"\n💡 Recommendation:")
            print(f"  → Test with FINER GRANULARITY around transition regions")
            print(f"  → Current spacing: {delta_values[1]-delta_values[0]:.3f} to {delta_values[-1]-delta_values[-2]:.3f}")
            print(f"  → Try denser grid: np.linspace({delta_values[0]}, {delta_values[-1]}, 50)")

            print("\n" + "!"*70)

        # ========== STEP 5: Visualize results ==========
        print("\n\nSTEP 5: Creating visualization")
        print("-" * 50)

        # Create and save the plot showing both curves and all δ* intersections
        visualize_delta_analysis(
            delta_values, V_clustered_mean, V_clustered_se,
            V_random_mean, V_random_se, delta_star, all_intersections
            )

        # Save detailed numerical results to CSV for further analysis
        results_summary = pd.DataFrame({
            'delta': delta_values,
            'V_clustered_mean': V_clustered_mean,
            'V_clustered_se': V_clustered_se,
            'V_random_mean': V_random_mean,
            'V_random_se': V_random_se,
            'difference': V_random_mean - V_clustered_mean  # Positive = random better
            })

        # Create output directory if it doesn't exist
        os.makedirs('../output_data', exist_ok=True)

        # Save to CSV
        output_file = '../output_data/delta_star_results.csv'
        results_summary.to_csv(output_file, index=False)
        print(f"\nDetailed results saved to: {output_file}")

        # Print completion message
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE!")
        print("="*70)

    return delta_star, results_summary


# Execute the analysis when script is run directly (not imported as a module)
if __name__ == "__main__":
    # Run the complete analysis and store the results
    delta_star, results = main()
