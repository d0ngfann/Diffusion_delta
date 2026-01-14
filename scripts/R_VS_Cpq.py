"""
Pair-wise Delta* Analysis: Random Regular vs Cycle-Power-Union-Random (C_{p,q/n})

This script implements a pair-wise comparison approach to find delta* (critical temporal
discount factor) for complex contagion simulations.

Key Methodology:
- Generate 100 independent network pairs (Random Regular vs C_{p,q/n})
- Each pair uses independent random seeds
- For each pair, determine:
  * Case A: C_{p,q/n} always dominates (in valid region below 60% threshold)
  * Case B: Random always dominates
  * Case C: Intersection exists (delta* found)
  * Case Invalid: No valid comparison region (rare)

Network Structures:
1. Random Regular Graph: All nodes have exactly degree k, randomly connected
2. C_{p,q/n}: Cycle-power-union-random graph
   - C_p: Each node connected to all nodes within p-hop distance on cycle (2p edges)
   - G_{n,q/n}: Erdős-Rényi random graph with edge probability q/n (expected q edges)
   - Expected total degree: 2p + q ≈ k

The 60% threshold (n*0.6) acts as a "goal line" - we only analyze behavior before
reaching this saturation level, as that's where interesting dynamics occur.
"""

import os
import math
import random
import numpy as np
import pandas as pd
import networkx as nx
from DH_helper_functions import (
    sirsic_wmem, run_simulation,
    seed_strat_one_nbrs, seed_strat_random, seed_strat_adj, seed_strat_path
)


###################
# Network Generation Functions

def create_cycle_power_union_random(n, p, q, seed=None):
    """
    Create a Cycle-Power-Union-Random graph C_{p,q/n}.

    The graph is defined as: C_{p,q/n} := C_p ∪ G_{n,q/n}

    - C_p: Cycle-power-p graph (each node connected to all nodes within p-hop distance)
    - G_{n,q/n}: Erdős-Rényi random graph with edge probability q/n

    Parameters:
        n: int
            Number of nodes
        p: int
            Cycle-power parameter (each node connects to all nodes within p hops)
            Results in exactly 2p lattice edges per node
        q: float
            Random component parameter (edge probability = q/n)
            Results in expected q random edges per node
        seed: int or None
            Random seed for reproducibility

    Returns:
        G: NetworkX graph
            The C_{p,q/n} graph with expected degree ≈ 2p + q

    Example:
        n=1000, p=2, q=4:
        - Each node has 4 lattice neighbors (2 on each side, within 2-hop distance)
        - Each node has ~4 random edges (expected)
        - Total expected degree: 8
    """
    # Initialize graph
    G = nx.Graph()
    G.add_nodes_from(range(n))

    rng = np.random.RandomState(seed)

    # 1. Add cycle-power-p edges (C_p)
    # Each node i connects to all nodes within p-hop distance on the cycle
    for i in range(n):
        for offset in range(1, p + 1):
            # Connect to left neighbor at distance 'offset'
            j_left = (i - offset) % n
            G.add_edge(i, j_left)

            # Connect to right neighbor at distance 'offset'
            j_right = (i + offset) % n
            G.add_edge(i, j_right)

    # 2. Add Erdős-Rényi random edges (G_{n,q/n})
    # Edge probability = q/n for each pair
    edge_prob = q / n

    for i in range(n):
        for j in range(i + 1, n):
            # Skip if edge already exists from cycle-power
            if not G.has_edge(i, j):
                # Add edge with probability q/n
                if rng.random() < edge_prob:
                    G.add_edge(i, j)

    return G


def create_random_regular(n, k, seed=None):
    """
    Create a random regular graph where all nodes have exactly degree k.

    Parameters:
        n: int
            Number of nodes
        k: int
            Degree (must be even for regular graph construction)
        seed: int or None
            Random seed for reproducibility

    Returns:
        G: NetworkX graph
            Random regular graph with all nodes having degree k

    Raises:
        ValueError: If k*n is odd (impossible to construct)
    """
    return nx.random_regular_graph(k, n, seed=seed)


###################
# Core Functions (Reused from DH_delta_star_pairwise.py)

def cumulative_to_new_infections(cumulative_series):
    """
    Convert cumulative adoption counts to per-timestep new infections.

    Parameters:
        cumulative_series: list of integers
            Cumulative number of adopters at each time step

    Returns:
        new_infections: list of integers
            Number of newly infected nodes at each time step
    """
    new_infections = [cumulative_series[0]]
    for t in range(1, len(cumulative_series)):
        new_at_t = cumulative_series[t] - cumulative_series[t-1]
        new_infections.append(new_at_t)
    return new_infections


def calculate_V(new_infections, delta):
    """
    Calculate the discounted cumulative activation size V.

    Formula: V = Σ(cₜ × δᵗ) for t = 1, 2, 3, ...

    Parameters:
        new_infections: list of integers
            Number of newly activated nodes at each time step
        delta: float between 0 and 1
            Temporal discount factor

    Returns:
        V: float
            Discounted cumulative activation size
    """
    V = 0.0
    for t in range(len(new_infections)):
        V += new_infections[t] * (delta ** (t + 1))
    return V


###################
# Simulation Functions

def simulate_single_network(G, network_type, k, n, thrshld, p1, p2, beta, seeds, pair_id):
    """
    Run a single simulation on one network with continuous thrshld and beta.
    This version implements a custom SIR-style simulation loop to handle
    stochastic, heterogeneous node parameters.
    """
    # Use a separate random generator for simulation to not interfere with network generation
    sim_rng = random.Random(pair_id)

    # 1. Assign heterogeneous integer thresholds to each node
    node_thresholds = {}
    floor_thrshld = math.floor(thrshld)
    prob_thrshld = thrshld - floor_thrshld
    for u in G.nodes():
        extra = 1 if sim_rng.random() < prob_thrshld else 0
        node_thresholds[u] = floor_thrshld + extra

    # 2. Initialize simulation state (S=Susceptible, I=Infected, R=Recovered/Adopted)
    states = {u: 'S' for u in G.nodes()}
    recovery_times = {}

    # 3. Seeding: Select 'seeds' nodes to start as infected
    # Use path seeding: ensures seeds form a connected linear path in BOTH networks
    # Set random seed for reproducibility
    random.seed(pair_id)
    initial_seeds = seed_strat_path(G, n, int(seeds))

    floor_beta = math.floor(beta)
    prob_beta = beta - floor_beta

    for u in initial_seeds:
        states[u] = 'I'
        # Determine influence duration and set recovery time
        duration = floor_beta + (1 if sim_rng.random() < prob_beta else 0)
        recovery_times[u] = 1 + duration  # Recover at the start of time t = 1 + duration

    # 4. Run simulation loop
    # Start with t=0 state (initial seeds)
    full_timeseries = [len(initial_seeds)]
    
    max_steps = n  # Stop simulation after n steps as a safeguard
    for t in range(1, max_steps):
        
        # A. Recovery Step: Nodes scheduled to recover at time t transition from I to R
        recovering_nodes = {u for u, rec_time in recovery_times.items() if rec_time == t}
        for u_rec in recovering_nodes:
            if states[u_rec] == 'I':
                states[u_rec] = 'R'

        # B. Infection Step
        newly_infected = set()
        susceptible_nodes = [u for u, state in states.items() if state == 'S']
        
        for s_node in susceptible_nodes:
            # Count infected neighbors
            infected_neighbors = [v for v in G.neighbors(s_node) if states[v] == 'I']
            num_infected_neighbors = len(infected_neighbors)

            if num_infected_neighbors == 0:
                continue

            # Determine adoption probability based on threshold
            prob = p2 if num_infected_neighbors >= node_thresholds[s_node] else p1
            
            # Attempt to infect the node
            if sim_rng.random() < prob:
                newly_infected.add(s_node)
        
        # C. Update states for newly infected nodes
        for u_inf in newly_infected:
            states[u_inf] = 'I'
            # Determine recovery time for the newly infected node
            duration = floor_beta + (1 if sim_rng.random() < prob_beta else 0)
            recovery_times[u_inf] = t + duration
        
        # D. Record statistics for this timestep
        current_adopted_count = len([u for u, state in states.items() if state in ['I', 'R']])
        full_timeseries.append(current_adopted_count)

        # E. Check for simulation end condition (no more infected nodes)
        if not any(state == 'I' for state in states.values()):
            break
            
    return full_timeseries


def compute_V_for_single_trial(full_timeseries, delta_values):
    """
    Compute V values across different delta values for a single trial.

    Parameters:
        full_timeseries: list
            Cumulative adoption time series
        delta_values: list or array
            Delta values to test

    Returns:
        V_array: numpy array
            V values for each delta
    """
    # Convert cumulative to new infections
    new_infections = cumulative_to_new_infections(full_timeseries)

    # Calculate V for each delta
    V_array = np.zeros(len(delta_values))
    for i, delta in enumerate(delta_values):
        V_array[i] = calculate_V(new_infections, delta)

    return V_array


def find_intersections_in_valid_region(delta_values, V_cpq, V_random,
                                       valid_mask, min_delta=0.0001):
    """
    Find all intersection points (delta*) within the valid region.

    Parameters:
        delta_values: array
            Delta values tested
        V_cpq: array
            V values for C_{p,q/n} network
        V_random: array
            V values for random network
        valid_mask: boolean array
            Mask indicating valid region (both below threshold)
        min_delta: float
            Minimum delta to consider (exclude delta ≈ 0)

    Returns:
        intersections: list of dicts
            Each dict contains:
                - 'delta_star': float
                - 'transition': str ('cpq_to_random' or 'random_to_cpq')
    """
    # Apply valid mask
    delta_valid = delta_values[valid_mask]
    V_cpq_valid = V_cpq[valid_mask]
    V_random_valid = V_random[valid_mask]

    # Check if we have enough points
    if len(delta_valid) < 2:
        return []

    # Calculate difference: positive = random better
    diff = V_random_valid - V_cpq_valid

    # Find sign changes
    sign_changes = np.where(np.diff(np.sign(diff)))[0]

    intersections = []
    for idx in sign_changes:
        # Linear interpolation
        x1, x2 = delta_valid[idx], delta_valid[idx + 1]
        y1, y2 = diff[idx], diff[idx + 1]

        # Find zero crossing
        delta_star = x1 - y1 * (x2 - x1) / (y2 - y1)

        # Check minimum delta condition
        if delta_star < min_delta:
            continue

        # Determine transition direction
        if y1 < 0 and y2 > 0:
            transition = 'cpq_to_random'
        else:
            transition = 'random_to_cpq'

        intersections.append({
            'delta_star': delta_star,
            'transition': transition
        })

    return intersections


def analyze_pair(pair_id, n, k, p, q, thrshld, p1, p2, beta, seeds, delta_values):
    """
    Analyze a single network pair (C_{p,q/n} vs Random Regular).

    Parameters:
        pair_id: int
            Identifier for this pair (used for random seed)
        n: int
            Number of nodes
        k: int
            Target degree for random graph
        p: int
            Cycle-power parameter for C_{p,q/n}
        q: float
            Random component parameter for C_{p,q/n}
        thrshld, p1, p2, beta, seeds: simulation parameters
        delta_values: array
            Delta values to test

    Returns:
        result: dict
            Analysis result with keys:
                - 'pair_id'
                - 'case': 'A', 'B', 'C', or 'Invalid'
                - 'delta_star': float or None
                - 'dominant_network': str or None
                - 'all_delta_stars': list
                - 'multiple_intersections': bool
                - 'reason': str (for Invalid cases)
    """
    # Generate C_{p,q/n} network
    G_cpq = create_cycle_power_union_random(n, p, q, seed=pair_id)

    # Generate Random Regular network
    G_random = create_random_regular(n, k, seed=pair_id)

    # Simulate C_{p,q/n} network
    timeseries_cpq = simulate_single_network(
        G_cpq, 'Cpq', k=k, n=n, thrshld=thrshld,
        p1=p1, p2=p2, beta=beta, seeds=seeds, pair_id=pair_id
    )

    # Simulate Random network
    timeseries_random = simulate_single_network(
        G_random, 'Random', k=k, n=n, thrshld=thrshld,
        p1=p1, p2=p2, beta=beta, seeds=seeds, pair_id=pair_id
    )

    # Compute V values for both networks
    V_cpq = compute_V_for_single_trial(timeseries_cpq, delta_values)
    V_random = compute_V_for_single_trial(timeseries_random, delta_values)

    # Define threshold (goal line)
    threshold = n * 0.6

    # Define valid region: both below threshold
    valid_mask = (V_cpq <= threshold) & (V_random <= threshold)

    # Check if we have a valid comparison region
    if not np.any(valid_mask):
        return {
            'pair_id': pair_id,
            'case': 'Invalid',
            'delta_star': None,
            'dominant_network': None,
            'all_delta_stars': [],
            'multiple_intersections': False,
            'reason': 'All values exceed 60% threshold'
        }

    # Find intersections within valid region
    intersections = find_intersections_in_valid_region(
        delta_values, V_cpq, V_random, valid_mask
    )

    # Case C: Intersection(s) found
    if len(intersections) > 0:
        return {
            'pair_id': pair_id,
            'case': 'C',
            'delta_star': intersections[0]['delta_star'],  # Primary = first
            'dominant_network': None,
            'all_delta_stars': [i['delta_star'] for i in intersections],
            'multiple_intersections': len(intersections) > 1,
            'reason': None
        }

    # No intersection: determine which network dominates
    # Extract valid portions
    V_cpq_valid = V_cpq[valid_mask]
    V_random_valid = V_random[valid_mask]

    # Check which network is consistently above the other
    # "Above" means better performance (higher V)
    if np.all(V_cpq_valid >= V_random_valid):
        # C_{p,q/n} always dominates
        return {
            'pair_id': pair_id,
            'case': 'A',
            'delta_star': None,
            'dominant_network': 'Cpq',
            'all_delta_stars': [],
            'multiple_intersections': False,
            'reason': None
        }
    elif np.all(V_random_valid >= V_cpq_valid):
        # Random always dominates
        return {
            'pair_id': pair_id,
            'case': 'B',
            'delta_star': None,
            'dominant_network': 'Random',
            'all_delta_stars': [],
            'multiple_intersections': False,
            'reason': None
        }
    else:
        # This shouldn't happen: curves cross but no intersection found
        # Fallback: check mean difference
        mean_diff = np.mean(V_cpq_valid - V_random_valid)
        if mean_diff > 0:
            dominant = 'Cpq'
            case = 'A'
        else:
            dominant = 'Random'
            case = 'B'

        return {
            'pair_id': pair_id,
            'case': case,
            'delta_star': None,
            'dominant_network': dominant,
            'all_delta_stars': [],
            'multiple_intersections': False,
            'reason': 'Determined by mean (edge case)'
        }


###################
# Main Function

def main(n=1000, k=8, p=2, q=4, thrshld=2, p1=0.1, p2=1.0, beta=5, seeds=None,
         n_pairs=100, delta_values=None, verbose=True):
    """
    Main function to run pair-wise delta* analysis: Random Regular vs C_{p,q/n}.

    Parameters:
        n: int
            Number of nodes (default: 1000)
        k: int
            Target degree for random graph (default: 8)
        p: int
            Cycle-power parameter for C_{p,q/n} (default: 2)
            Results in 2p lattice edges per node
        q: float
            Random component parameter for C_{p,q/n} (default: 4)
            Results in expected q random edges per node
        thrshld: int
            Social reinforcement threshold (default: 2)
        p1: float
            Base adoption probability (default: 0.1)
        p2: float
            Reinforced adoption probability (default: 1.0)
        beta: int
            Time of influence (default: 5)
        seeds: int or None
            Number of initial seeds. If None, uses thrshld (default: None)
        n_pairs: int
            Number of network pairs to generate (default: 100)
        delta_values: list or None
            Delta values to test. If None, uses np.linspace(0, 1, 201)
        verbose: bool
            Print progress (default: True)

    Returns:
        results_df: DataFrame
            Results for all pairs

    Note:
        Expected degree for C_{p,q/n}: 2p + q
        If you want C_{p,q/n} to have similar degree as Random, set 2p + q ≈ k
    """
    # Set defaults
    if seeds is None:
        seeds = math.ceil(thrshld)
    if delta_values is None:
        delta_values = np.linspace(0, 1, 201)

    # Calculate expected degrees
    expected_cpq_degree = 2 * p + q

    # Print header
    if verbose:
        print("="*70)
        print("PAIR-WISE DELTA* ANALYSIS: Random Regular vs C_{p,q/n}")
        print("="*70)
        print(f"\nNetwork Parameters:")
        print(f"  Number of nodes (n): {n}")
        print(f"  Random Regular Graph:")
        print(f"    - Degree (k): {k} (exact)")
        print(f"  C_{{p,q/n}} Graph:")
        print(f"    - Cycle-power (p): {p} → {2*p} lattice edges per node")
        print(f"    - Random component (q): {q} → ~{q} random edges per node (expected)")
        print(f"    - Expected total degree: 2p + q = {expected_cpq_degree}")

        if abs(expected_cpq_degree - k) > 0.5:
            print(f"\n  ⚠ WARNING: Degree mismatch!")
            print(f"    Random: k={k} (exact)")
            print(f"    C_{{p,q/n}}: ~{expected_cpq_degree} (expected)")
            print(f"    Consider adjusting p and q so that 2p + q ≈ {k}")

        print(f"\nDiffusion Parameters:")
        print(f"  Threshold (i): {thrshld}")
        print(f"  Base adoption (p1): {p1}")
        print(f"  Reinforced adoption (p2): {p2}")
        print(f"  Time of influence (beta): {beta}")
        print(f"  Initial seeds: {seeds}")

        print(f"\nExperiment Parameters:")
        print(f"  Number of pairs: {n_pairs}")
        print(f"  Delta values: {len(delta_values)} points from 0 to 1")
        print(f"  Threshold (goal line): {n * 0.6:.0f} nodes (60% of {n})")
        print(f"\n{'='*70}\n")

    # Run analysis for each pair
    results = []
    for pair_id in range(n_pairs):
        if verbose:
            print(f"Processing pair {pair_id + 1}/{n_pairs}...", end='\r')

        result = analyze_pair(
            pair_id=pair_id,
            n=n, k=k, p=p, q=q,
            thrshld=thrshld, p1=p1, p2=p2,
            beta=beta, seeds=seeds, delta_values=delta_values
        )
        results.append(result)

    if verbose:
        print(f"\nCompleted {n_pairs} pairs!")

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Print summary
    if verbose:
        print(f"\n{'='*70}")
        print("SUMMARY")
        print("="*70)

        case_counts = results_df['case'].value_counts()
        print(f"\nCase distribution:")
        for case, count in case_counts.items():
            pct = count / n_pairs * 100
            if case == 'A':
                print(f"  Case A (C_{{p,q/n}} dominates): {count} ({pct:.1f}%)")
            elif case == 'B':
                print(f"  Case B (Random dominates): {count} ({pct:.1f}%)")
            elif case == 'C':
                print(f"  Case C (Intersection exists): {count} ({pct:.1f}%)")
            elif case == 'Invalid':
                print(f"  Case Invalid: {count} ({pct:.1f}%)")

        # Delta* statistics (for Case C only)
        delta_star_values = results_df[results_df['case'] == 'C']['delta_star'].dropna()
        if len(delta_star_values) > 0:
            print(f"\nDelta* statistics (Case C only, n={len(delta_star_values)}):")
            print(f"  Mean: {delta_star_values.mean():.4f}")
            print(f"  Std: {delta_star_values.std():.4f}")
            print(f"  SE: {delta_star_values.sem():.4f}")
            print(f"  Median: {delta_star_values.median():.4f}")
            print(f"  Min: {delta_star_values.min():.4f}")
            print(f"  Max: {delta_star_values.max():.4f}")

            # 95% CI (if enough samples)
            if len(delta_star_values) >= 2:
                from scipy import stats
                ci = stats.t.interval(0.95, len(delta_star_values)-1,
                                     loc=delta_star_values.mean(),
                                     scale=delta_star_values.sem())
                print(f"  95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]")
        else:
            print(f"\nNo delta* values found (no Case C pairs)")

        # Multiple intersections
        multiple_count = results_df['multiple_intersections'].sum()
        if multiple_count > 0:
            print(f"\n⚠ Warning: {multiple_count} pairs had multiple intersections (non-monotonic)")

        print(f"\n{'='*70}\n")

    # Save to CSV
    os.makedirs('../output_data', exist_ok=True)
    output_file = f'../output_data/R_VS_Cpq_n{n}_k{k}_p{p}_q{q}_results.csv'
    results_df.to_csv(output_file, index=False)
    if verbose:
        print(f"Results saved to: {output_file}\n")

    return results_df


if __name__ == "__main__":
    # Run the pair-wise analysis
    # Customize parameters here:
    results = main(
        n=2000,      # Number of nodes
        k=8,         # Degree for Random Regular graph
        p=2,         # Cycle-power parameter (2p lattice edges)
        q=4,         # Random component (expected q random edges)
        thrshld=2,   # Social reinforcement threshold
        p1=0.1,      # Base adoption probability
        p2=1.0,      # Reinforced adoption probability
        beta=5,      # Time of influence
        seeds=None,  # Initial seeds (None = use thrshld)
        n_pairs=100, # Number of network pairs
        delta_values=None,  # Delta values to test (None = 201 points from 0 to 1)
        verbose=True
    )
