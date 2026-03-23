"""
Beta=Infinity (SI Model) Network Simulation for Delta-Star Analysis

Simulates diffusion on C_{p,q/n} vs Random Regular networks under the
SI model (once infected, always infectious -- beta = infinity).

Key simplifications from beta=infinity:
  - No recovery: infected nodes stay infectious forever
  - No recovery_times tracking
  - State reduces to: susceptible or infected

Adoption rule: Step-function threshold (matches theoretical derivation)
  Pi(m) = 0   if m = 0
  Pi(m) = p1  if 1 <= m < i
  Pi(m) = p2  if m >= i
"""

import random
import math
import numpy as np
import networkx as nx


# ================================================================
# Network Generation
# ================================================================

def create_cycle_power_union_random(n, p, q, seed=None):
    """
    Create C_{p,q/n} = C_p union G_{n,q/n}.

    C_p: each node connects to all nodes within p-hop distance on cycle (2p edges)
    G_{n,q/n}: Erdos-Renyi with edge probability q/n (~q random edges per node)
    Expected degree: 2p + q
    """
    G = nx.Graph()
    G.add_nodes_from(range(n))
    rng = np.random.RandomState(seed)

    # Lattice edges (C_p)
    for i in range(n):
        for offset in range(1, p + 1):
            G.add_edge(i, (i - offset) % n)
            G.add_edge(i, (i + offset) % n)

    # Random edges (G_{n,q/n})
    edge_prob = q / n
    for i in range(n):
        for j in range(i + 1, n):
            if not G.has_edge(i, j):
                if rng.random() < edge_prob:
                    G.add_edge(i, j)

    return G


def create_random_regular(n, k, seed=None):
    """Create a random regular graph where all nodes have exactly degree k."""
    return nx.random_regular_graph(k, n, seed=seed)


# ================================================================
# Seeding Strategy
# ================================================================

def seed_strat_path(G, n, seeds):
    """
    Pick seeds that form a connected linear path.
    Works for both structured (Cpq) and random networks.
    """
    start = random.randint(0, n - 1)
    seed_list = [start]
    current = start
    visited = {start}

    attempts = 0
    max_attempts = n * 10

    while len(seed_list) < seeds and attempts < max_attempts:
        attempts += 1
        neighbors = [nb for nb in G.neighbors(current) if nb not in visited]

        if neighbors:
            current = random.choice(neighbors)
            seed_list.append(current)
            visited.add(current)
        else:
            found_next = False
            for seed_node in seed_list:
                neighbors = [nb for nb in G.neighbors(seed_node)
                             if nb not in visited]
                if neighbors:
                    current = random.choice(neighbors)
                    seed_list.append(current)
                    visited.add(current)
                    found_next = True
                    break
            if not found_next:
                remaining = [i for i in range(n) if i not in visited]
                if remaining:
                    current = random.choice(remaining)
                    seed_list.append(current)
                    visited.add(current)
                else:
                    break

    return seed_list[:seeds]


# ================================================================
# SI Simulation (beta = infinity)
# ================================================================

def simulate_si(G, n, threshold, p1, p2, seeds_list):
    """
    Simulate SI diffusion on network G with beta=infinity.

    Adoption rule: Step-function threshold (matches theory)
      Pi(m) = 0   if m = 0
      Pi(m) = p1  if 1 <= m < threshold
      Pi(m) = p2  if m >= threshold

    Parameters:
        G: networkx graph
        n: number of nodes
        threshold: adoption threshold (integer)
        p1: base adoption probability
        p2: reinforced adoption probability
        seeds_list: list of initial seed node indices

    Returns:
        cumulative_timeseries: list of cumulative infected counts [c_0, c_1, ...]
    """
    i = int(math.ceil(threshold))

    # Initialize: all susceptible except seeds
    infected = set(seeds_list)
    susceptible = set(range(n)) - infected
    cumulative_timeseries = [len(infected)]

    max_steps = n
    for t in range(1, max_steps):
        if not susceptible:
            break

        newly_infected = set()
        for s_node in susceptible:
            # Count infected neighbors (beta=inf: ALL infected are infectious)
            m = sum(1 for nb in G.neighbors(s_node) if nb in infected)

            if m == 0:
                continue

            # Step-function adoption probability
            if m < i:
                prob = p1
            else:
                prob = p2

            if random.random() < prob:
                newly_infected.add(s_node)

        if not newly_infected:
            break

        infected |= newly_infected
        susceptible -= newly_infected
        cumulative_timeseries.append(len(infected))

    return cumulative_timeseries


# ================================================================
# V(delta) and delta* computation
# ================================================================

def cumulative_to_new_infections(cumulative):
    """Convert cumulative counts to per-timestep new infections."""
    new_inf = [cumulative[0]]
    for t in range(1, len(cumulative)):
        new_inf.append(cumulative[t] - cumulative[t - 1])
    return new_inf


def calculate_V(new_infections, delta):
    """V(delta) = sum_{t=0}^{T} c_t * delta^t."""
    V = 0.0
    for t, c_t in enumerate(new_infections):
        V += c_t * (delta ** (t + 1))
    return V


def find_delta_star_from_V(V_cpq, V_random, delta_values, min_delta=0.0001):
    """
    Find delta* where V_cpq(delta) = V_random(delta).

    Returns:
        dict with 'case', 'delta_star', 'all_delta_stars'
    """
    diff = np.array(V_cpq) - np.array(V_random)

    # Find sign changes
    intersections = []
    for idx in range(len(diff) - 1):
        if diff[idx] * diff[idx + 1] < 0:
            x1, x2 = delta_values[idx], delta_values[idx + 1]
            y1, y2 = diff[idx], diff[idx + 1]
            ds = x1 - y1 * (x2 - x1) / (y2 - y1)
            if ds >= min_delta:
                intersections.append(ds)

    if intersections:
        return {
            'case': 'C',
            'delta_star': intersections[0],
            'all_delta_stars': intersections,
        }

    # No intersection: determine dominance
    if np.all(diff >= 0):
        return {'case': 'A', 'delta_star': None, 'all_delta_stars': []}
    elif np.all(diff <= 0):
        return {'case': 'B', 'delta_star': None, 'all_delta_stars': []}
    else:
        mean_diff = np.mean(diff)
        return {
            'case': 'A' if mean_diff > 0 else 'B',
            'delta_star': None,
            'all_delta_stars': [],
        }


# ================================================================
# Full pair analysis
# ================================================================

def analyze_pair(pair_id, n, k, p, q, threshold, p1, p2, delta_values):
    """
    Analyze a single network pair: C_{p,q/n} vs Random Regular.
    Uses beta=infinity (SI model).

    Returns:
        dict with pair_id, case, delta_star, etc.
    """
    # Generate networks
    G_cpq = create_cycle_power_union_random(n, p, q, seed=pair_id)
    G_random = create_random_regular(n, k, seed=pair_id)

    seeds = int(math.ceil(threshold))

    # Seed both networks with same random seed for fairness
    random.seed(pair_id)
    seeds_cpq = seed_strat_path(G_cpq, n, seeds)

    random.seed(pair_id)
    seeds_random = seed_strat_path(G_random, n, seeds)

    # Simulate (beta=infinity)
    ts_cpq = simulate_si(G_cpq, n, threshold, p1, p2, seeds_cpq)
    ts_random = simulate_si(G_random, n, threshold, p1, p2, seeds_random)

    # Convert to new infections
    new_cpq = cumulative_to_new_infections(ts_cpq)
    new_random = cumulative_to_new_infections(ts_random)

    # Compute V for each delta
    V_cpq = [calculate_V(new_cpq, d) for d in delta_values]
    V_random = [calculate_V(new_random, d) for d in delta_values]

    # Find delta*
    ds_result = find_delta_star_from_V(V_cpq, V_random, delta_values)

    return {
        'pair_id': pair_id,
        'case': ds_result['case'],
        'delta_star': ds_result['delta_star'],
        'all_delta_stars': ds_result['all_delta_stars'],
        'total_infected_cpq': ts_cpq[-1],
        'total_infected_random': ts_random[-1],
        'timesteps_cpq': len(ts_cpq),
        'timesteps_random': len(ts_random),
    }


# ================================================================
# CLI
# ================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='SI simulation (beta=inf): C_{p,q/n} vs Random Regular')
    parser.add_argument('--p', type=int, required=True)
    parser.add_argument('--q', type=float, default=None)
    parser.add_argument('--p1', type=float, default=0.1)
    parser.add_argument('--p2', type=float, default=1.0)
    parser.add_argument('--threshold', type=float, default=2.0)
    parser.add_argument('--N', type=int, default=2000)
    parser.add_argument('--n-pairs', type=int, default=10)
    args = parser.parse_args()

    q = args.q if args.q is not None else float(8 - 2 * args.p)
    k = 8
    delta_values = np.linspace(0, 1, 201)

    print("=" * 60)
    print(f"SI Simulation: C_{{{args.p},{q:.0f}/n}} vs Random Regular")
    print(f"  beta = infinity, N={args.N}")
    print(f"  p1={args.p1}, p2={args.p2}, threshold={args.threshold}")
    print(f"  pairs={args.n_pairs}")
    print("=" * 60)

    cases = {'A': 0, 'B': 0, 'C': 0}
    delta_stars = []

    for pid in range(args.n_pairs):
        print(f"\r  Pair {pid + 1}/{args.n_pairs}...", end='', flush=True)
        result = analyze_pair(pid, args.N, k, args.p, q,
                              args.threshold, args.p1, args.p2, delta_values)
        cases[result['case']] = cases.get(result['case'], 0) + 1
        if result['delta_star'] is not None:
            delta_stars.append(result['delta_star'])

    print(f"\n\nResults:")
    print(f"  Case A (C dom): {cases.get('A', 0)}")
    print(f"  Case B (R dom): {cases.get('B', 0)}")
    print(f"  Case C (cross): {cases.get('C', 0)}")
    if delta_stars:
        print(f"  delta* mean: {np.mean(delta_stars):.4f}")
        print(f"  delta* std:  {np.std(delta_stars):.4f}")
    print("=" * 60)
