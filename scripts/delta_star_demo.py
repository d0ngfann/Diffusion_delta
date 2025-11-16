# δ* (delta-star) 데모 스크립트 - 모의 데이터로 개념 시연
"""
This demo script demonstrates the δ* concept using mock simulation data.
No external dependencies required - runs with standard Python only.

This allows you to understand the implementation without setting up
the full simulation environment.
"""

import random
random.seed(42)

###################
# Core Functions (same as in delta_star_analysis.py)

def cumulative_to_new_infections(cumulative_series):
    """
    Convert cumulative adoption counts to per-timestep new infections.

    Example:
        Input:  [5, 8, 12, 15, 15, 15] (cumulative)
        Output: [5, 3, 4, 3, 0, 0]     (new per timestep)
    """
    new_infections = [cumulative_series[0]]
    for t in range(1, len(cumulative_series)):
        new_at_t = cumulative_series[t] - cumulative_series[t-1]
        new_infections.append(new_at_t)
    return new_infections


def calculate_V(new_infections, delta):
    """
    Calculate discounted cumulative activation size V.

    Formula: V = c₁×δ¹ + c₂×δ² + c₃×δ³ + ... + cₜ×δᵗ

    Parameters:
        new_infections: List of newly infected at each timestep
        delta: Temporal discount factor (0 to 1)

    Returns:
        V: Discounted cumulative activation size
    """
    V = 0.0
    for t in range(len(new_infections)):
        V += new_infections[t] * (delta ** (t + 1))
    return V


###################
# Mock Data Generation

def generate_mock_timeseries(network_type='clustered', n=2000, trials=20):
    """
    Generate mock time series data mimicking real simulation results.

    Clustered networks: Fast initial spread, plateaus early
    Random networks: Slower start, continues spreading longer

    Returns:
        List of cumulative time series (one per trial)
    """
    all_trials = []

    for trial in range(trials):
        if network_type == 'clustered':
            # Clustered: Quick burst early due to reinforcement
            max_spread = random.randint(int(0.6*n), int(0.8*n))
            # Rapid early adoption
            t0 = 5
            t1 = t0 + random.randint(80, 120)
            t2 = t1 + random.randint(150, 200)
            t3 = t2 + random.randint(100, 150)
            # Then plateau
            timeseries = [t0, t1, t2, t3] + [max_spread] * 10

        else:  # random network
            # Random: Slower start but continues longer
            max_spread = random.randint(int(0.5*n), int(0.7*n))
            # Slower but sustained growth
            t0 = 5
            t1 = t0 + random.randint(40, 60)
            t2 = t1 + random.randint(80, 120)
            t3 = t2 + random.randint(120, 180)
            t4 = t3 + random.randint(100, 150)
            t5 = t4 + random.randint(80, 120)
            # Longer tail
            timeseries = [t0, t1, t2, t3, t4, t5] + [max_spread] * 8

        all_trials.append(timeseries)

    return all_trials


def compute_V_for_trials(trials_data, delta_values):
    """
    Compute V for all trials across different delta values.

    Returns:
        means: List of mean V values for each delta
        ses: List of standard errors for each delta
    """
    n_trials = len(trials_data)
    n_deltas = len(delta_values)

    # Matrix: rows = trials, cols = delta values
    V_matrix = []

    for trial_cumulative in trials_data:
        trial_V = []
        new_infections = cumulative_to_new_infections(trial_cumulative)

        for delta in delta_values:
            V = calculate_V(new_infections, delta)
            trial_V.append(V)

        V_matrix.append(trial_V)

    # Calculate means and standard errors
    means = []
    ses = []

    for delta_idx in range(n_deltas):
        values = [V_matrix[trial][delta_idx] for trial in range(n_trials)]
        mean_V = sum(values) / len(values)
        variance = sum((v - mean_V)**2 for v in values) / len(values)
        se = (variance ** 0.5) / (len(values) ** 0.5)

        means.append(mean_V)
        ses.append(se)

    return means, ses


def find_delta_star(delta_values, V_clustered_mean, V_random_mean):
    """
    Find δ* where clustered and random networks have equal performance.
    """
    diff = [V_random_mean[i] - V_clustered_mean[i] for i in range(len(delta_values))]

    # Find sign change
    for i in range(len(diff) - 1):
        if (diff[i] < 0 and diff[i+1] > 0) or (diff[i] > 0 and diff[i+1] < 0):
            # Linear interpolation
            x1, x2 = delta_values[i], delta_values[i + 1]
            y1, y2 = diff[i], diff[i + 1]
            delta_star = x1 - y1 * (x2 - x1) / (y2 - y1)
            return delta_star

    return None


def create_ascii_plot(delta_values, V_clustered, V_random, delta_star=None):
    """
    Create a simple ASCII plot to visualize results.
    """
    print("\n" + "="*70)
    print("VISUALIZATION: V vs δ")
    print("="*70)

    # Normalize values for display
    all_vals = V_clustered + V_random
    min_val = min(all_vals)
    max_val = max(all_vals)
    scale = 50 / (max_val - min_val) if max_val > min_val else 1

    print(f"\n{'δ':>6} | {'V_clustered':>12} | {'V_random':>12} | Visual Comparison")
    print("-" * 70)

    for i, delta in enumerate(delta_values):
        v_c = V_clustered[i]
        v_r = V_random[i]

        # Normalize to 0-50 range for ASCII plot
        c_pos = int((v_c - min_val) * scale)
        r_pos = int((v_r - min_val) * scale)

        # Create visual bars
        bar = [' '] * 52
        bar[c_pos] = 'C'
        bar[r_pos] = 'R'
        visual = ''.join(bar)

        winner = "C" if v_c > v_r else "R"
        mark = " *" if delta_star and abs(delta - delta_star) < 0.05 else ""

        print(f"{delta:>6.2f} | {v_c:>12.2f} | {v_r:>12.2f} | {visual} {winner}{mark}")

    print("\n  Legend: C = Clustered, R = Random, * = near δ*")
    print("="*70)


###################
# Main Demo

def main():
    """
    Run the complete δ* demonstration with mock data.
    """
    print("="*70)
    print("δ* (DELTA-STAR) DEMONSTRATION")
    print("Using mock simulation data to demonstrate the concept")
    print("="*70)

    # Parameters
    n = 2000  # network size
    trials = 20
    delta_values = [0.5, 0.7, 0.9, 0.95, 0.99]

    print("\nParameters:")
    print(f"  Network size: {n} nodes")
    print(f"  Trials per network: {trials}")
    print(f"  Delta values tested: {delta_values}")

    # Generate mock data
    print("\n" + "-"*70)
    print("Generating mock simulation data...")
    print("-"*70)

    clustered_data = generate_mock_timeseries('clustered', n, trials)
    random_data = generate_mock_timeseries('random', n, trials)

    print(f"  Clustered network: {len(clustered_data)} trials generated")
    print(f"  Random network: {len(random_data)} trials generated")

    # Show example time series
    print("\nExample time series (Trial 0):")
    print(f"  Clustered (cumulative): {clustered_data[0]}")
    print(f"  Random (cumulative):    {random_data[0]}")

    c_new = cumulative_to_new_infections(clustered_data[0])
    r_new = cumulative_to_new_infections(random_data[0])
    print(f"\n  Clustered (new per step): {c_new}")
    print(f"  Random (new per step):    {r_new}")

    # Calculate V for different deltas
    print("\n" + "-"*70)
    print("Computing V (discounted cumulative activation)...")
    print("-"*70)

    V_clustered_mean, V_clustered_se = compute_V_for_trials(clustered_data, delta_values)
    V_random_mean, V_random_se = compute_V_for_trials(random_data, delta_values)

    # Find delta*
    print("\n" + "-"*70)
    print("Finding δ* (critical threshold)...")
    print("-"*70)

    delta_star = find_delta_star(delta_values, V_clustered_mean, V_random_mean)

    if delta_star:
        print(f"\n*** δ* = {delta_star:.4f} ***\n")
        print("Interpretation:")
        print(f"  • When δ < {delta_star:.3f}: Clustered networks WIN")
        print(f"    → Early spread matters more (reinforcement effect)")
        print(f"  • When δ > {delta_star:.3f}: Random networks WIN")
        print(f"    → Long-term reach matters more")
    else:
        print("\nNo intersection found in the tested range")

    # Results table
    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)
    print(f"\n{'δ':>8} | {'V_clustered':>15} | {'V_random':>15} | {'Difference':>15} | {'Winner':>10}")
    print("-" * 80)

    for i, delta in enumerate(delta_values):
        diff = V_random_mean[i] - V_clustered_mean[i]
        winner = "Random" if diff > 0 else "Clustered"
        print(f"{delta:>8.2f} | {V_clustered_mean[i]:>15.2f} | {V_random_mean[i]:>15.2f} | "
              f"{diff:>+15.2f} | {winner:>10}")

    print("="*80)

    # ASCII visualization
    create_ascii_plot(delta_values, V_clustered_mean, V_random_mean, delta_star)

    # Example V calculation
    print("\n" + "="*70)
    print("EXAMPLE: How V is calculated")
    print("="*70)

    example_new = c_new[:5]  # First 5 timesteps
    example_delta = 0.9

    print(f"\nGiven new infections: {example_new}")
    print(f"And delta = {example_delta}")
    print(f"\nV = ", end="")

    terms = []
    for t in range(len(example_new)):
        term = example_new[t] * (example_delta ** (t+1))
        terms.append(f"{example_new[t]}×{example_delta}^{t+1}")

    print(" + ".join(terms))

    V_example = calculate_V(example_new, example_delta)
    print(f"  = {V_example:.2f}")

    print("\n" + "="*70)
    print("DEMO COMPLETE!")
    print("="*70)
    print("\nNext steps:")
    print("  1. Review delta_star_analysis.py for the full implementation")
    print("  2. Set up the conda environment (see README)")
    print("  3. Run real simulations with: python delta_star_analysis.py")


if __name__ == "__main__":
    main()
