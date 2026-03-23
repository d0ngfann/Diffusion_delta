"""
Simulation Sweep: p2 vs threshold for beta=infinity (SI model)

Runs network simulations for C_{p,q/n} vs Random Regular across a grid
of (p2, threshold) values. Beta is fixed at infinity (SI model).

Usage:
    cd scripts/beta_infinite_delta_star
    python sweep_simulation/sweep_p2_vs_thrshld.py --p 1 --q 6 --n-pairs 500
    python sweep_simulation/sweep_p2_vs_thrshld.py --p 2 --q 4 --n-pairs 500
    python sweep_simulation/sweep_p2_vs_thrshld.py --p 3 --q 2 --n-pairs 500
"""

import os
import sys
import time
import math
import argparse
import itertools
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.simulation import analyze_pair

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    class tqdm:
        def __init__(self, iterable, desc="", unit="", total=None):
            self.iterable = iterable
        def __iter__(self):
            return iter(self.iterable)
        @staticmethod
        def write(msg):
            print(msg)


def generate_parameter_combinations(p, q):
    """Generate all (p2, threshold) combinations."""
    n_values = [2000]
    k_values = [8]
    p1_values = [0.1]
    p2_values = np.arange(0.50, 1.01, 0.01)
    thrshld_values = np.arange(2.0, 5.1, 0.1)

    valid_combos = []
    combo_id = 1
    for n, k, thrshld, p1, p2 in itertools.product(
            n_values, k_values, thrshld_values, p1_values, p2_values):
        if p1 <= p2:
            valid_combos.append({
                'combination_id': combo_id,
                'n': int(n), 'k': int(k),
                'p': int(p), 'q': float(q),
                'thrshld': round(float(thrshld), 2),
                'p1': float(p1), 'p2': round(float(p2), 4),
                'seeds': int(math.ceil(thrshld)),
            })
            combo_id += 1
    return valid_combos


def run_single_combination(params, n_pairs=100, delta_values=None):
    """Run simulation for one (p2, threshold) combination."""
    if delta_values is None:
        delta_values = np.linspace(0, 1, 201)
    try:
        pair_results = []
        for pair_id in range(n_pairs):
            result = analyze_pair(
                pair_id=pair_id,
                n=params['n'], k=params['k'],
                p=params['p'], q=params['q'],
                threshold=params['thrshld'],
                p1=params['p1'], p2=params['p2'],
                delta_values=delta_values,
            )
            pair_results.append({
                'combination_id': params['combination_id'],
                'pair_id': pair_id,
                **params, **result,
            })

        case_counts = pd.Series(
            [r['case'] for r in pair_results]).value_counts()
        delta_stars = [
            r['delta_star'] for r in pair_results
            if r.get('case') == 'C' and r.get('delta_star') is not None
        ]

        summary = {
            **params,
            'case_A_count': case_counts.get('A', 0),
            'case_B_count': case_counts.get('B', 0),
            'case_C_count': case_counts.get('C', 0),
            'case_Invalid_count': case_counts.get('Invalid', 0),
            'delta_star_mean': np.mean(delta_stars) if delta_stars else None,
            'delta_star_std': (np.std(delta_stars, ddof=1)
                               if len(delta_stars) > 1 else None),
            'status': 'success',
        }
        return summary, pair_results
    except Exception as e:
        return {**params, 'status': 'error', 'error_message': str(e)}, []


def run_parameter_sweep(combinations, n_pairs=100, n_workers=8):
    """Run sweep across all parameter combinations in parallel."""
    n_combos = len(combinations)
    summary_results, full_results = [], []
    print(f"\nStarting sweep: {n_combos} combinations, "
          f"{n_pairs} pairs each, {n_workers} workers...")
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_params = {
            executor.submit(run_single_combination, params, n_pairs): params
            for params in combinations
        }
        progress = tqdm(as_completed(future_to_params),
                        total=n_combos, desc="Progress")
        for future in progress:
            summary, pairs = future.result()
            summary_results.append(summary)
            full_results.extend(pairs)

    total_time = time.time() - start_time
    print(f"\nSweep complete! Total time: {timedelta(seconds=int(total_time))}")
    return pd.DataFrame(summary_results), pd.DataFrame(full_results), total_time


def save_results(summary_df, full_df, timestamp, p, q, output_dir=None):
    """Save results to CSV files in the results directory."""
    if output_dir is None:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        network_label = f'p{p}q{int(q)}'
        output_dir = os.path.join(base, 'results', 'simulation', network_label)

    os.makedirs(output_dir, exist_ok=True)
    suffix = f"p{p}_q{int(q)}_beta_inf_p2_vs_thrshld"

    summary_file = os.path.join(
        output_dir, f'{timestamp}_sweep_summary_{suffix}.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSummary saved: {summary_file}")

    if not full_df.empty:
        full_file = os.path.join(
            output_dir, f'{timestamp}_sweep_full_{suffix}.csv')
        full_df.to_csv(full_file, index=False)
        print(f"Full results: {full_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Simulation sweep: p2 vs thrshld (beta=inf SI model)')
    parser.add_argument('--p', type=int, required=True,
                        help='Cycle-power parameter (1, 2, or 3)')
    parser.add_argument('--q', type=float, default=None,
                        help='Random component (default: 8-2p)')
    parser.add_argument('--n-pairs', type=int, default=500,
                        help='Network pairs per combination (default: 500)')
    parser.add_argument('--n-workers', type=int, default=None,
                        help='Parallel workers (default: SLURM or 8)')
    args = parser.parse_args()

    q = args.q if args.q is not None else float(8 - 2 * args.p)
    n_workers = args.n_workers or int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print("=" * 70)
    print(f"SIMULATION SWEEP: p2 vs thrshld (beta=infinity)")
    print(f"  Network: C_{{{args.p},{q:.0f}/n}} vs Random Regular (k=8)")
    print(f"  Timestamp: {timestamp}")
    print(f"  n_pairs={args.n_pairs}, n_workers={n_workers}")
    print("=" * 70)

    combinations = generate_parameter_combinations(args.p, q)
    print(f"  {len(combinations)} parameter combinations")

    summary_df, full_df, _ = run_parameter_sweep(
        combinations, args.n_pairs, n_workers)

    if not summary_df.empty:
        save_results(summary_df, full_df, timestamp, args.p, q)

    # Summary statistics
    n_A = (summary_df.get('case_A_count', pd.Series([0])) > 0).sum()
    n_B = (summary_df.get('case_B_count', pd.Series([0])) > 0).sum()
    n_C = (summary_df.get('case_C_count', pd.Series([0])) > 0).sum()
    print(f"\nCase distribution (combinations with >0 pairs):")
    print(f"  A={n_A}, B={n_B}, C={n_C}")
    print("\nSWEEP COMPLETE!")


if __name__ == "__main__":
    main()
