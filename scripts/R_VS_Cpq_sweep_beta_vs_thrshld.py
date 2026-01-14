# python3 R_VS_Cpq_sweep_beta_vs_thrshld.py --p 2 --q 4 --p2 1.0 --n-pairs 100
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
from R_VS_Cpq import analyze_pair

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    class tqdm:
        def __init__(self, iterable, desc="", unit=""):
            self.iterable = iterable
            self.desc = desc
        def __iter__(self):
            return iter(self.iterable)
        @staticmethod
        def write(msg):
            print(msg)

def generate_parameter_combinations(p, q, p2):
    """
    Generate all valid parameter combinations for the beta vs thrshld sweep.

    Parameters:
        p: int - Cycle-power parameter
        q: float - Random component parameter
        p2: float - Reinforced adoption probability (fixed for this sweep)
    """
    n_values = [2000]
    k_values = [8]
    p1_values = [0.01]
    p2_values = [p2]
    beta_values = np.arange(1.0, 10.2, 0.2)
    thrshld_values = np.arange(2.0, 5.1, 0.1)
    all_combos = list(itertools.product(
        n_values, k_values, thrshld_values, p1_values, p2_values, beta_values
    ))
    valid_combos = []
    combination_id = 1
    for n, k, thrshld, p1, p2, beta in all_combos:
        if p1 <= p2:
            valid_combos.append({
                'combination_id': combination_id,
                'n': int(n),
                'k': int(k),
                'p': int(p),
                'q': float(q),
                'thrshld': float(thrshld),
                'p1': float(p1),
                'p2': float(p2),
                'beta': float(beta),
                'seeds': int(math.ceil(thrshld))
            })
            combination_id += 1
    return valid_combos

def run_single_combination(params, n_pairs=100, delta_values=None, verbose=False):
    if delta_values is None:
        delta_values = np.linspace(0, 1, 201)
    try:
        pair_results = []
        for pair_id in range(n_pairs):
            result = analyze_pair(
                pair_id=pair_id, n=params['n'], k=params['k'], p=params['p'],
                q=params['q'], thrshld=params['thrshld'], p1=params['p1'],
                p2=params['p2'], beta=params['beta'], seeds=params['seeds'],
                delta_values=delta_values
            )
            pair_results.append({'combination_id': params['combination_id'], 'pair_id': pair_id, **params, **result})
        
        case_counts = pd.Series([r['case'] for r in pair_results]).value_counts()
        delta_stars = [r['delta_star'] for r in pair_results if r.get('case') == 'C' and r.get('delta_star') is not None]
        
        summary = {
            **params,
            'case_A_count': case_counts.get('A', 0), 'case_B_count': case_counts.get('B', 0),
            'case_C_count': case_counts.get('C', 0), 'case_Invalid_count': case_counts.get('Invalid', 0),
            'delta_star_mean': np.mean(delta_stars) if delta_stars else None,
            'delta_star_std': np.std(delta_stars, ddof=1) if len(delta_stars) > 1 else None,
            'status': 'success'
        }
        return summary, pair_results
    except Exception as e:
        return {**params, 'status': 'error', 'error_message': str(e)}, []

def run_parameter_sweep(combinations, n_pairs=100, n_workers=8):
    n_combos = len(combinations)
    summary_results, full_results = [], []
    print(f"\nStarting sweep: {n_combos} combinations, {n_pairs} pairs each...")
    start_time = time.time()
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_params = {executor.submit(run_single_combination, params, n_pairs): params for params in combinations}
        
        progress = tqdm(as_completed(future_to_params), total=n_combos, desc="Progress")
        for future in progress:
            summary, pairs = future.result()
            summary_results.append(summary)
            full_results.extend(pairs)
    
    total_time = time.time() - start_time
    print(f"\nSweep complete! Total time: {timedelta(seconds=int(total_time))}")
    return pd.DataFrame(summary_results), pd.DataFrame(full_results), total_time

def save_results(summary_df, full_df, timestamp, p, q, p2, output_dir='../output_data'):
    os.makedirs(output_dir, exist_ok=True)
    suffix = f"p{p}_q{q}_beta_vs_thrshld_p2_{p2}"
    summary_file = os.path.join(output_dir, f'{timestamp}_sweep_summary_{suffix}.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSummary saved to: {summary_file}")

    if not full_df.empty:
        full_file = os.path.join(output_dir, f'{timestamp}_sweep_full_{suffix}.csv')
        full_df.to_csv(full_file, index=False)
        print(f"Full results saved to: {full_file}")

def main():
    parser = argparse.ArgumentParser(description='Sweep beta vs thrshld for R vs C_p,q/n analysis')
    parser.add_argument('--p', type=int, required=True, help='Cycle-power parameter')
    parser.add_argument('--q', type=float, required=True, help='Random component')
    parser.add_argument('--p2', type=float, required=True, help='Reinforced adoption probability (fixed for this sweep)')
    parser.add_argument('--n-pairs', type=int, default=500, help='Pairs per combination')
    parser.add_argument('--n-workers', type=int, default=None, help='Num parallel workers')
    args = parser.parse_args()

    n_workers = args.n_workers or int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print("="*70)
    print("DELTA-STAR SWEEP: beta vs thrshld")
    print(f"Timestamp: {timestamp}")
    print(f"Fixed p2: {args.p2}")

    combinations = generate_parameter_combinations(args.p, args.q, args.p2)
    summary_df, full_df, _ = run_parameter_sweep(combinations, args.n_pairs, n_workers)

    if not summary_df.empty:
        save_results(summary_df, full_df, timestamp, args.p, args.q, args.p2)
    print("\nSWEEP COMPLETE!")

if __name__ == "__main__":
    main()

