"""
Theoretical Sweep: p2 vs threshold for beta=infinity (SI model)

Computes closed-form delta* for C_{p,q/n} vs Random Regular across a grid
of (p2, threshold) values using the MGF and integral approaches.

Much faster than simulation (no network generation or stochastic trials).

Usage:
    cd scripts/beta_infinite_delta_star
    python sweep_theory/sweep_p2_vs_thrshld.py --p 1
    python sweep_theory/sweep_p2_vs_thrshld.py --p 2
    python sweep_theory/sweep_p2_vs_thrshld.py --p 3
"""

import os
import sys
import time
import argparse
from datetime import datetime, timedelta
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.theory import compute_sweep_p2_vs_thrshld


def save_results(df, timestamp, p, q, p1, output_dir=None):
    """Save theoretical sweep results."""
    if output_dir is None:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        network_label = f'p{p}q{int(q)}'
        output_dir = os.path.join(base, 'results', 'theory', network_label)

    os.makedirs(output_dir, exist_ok=True)
    suffix = f"p{p}_q{int(q)}_beta_inf_p2_vs_thrshld"
    filename = f'{timestamp}_theory_sweep_{suffix}_p1_{p1}.csv'
    filepath = os.path.join(output_dir, filename)
    df.to_csv(filepath, index=False)
    print(f"Saved: {filepath}")
    return filepath


def main():
    parser = argparse.ArgumentParser(
        description='Theoretical sweep: p2 vs thrshld (beta=inf)')
    parser.add_argument('--p', type=int, required=True,
                        help='Cycle-power parameter (1, 2, or 3)')
    parser.add_argument('--q', type=float, default=None,
                        help='Random component (default: 8-2p)')
    parser.add_argument('--p1', type=float, default=0.1,
                        help='Base adoption probability')
    parser.add_argument('--N', type=int, default=2000,
                        help='Network size')
    args = parser.parse_args()

    q = args.q if args.q is not None else float(8 - 2 * args.p)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print("=" * 70)
    print(f"THEORETICAL SWEEP: p2 vs thrshld (beta=infinity)")
    print(f"  Network: C_{{{args.p},{q:.0f}/n}} vs Random Regular (k=8)")
    print(f"  p1={args.p1}, N={args.N}")
    print(f"  Timestamp: {timestamp}")
    print("=" * 70)

    start_time = time.time()

    p2_values = np.arange(0.50, 1.01, 0.01)
    thrshld_values = np.arange(2.0, 5.1, 0.1)

    total = len(p2_values) * len(thrshld_values)
    print(f"  {total} parameter combinations")
    print("  Computing...", end='', flush=True)

    results = compute_sweep_p2_vs_thrshld(
        args.p, q, p1=args.p1, N=args.N,
        p2_values=p2_values, thrshld_values=thrshld_values,
    )

    elapsed = time.time() - start_time
    print(f" done ({elapsed:.1f}s)")

    df = pd.DataFrame(results)

    # Format output columns to match simulation sweep format
    df_out = pd.DataFrame({
        'combination_id': range(1, len(df) + 1),
        'n': df['N'],
        'k': 8,
        'p': df['p'],
        'q': df['q'],
        'thrshld': df['thrshld'],
        'p1': df['p1'],
        'p2': df['p2'],
        'beta': 'inf',
        'seeds': df['seeds'],
        'case_A_count': (df['case'] == 'A').astype(int),
        'case_B_count': (df['case'] == 'B').astype(int),
        'case_C_count': (df['case'] == 'C').astype(int),
        'case_Invalid_count': (df['case'] == 'Invalid').astype(int),
        'delta_star_mean': df['delta_star'],
        'delta_star_mgf': df['delta_star_mgf'],
        'delta_star_integral': df['delta_star_integral'],
        'g_R_half': df['g_R_half'],
        'g_C_half': df['g_C_half'],
        'r_R': df['r_R'],
        'r_C': df['r_C'],
        'r_H': df['r_H'],
        'status': 'success',
    })

    save_results(df_out, timestamp, args.p, q, args.p1)

    # Summary
    n_A = (df['case'] == 'A').sum()
    n_B = (df['case'] == 'B').sum()
    n_C = (df['case'] == 'C').sum()
    print(f"\nCase distribution: A={n_A}, B={n_B}, C={n_C}")
    if n_C > 0:
        ds_vals = df.loc[df['case'] == 'C', 'delta_star'].dropna()
        print(f"  delta* range: [{ds_vals.min():.4f}, {ds_vals.max():.4f}]")
        print(f"  delta* mean:  {ds_vals.mean():.4f}")

    print(f"\nTotal time: {timedelta(seconds=int(elapsed))}")
    print("DONE!")


if __name__ == "__main__":
    main()
