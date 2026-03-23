"""
Plotting utilities for beta=infinity delta-star analysis.

Generates heatmaps of delta* across (p2, threshold) parameter space,
with optional theoretical contour overlay.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')


NETWORK_LABELS = {
    (1, 6): r'$C_{1,6/n}$ (p1q6)',
    (2, 4): r'$C_{2,4/n}$ (p2q4)',
    (3, 2): r'$C_{3,2/n}$ (p3q2)',
}


def load_sweep_results(filepath):
    """Load sweep summary CSV."""
    return pd.read_csv(filepath)


def create_heatmap(df, x_col, y_col, value_col='delta_star',
                   title=None, ax=None, vmin=0, vmax=1, cmap='RdYlBu_r'):
    """
    Create a heatmap from sweep results.

    Parameters:
        df: DataFrame with columns x_col, y_col, value_col
        x_col: column for x-axis (e.g., 'p2')
        y_col: column for y-axis (e.g., 'thrshld')
        value_col: column for color values
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    pivot = df.pivot_table(index=y_col, columns=x_col, values=value_col)

    im = ax.pcolormesh(
        pivot.columns, pivot.index, pivot.values,
        cmap=cmap, vmin=vmin, vmax=vmax, shading='auto'
    )

    ax.set_xlabel(x_col, fontsize=12)
    ax.set_ylabel(y_col, fontsize=12)
    if title:
        ax.set_title(title, fontsize=14)

    return im


def create_case_heatmap(df, x_col, y_col, title=None, ax=None):
    """
    Create a heatmap showing case classification (A/B/C).

    Uses simulation data with case_A_count, case_B_count, case_C_count columns.
    Color encodes the dominant case.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Compute dominant case fraction
    total = (df.get('case_A_count', 0) + df.get('case_B_count', 0) +
             df.get('case_C_count', 0))
    total = total.replace(0, 1)  # avoid division by zero

    # Encode: A=-1, B=1, C=0 (weighted average)
    case_score = (
        -1 * df.get('case_A_count', 0) / total +
        1 * df.get('case_B_count', 0) / total
    )

    df_plot = df.copy()
    df_plot['case_score'] = case_score

    pivot = df_plot.pivot_table(index=y_col, columns=x_col, values='case_score')

    im = ax.pcolormesh(
        pivot.columns, pivot.index, pivot.values,
        cmap='RdBu', vmin=-1, vmax=1, shading='auto'
    )

    ax.set_xlabel(x_col, fontsize=12)
    ax.set_ylabel(y_col, fontsize=12)
    if title:
        ax.set_title(title, fontsize=14)

    return im


def overlay_theory_contours(ax, theory_df, x_col, y_col,
                            contour_levels=None, color='black'):
    """
    Overlay theoretical delta* contour lines on a heatmap.

    Parameters:
        ax: matplotlib axes
        theory_df: DataFrame with theoretical results
        x_col, y_col: column names for axes
        contour_levels: list of delta* values for contour lines
    """
    case_C = theory_df[theory_df['case'] == 'C'].copy()
    if case_C.empty:
        return

    ds_col = 'delta_star'
    if ds_col not in case_C.columns:
        ds_col = 'delta_star_mgf'

    pivot = case_C.pivot_table(index=y_col, columns=x_col, values=ds_col)

    if pivot.empty:
        return

    X, Y = np.meshgrid(pivot.columns, pivot.index)
    Z = pivot.values

    if contour_levels is None:
        valid = Z[~np.isnan(Z)]
        if len(valid) == 0:
            return
        contour_levels = np.linspace(valid.min(), valid.max(), 6)[1:-1]

    cs = ax.contour(
        X, Y, Z, levels=contour_levels, colors=color,
        linewidths=1.5, linestyles='--'
    )

    ax.clabel(cs, inline=True, fontsize=8, fmt='%.2f',
              colors=color)


def plot_comparison(sim_df, theory_df, p, q, x_col='p2', y_col='thrshld',
                    output_path=None):
    """
    Create side-by-side comparison: simulation heatmap + theory overlay.
    """
    label = NETWORK_LABELS.get((p, q), f'C_{{{p},{q}/n}}')

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Left: simulation case distribution
    if sim_df is not None and not sim_df.empty:
        im1 = create_case_heatmap(
            sim_df, x_col, y_col,
            title=f'{label} - Case Distribution (Sim)',
            ax=axes[0]
        )
        plt.colorbar(im1, ax=axes[0], label='B(-1) / C(0) / A(+1)')

    # Right: delta* heatmap with theory overlay
    ds_col = 'delta_star_mean' if 'delta_star_mean' in (
        sim_df.columns if sim_df is not None else []) else 'delta_star'

    source_df = sim_df if sim_df is not None else theory_df
    im2 = create_heatmap(
        source_df, x_col, y_col, value_col=ds_col,
        title=f'{label} - delta* (beta=inf)',
        ax=axes[1]
    )
    plt.colorbar(im2, ax=axes[1], label='delta*')

    if theory_df is not None and not theory_df.empty:
        overlay_theory_contours(axes[1], theory_df, x_col, y_col)

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close()
    return fig


def plot_theory_only(theory_df, p, q, x_col='p2', y_col='thrshld',
                     output_path=None):
    """
    Plot theoretical delta* heatmap only (no simulation data needed).
    """
    label = NETWORK_LABELS.get((p, q), f'C_{{{p},{q}/n}}')

    fig, ax = plt.subplots(figsize=(10, 8))

    im = create_heatmap(
        theory_df, x_col, y_col, value_col='delta_star',
        title=f'{label} - Theoretical delta* (beta=inf)',
        ax=ax
    )
    plt.colorbar(im, ax=ax, label='delta*')
    overlay_theory_contours(ax, theory_df, x_col, y_col)

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_path}")

    plt.close()
    return fig
