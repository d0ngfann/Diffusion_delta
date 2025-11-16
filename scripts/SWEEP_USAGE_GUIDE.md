# Delta-Star Parameter Sweep - Usage Guide

## Overview

The `delta_star_sweep.py` script performs a systematic parameter sweep to find delta-star (δ*) values across multiple combinations of network and diffusion parameters.

## Quick Start

```bash
cd scripts
python3 delta_star_sweep.py
```

The script will:
1. Generate all valid parameter combinations (513 total)
2. Run simulations for each combination
3. Save results to CSV
4. Generate summary statistics
5. Create visualization plots

**Estimated runtime**: ~43 minutes (assuming 5 seconds per combination)

## Parameter Grid

The sweep tests combinations of:

- **k (network degree)**: [4, 8, 12]
- **thrshld (threshold)**: [2, 3, 4]
- **p1 (base adoption)**: [0.1, 0.3, 0.5, 0.7, 0.9]
- **p2 (reinforced adoption)**: [0.3, 0.5, 0.7, 0.9, 1.0]
- **beta (time of influence)**: [1, 3, 5]

**Constraints**:
- Only tests combinations where `p1 ≤ p2` (physical constraint)
- `seeds` always equals `thrshld` (recommended practice)

**Total combinations**: 513 valid combinations

## Output Files

### 1. Results CSV (`output_data/sweep_results.csv`)

Contains one row per parameter combination with columns:
- `combination_id`: Sequential ID number
- `k`, `thrshld`, `p1`, `p2`, `beta`, `seeds`: Parameter values
- `primary_delta_star`: The first meaningful δ* found (or None)
- `num_intersections`: Number of meaningful intersections found
- `status`: "found", "not_found", "multiple", or "error"
- `intersection_1` through `intersection_5`: Additional intersection values

### 2. Summary Statistics (`output_data/sweep_summary.txt`)

Text file containing:
- Overall statistics (success rate, etc.)
- Delta-star statistics when found (mean, median, std, min, max)
- Statistics grouped by each parameter
- Success rate by parameter

### 3. Visualization Plots (`output_plots/sweep_plots/`)

Four plots are generated:

**a) Heatmap: k vs thrshld** (`heatmap_k_vs_thrshld.png`)
- Y-axis: Network degree (k)
- X-axis: Social reinforcement threshold (i)
- Color: Average delta-star value
- Shows how network structure and threshold affect δ*

**b) Scatter: p1 vs p2** (`scatter_p1_vs_p2.png`)
- X-axis: Base adoption probability (p₁)
- Y-axis: Reinforced adoption probability (p₂)
- Color: Delta-star value
- Size: Beta value (time of influence)
- Shows how adoption rates affect δ*

**c) Box Plot: beta distribution** (`boxplot_beta.png`)
- X-axis: Beta values (1, 3, 5)
- Y-axis: Delta-star distribution
- Shows how time of influence affects δ* variability

**d) Pie Chart: status distribution** (`pie_chart_status.png`)
- Shows proportion of combinations with different outcomes:
  - Green: Single intersection found
  - Orange: Multiple intersections
  - Red: No intersection found
  - Gray: Errors

## Console Output

During execution, you'll see progress like:

```
Starting parameter sweep...
Total combinations to test: 513
Estimated time: 0:42:45

Progress: [===>    ] 25/513 (4.9%)
✓ [25/513] k=4, i=2, p1=0.3, p2=0.5, β=1 → δ*=0.4523 | Elapsed: 0:02:05 | Remaining: 0:40:40
○ [26/513] k=4, i=2, p1=0.3, p2=0.5, β=3 → δ*=None | Elapsed: 0:02:10 | Remaining: 0:40:35
```

**Status symbols**:
- `✓`: Delta-star found successfully
- `○`: No intersection found
- `✗`: Error occurred

## Customizing the Sweep

### Modifying Parameter Grids

Edit `delta_star_sweep.py`, function `generate_parameter_combinations()`:

```python
# Change these lines (around line 57-61)
k_values = [4, 8, 12]              # Add/remove values
thrshld_values = [2, 3, 4]         # Add/remove values
p1_values = [0.1, 0.3, 0.5, 0.7, 0.9]  # Add/remove values
p2_values = [0.3, 0.5, 0.7, 0.9, 1.0]   # Add/remove values
beta_values = [1, 3, 5]            # Add/remove values
```

### Running a Subset for Testing

If you want to test with just a few combinations first:

```python
# Edit the main() function around line 429
combinations = generate_parameter_combinations()
combinations = combinations[:10]  # Test with first 10 only
```

### Changing Number of Trials

Edit `run_single_combination()` function (line 124):

```python
delta_star, results_summary, all_intersections = run_delta_star_analysis(
    ...
    trials=50,  # Change this number
    ...
)
```

**Note**: Fewer trials = faster but less reliable. More trials = slower but more robust.

## Advanced Usage

### Running Programmatically

```python
import delta_star_sweep as sweep

# Generate combinations
combos = sweep.generate_parameter_combinations()

# Run just specific combinations
subset = combos[0:5]  # First 5 combinations
results_df = sweep.run_parameter_sweep(subset)

# Save and visualize
sweep.save_results(results_df)
sweep.generate_summary_statistics(results_df)
sweep.create_visualizations(results_df)
```

### Analyzing Results

```python
import pandas as pd

# Load results
df = pd.read_csv('../output_data/sweep_results.csv')

# Filter to successful cases
df_found = df[df['status'] == 'found']

# Find parameter combinations with delta star < 0.5
early_regime = df_found[df_found['primary_delta_star'] < 0.5]

# Group by threshold and find average delta star
avg_by_threshold = df_found.groupby('thrshld')['primary_delta_star'].mean()
print(avg_by_threshold)

# Find combinations with multiple intersections
complex_cases = df[df['status'] == 'multiple']
print(f"Found {len(complex_cases)} cases with non-monotonic behavior")
```

## Troubleshooting

### "ModuleNotFoundError: No module named 'tqdm'"

The script will work without tqdm, but you won't get a fancy progress bar. To install:
```bash
pip install tqdm
```

### "ModuleNotFoundError: No module named 'seaborn'"

The script will create simpler visualizations. To install for better plots:
```bash
pip install seaborn
```

### Script Running Slowly

Each combination takes ~5 seconds with default settings (50 trials). To speed up:
- Reduce number of trials (less reliable)
- Reduce number of delta values tested
- Run a subset of combinations

### Out of Memory

If running many combinations:
- Process in batches (run subsets separately)
- Reduce number of trials
- Reduce network size (though n=k*250 is already small)

## Expected Results

Based on the parameter ranges:
- **High success rate** expected for combinations where p1 << p2 (strong reinforcement)
- **Low success rate** expected when p1 ≈ p2 (weak reinforcement)
- **Multiple intersections** more likely with intermediate beta values
- **No intersection** possible when one network dominates across all δ values

## Files Created

```
output_data/
├── sweep_results.csv          # Full results table
└── sweep_summary.txt           # Summary statistics

output_plots/sweep_plots/
├── heatmap_k_vs_thrshld.png   # Heatmap
├── scatter_p1_vs_p2.png       # Scatter plot
├── boxplot_beta.png           # Box plot
└── pie_chart_status.png       # Pie chart
```

## Performance Tips

1. **Run overnight**: 513 combinations × 5 sec ≈ 43 minutes
2. **Test first**: Run with `combinations[:10]` to verify everything works
3. **Use screen/tmux**: For long runs on servers
4. **Monitor progress**: Check the console output periodically
5. **Save partial results**: The script saves after completing all combinations

## Next Steps After Running

1. **Examine summary.txt** for overall patterns
2. **Load CSV in Excel/R** for custom analysis
3. **Look at visualizations** to identify trends
4. **Investigate interesting cases**:
   - Multiple intersections (non-monotonic)
   - Very low/high δ* values
   - Parameter combinations where no intersection found

## Citation

If using this analysis in publications, cite:
```
Delta-Star Parameter Sweep Analysis
[Your Lab/Institution]
Generated using Claude Code
```

---

Last updated: 2025-11-11
