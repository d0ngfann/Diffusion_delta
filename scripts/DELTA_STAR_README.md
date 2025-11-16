# δ* (Delta-Star) Analysis - Finding the Critical Temporal Discount Threshold

## Overview

This analysis introduces a **temporal discounting perspective** to complex contagion simulations. Instead of only comparing final adoption rates between clustered and random networks, we examine how the **timing of adoptions** affects network performance.

### Key Concept: Temporal Discounting

The core idea is that adoptions at different time points might have different values:
- **Early adoptions** might be more valuable (e.g., getting a product to market quickly)
- **Later adoptions** might be less valuable (diminishing returns over time)

We control this preference using a **temporal discount factor δ** (delta):
- **δ ≈ 0**: Only immediate spread matters
- **δ ≈ 1**: All time steps valued equally (equivalent to total final adoption)

### The Critical Threshold: δ*

**δ* (delta-star)** is the specific value of δ where clustered and random networks perform **exactly equally**.

- **Below δ*** (δ < δ*): **Clustered networks win**
  - Social reinforcement drives fast early spread
  - When early adoptions are valued more, clustering helps

- **Above δ*** (δ > δ*): **Random networks win**
  - Broad reach enables sustained long-term diffusion
  - When all adoptions count equally, reach matters more

---

## Quick Start: Run the Demo

The demo script demonstrates the δ* concept using **mock data** and requires **no external dependencies**:

```bash
cd scripts
python3 delta_star_demo.py
```

This will:
1. Generate realistic mock time-series data for both network types
2. Calculate V (discounted cumulative activation) for different δ values
3. Find δ* where the networks perform equally
4. Display results in an ASCII visualization

**Expected output:**
- Results table showing V values for each δ
- ASCII plot comparing network performance
- Example calculation showing how V is computed

---

## Full Analysis: Real Simulations

### Step 1: Environment Setup

The full analysis requires a Python environment with:
- numpy
- pandas
- matplotlib
- networkx

**Option A: Use existing conda environment (on SLURM cluster)**
```bash
module load anaconda3/2022.01
source activate centola_riff_env
```

**Option B: Create new conda environment (local machine)**
```bash
conda create -n delta_star python=3.7 numpy pandas matplotlib networkx
conda activate delta_star
```

**Option C: Use pip**
```bash
pip install numpy pandas matplotlib networkx
```

### Step 2: Run Full Analysis

```bash
cd scripts
python delta_star_analysis.py
```

This will:
1. Run 20 trials on clustered networks (perc_rewire=0)
2. Run 20 trials on random networks (perc_rewire=1)
3. Calculate V for δ = [0.5, 0.7, 0.9, 0.95, 0.99]
4. Find δ* through interpolation
5. Generate publication-quality visualization
6. Save results to `output_data/delta_star_results.csv`
7. Save plot to `output_plots/delta_star_analysis.png`

**Default parameters:**
- k = 8 (network degree)
- n = 2000 (nodes)
- i = 2 (threshold)
- p1 = 0.3 (base rate)
- p2 = 0.6 (reinforced rate)
- beta = 1 (time of influence)
- trials = 20 per network type

### Step 3: Customize Parameters

Edit `delta_star_analysis.py` and modify the parameters in the `main()` function:

```python
def main():
    # Simulation parameters
    k = 8              # Network degree
    n = k * 250        # Number of nodes
    thrshld = 2        # Social reinforcement threshold
    p1 = 0.3           # Base adoption probability
    p2 = 0.6           # Reinforced adoption probability
    beta = 1           # Time of influence
    trials = 20        # Trials per network type
    seeds = 5          # Initial seeds

    # Delta values to test
    delta_values = [0.5, 0.7, 0.9, 0.95, 0.99]

    # ... rest of code
```

To test finer granularity near δ*, you can expand the delta range:

```python
delta_values = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99]
```

---

## Understanding the Output

### The V Metric

**V = c₁×δ + c₂×δ² + c₃×δ³ + ... + cₜ×δᵗ**

Where:
- **cₜ** = number of newly activated nodes at time step t
- **δ** = temporal discount factor

**Intuition:**
- Each new adoption is weighted by δᵗ
- Earlier adoptions (small t) receive higher weight when δ < 1
- As δ → 1, all adoptions weighted equally (V → total final adoption)

### Example Calculation

Suppose:
- New infections = [5, 10, 20, 15, 5, 0, 0, ...]
- δ = 0.9

Then:
```
V = 5×(0.9¹) + 10×(0.9²) + 20×(0.9³) + 15×(0.9⁴) + 5×(0.9⁵)
  = 5×0.9 + 10×0.81 + 20×0.729 + 15×0.6561 + 5×0.59049
  = 4.5 + 8.1 + 14.58 + 9.8415 + 2.95245
  = 39.97
```

### Interpreting δ*

If you find **δ* = 0.92**, this means:

1. **For immediate/short-term goals (δ < 0.92):**
   - Use **clustered networks**
   - Social reinforcement enables rapid early cascades
   - Examples: viral marketing campaigns, time-sensitive behaviors

2. **For long-term/sustained diffusion (δ > 0.92):**
   - Use **random networks**
   - Broader reach enables continued spread
   - Examples: knowledge diffusion, lasting behavior change

3. **The closer δ* is to 1:**
   - The more dominant random networks are overall
   - Clustered networks only win when early speed is highly valued

4. **The closer δ* is to 0:**
   - The more dominant clustered networks are
   - Random networks only win in the very long run

---

## Implementation Details

### Files Created

1. **`delta_star_analysis.py`**
   - Full implementation with real simulations
   - Requires numpy, pandas, matplotlib, networkx
   - Generates publication-quality figures
   - Saves detailed CSV results

2. **`delta_star_demo.py`**
   - Demonstration with mock data
   - No dependencies required (pure Python)
   - Educational/conceptual illustration
   - ASCII visualization

3. **`DELTA_STAR_README.md`** (this file)
   - Complete documentation
   - Usage instructions
   - Theoretical background

### Key Functions

**`cumulative_to_new_infections(cumulative_series)`**
- Converts cumulative adoption counts to per-timestep new infections
- Input: [5, 8, 12, 15, ...] (cumulative)
- Output: [5, 3, 4, 3, ...] (new per step)

**`calculate_V(new_infections, delta)`**
- Computes discounted cumulative activation
- Core metric for the analysis
- Formula: V = Σ cₜ × δᵗ

**`run_delta_simulations(...)`**
- Runs simulations for both network types
- Uses existing `run_simulation()` from helper_functions.py
- Sets `full_series=True` to capture time-series data

**`compute_V_across_deltas(results_df, delta_values)`**
- Calculates V for all trials and delta values
- Returns matrix: rows = trials, columns = deltas
- Used to compute means and standard errors

**`find_delta_star(delta_values, V_clustered, V_random)`**
- Finds intersection point via linear interpolation
- Returns exact δ* value
- Returns None if no intersection in tested range

**`visualize_delta_analysis(...)`**
- Creates publication-quality matplotlib figure
- Shows V vs. δ for both network types
- Marks δ* with vertical line
- Includes error bands (standard error)

---

## Extending the Analysis

### Testing Different Parameters

You can run the analysis across different conditions:

```python
# Test different thresholds
for thrshld in [2, 3, 4]:
    cf_c, cf_r = run_delta_simulations(thrshld=thrshld, ...)
    # ... analyze ...

# Test different p2/p1 ratios
for p2 in [0.4, 0.6, 0.8]:
    cf_c, cf_r = run_delta_simulations(p1=0.3, p2=p2, ...)
    # ... analyze ...
```

### Batch Processing

To systematically explore parameter space:

```python
import itertools

# Define parameter grid
k_values = [4, 6, 8]
thrshld_values = [2, 3, 4]
p2_values = [0.4, 0.6, 0.8]

results = []
for k, thrshld, p2 in itertools.product(k_values, thrshld_values, p2_values):
    print(f"Testing: k={k}, i={thrshld}, p2={p2}")

    cf_c, cf_r = run_delta_simulations(k=k, thrshld=thrshld, p2=p2)
    V_c_mean, _ = compute_V_across_deltas(cf_c, delta_values)
    V_r_mean, _ = compute_V_across_deltas(cf_r, delta_values)

    delta_star = find_delta_star(delta_values, V_c_mean, V_r_mean)

    results.append({
        'k': k,
        'thrshld': thrshld,
        'p2': p2,
        'delta_star': delta_star
    })

# Save batch results
pd.DataFrame(results).to_csv('delta_star_batch_results.csv')
```

### Alternative Delta Definitions

The current implementation uses δᵗ⁺¹ (seeds at t=0 count as "time 1"). You could also try:

**Option 1: δᵗ (seeds at t=0 get full weight)**
```python
V += new_infections[t] * (delta ** t)
```

**Option 2: Exponential discounting with rate λ**
```python
V += new_infections[t] * np.exp(-lambda_rate * t)
```

**Option 3: Hyperbolic discounting**
```python
V += new_infections[t] / (1 + k * t)
```

---

## Theoretical Connections

### Relation to Existing Metrics

**When δ = 1:**
- V = total final adoption (max_spread_norm)
- Equivalent to existing analysis

**When δ → 0:**
- V ≈ c₁×δ (only first timestep matters)
- Measures immediate spread velocity

**Intermediate δ:**
- Balances early speed vs. final reach
- Captures the trade-off between reinforcement and reach

### Why This Matters

The δ* analysis reveals:

1. **Context-dependent optimal networks**
   - No single network is universally "best"
   - Depends on temporal preferences

2. **Mechanism insights**
   - Clustered networks: fast start, early saturation
   - Random networks: slower start, sustained growth

3. **Practical implications**
   - Design interventions based on time horizons
   - Marketing (low δ) vs. education (high δ)

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'numpy'"
→ Install required packages (see Environment Setup)

### "No intersection found in the given delta range"
→ Try expanding the delta_values range or adjusting simulation parameters

### "Warning: V values are very small/large"
→ This is normal; V scales with network size and delta values

### Simulations taking too long
→ Reduce `trials` from 20 to 10, or reduce network size `n`

### Results seem noisy
→ Increase `trials` from 20 to 50 or 100 for more stable estimates

---

## Citation

If you use this analysis in your research, please cite:

```
[Your paper citation here once published]

Delta-star (δ*) analysis: A temporal discounting framework for
evaluating network performance in complex contagion simulations.
```

---

## Contact

For questions or issues:
- Check existing simulation documentation in `CLAUDE.md`
- Review `helper_functions.py` for simulation details
- Open an issue on the repository

---

## Version History

- **v1.0** (2025): Initial implementation
  - Core δ* analysis framework
  - Demo and full analysis scripts
  - Visualization tools
