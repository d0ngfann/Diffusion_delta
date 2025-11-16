# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains replication code for a complex contagion simulation study. The codebase implements network-based diffusion models to study how behaviors spread through social networks under different conditions (network topology, adoption thresholds, time of influence, etc.).

## Core Architecture

### Simulation Pipeline
The simulation workflow follows a three-stage process:

1. **Configuration Generation** (`scripts/get_batch_files.r`): Generates batch configuration files (.conf) with parameter combinations for simulations
2. **Simulation Execution** (`scripts/simulation.py`): Runs the core diffusion model with specified parameters
3. **Analysis & Visualization** (`scripts/figures_code/*.r`): Processes simulation output and generates figures

### Key Python Components

**`scripts/helper_functions.py`**: Core simulation engine containing:
- `sirsic_wmem()`: Main diffusion model implementing SIR-like dynamics with memory/time of influence (beta parameter)
- `run_simulation()`: Wrapper that executes multiple trials and post-processes results
- `main()`: Top-level function that loops over parameter ranges and saves outputs
- Network construction functions (e.g., `moore_lattice()`)
- Seeding strategies: `seed_strat_random()`, `seed_strat_adj()`, `seed_strat_one_nbrs()`

**`scripts/simulation.py`**: Entry point script that:
- Parses command-line arguments for simulation parameters
- Constructs networks (Watts-Strogatz "WS", Moore lattice "MR", Hexagonal "HX")
- Calls `main()` from helper_functions to run simulations
- Saves trial-level results to CSV files in `output_data/cf/` directory

### Key R Components

**`scripts/figures_code/helper_functions.r`**: Plotting and analysis utilities including:
- `get_simulation_conf_file()`: Generates parameter grids for batch runs
- `get_long_file()`, `save_long_file()`: Data aggregation across simulation files
- Plotting functions for heatmaps and difference plots
- `get_p2_T()`: Theoretical threshold calculations

## Common Commands

### Running Simulations

**Generate batch configuration files:**
```r
cd scripts
Rscript get_batch_files.r
```
This creates .conf files (e.g., `batch_main_perc0.conf`, `batch_main_perc1.conf`) with simulation commands.

**Run a single simulation manually:**
```bash
cd scripts
python simulation.py --trials 100 --G_name 'WS' --n 2000 --k 8 --b 1 --thrshld 2 --p1 0.1 --perc 0 --seed_strat 'one_nbrs' --sig_thresh_sd 0 --rand 'ES'
```

**Run simulations via SLURM (on HPC):**
```bash
sbatch --array=0-[ARRAY_LENGTH] run_simulation.sh
```
Edit `run_simulation.sh` to change which .conf file is used (line 18).

### Generating Figures

Navigate to `scripts/figures_code/` and run the appropriate R script:
```bash
cd scripts/figures_code
Rscript main_fig1.r
Rscript main_fig2_si_EF.r
Rscript main_fig3_si_C.r
# etc.
```

Plots save to `output_plots/` directory.

## Important Parameters

- **G_name**: Network type - "WS" (Watts-Strogatz ring), "MR" (Moore lattice), "HX" (Hexagonal)
- **k**: Network degree (all nodes have same degree)
- **n**: Number of nodes (typically k*250)
- **b (beta)**: Time of influence - how long adopters remain influential (0 = infinite)
- **thrshld (i)**: Social reinforcement threshold - exposures needed to adopt at p2 rate
- **p1**: Base adoption probability (below threshold)
- **p2**: Reinforced adoption probability (at/above threshold)
- **perc**: Proportion of edges to rewire (0 = clustered, 1 = random network)
- **seed_strat**: Initial seeding strategy ("random", "adj", "one_nbrs")
- **sig_thresh_sd**: Heterogeneity in individual adoption (0 = homogeneous)
- **rand**: Randomization method (typically "ES" for edge swap)

## File Organization

- **`scripts/`**: Main simulation code (Python) and configuration generation (R)
- **`scripts/figures_code/`**: Analysis and plotting scripts (R)
- **`output_data/cf/`**: Simulation trial-level outputs (CSV files)
- **`output_plots/`**: Generated figures
- **`output_job/`**: SLURM job logs (when running on HPC)

## Key Concepts

**Adoption Trajectory**: The probability function `prob_list` defines adoption rates based on number of distinct influential neighbors. Typically: [0, p1, p1, ..., p2, p2, ...] where threshold determines the transition point.

**Memory/Time of Influence (beta)**: Controls how long adopted nodes remain influential. Beta=1 means only newly adopted nodes in the previous timestep are influential. Beta=0 means all adopters remain influential throughout.

**Network Rewiring**: The `perc` parameter controls interpolation between clustered lattices (perc=0) and random networks (perc=1) using double edge swaps while preserving degree distribution.

## Output Data Structure

Trial-level CSV files contain columns:
- `max_spread`, `max_spread_norm`: Final adoption counts and proportions
- `time_to_spread`, `time_to_X_spread`: Time to saturation at various thresholds
- Network and model parameters (k, n_nodes, thrshld, p1, p2, etc.)
- `seed`: Trial number/random seed
