# Complex Contagion Network Comparison: Random Regular vs. Cycle-Power-Union-Random

This repository contains simulation code for comparing diffusion dynamics between Random Regular graphs and Cycle-Power-Union-Random (C_{p,q/n}) hybrid networks using delta-star analysis.

## Overview

This research investigates how network structure (clustering vs. shortcuts) affects complex contagion diffusion under different parameter conditions including:
- Adoption probabilities (p2)
- Social reinforcement thresholds
- Temporal memory effects (beta)

**Key Innovation**: Uses linear path seeding (`seed_strat_path`) to create fair comparisons between different network topologies.

## Related Publication

This work builds upon:
> Guilbeault, D., Becker, J., & Centola, D. (2018). *Complex contagions: A decade in review*. Complex Spreading Phenomena in Social Systems, 3-25.

See: `Diffusion of complex contagions is shaped by a trade-off between reach and reinforcement.pdf`

## Repository Structure

```
.
├── scripts/                          # Core simulation code
│   ├── DH_helper_functions.py        # Simulation engine & seeding strategies
│   ├── R_VS_Cpq.py                   # Network generation & analysis functions
│   ├── R_VS_Cpq_sweep_p2_vs_thrshld.py      # p2 vs threshold sweep
│   ├── R_VS_Cpq_sweep_beta_vs_thrshld.py    # beta vs threshold sweep
│   ├── R_VS_Cpq_sweep_p2_vs_beta.py         # p2 vs beta sweep
│   ├── plot_deltastar_heatmaps.py           # Visualization (p1=0.1)
│   └── plot_deltastar_heatmaps_p1_0.01.py   # Visualization (p1=0.01)
├── sbatch/                           # SLURM job scripts (40 experiments)
├── output_data/                      # Simulation results
│   ├── 20260101/                     # January 2026 results
│   └── 20260108_p1_0.01/            # p1=0.01 parameter exploration
├── output_plots/                     # Generated visualizations
├── report.tex                        # Experimental design & results
├── requirements.txt                  # Python dependencies
└── README.md                         # This file
```

## Installation

### Prerequisites
- Python 3.8+
- HPC cluster with SLURM (for large-scale experiments)

### Setup

```bash
# Clone the repository
git clone https://github.com/d0ngfann/diffusion_network_RvsCpq.git
cd diffusion_network_RvsCpq

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Network Types

### 1. Random Regular Graph (Baseline)
- All nodes have exactly degree k=8
- Randomly connected, no inherent clustering

### 2. Cycle-Power-Union-Random (C_{p,q/n})
Four configurations tested:
- **p4q0**: Pure lattice (8 lattice edges, ~0 random)
- **p3q2**: Lattice-heavy (6 lattice, ~2 random)
- **p2q4**: Balanced (4 lattice, ~4 random)
- **p1q6**: Random-heavy (2 lattice, ~6 random)

Construction: C_{p,q/n} = C_p ∪ G_{n,q/n}
- C_p: Cycle-power-p graph (2p lattice edges per node)
- G_{n,q/n}: Erdős-Rényi random graph (edge probability q/n)

## Usage

### Single Experiment

Run a single parameter sweep:

```bash
cd scripts
python R_VS_Cpq_sweep_p2_vs_thrshld.py --p 2 --q 4 --beta 5.0
```

Optional arguments:
- `--n-pairs 500`: Number of network pairs (default: 500)
- `--n-workers 8`: Parallel workers (default: from SLURM or 8)

### Complete Experimental Suite

Run all 40 experiments (p2 vs threshold, beta vs threshold, p2 vs beta):

```bash
cd sbatch
bash submit_all_jobs.sh
```

Individual experiments can be submitted:
```bash
sbatch exp1_p2thrshld_p2q4_beta5.sbatch
sbatch exp2_betathrshld_p1q6_p2_0.8.sbatch
sbatch exp3_p2beta_p4q0_thrshld3.sbatch
```

### Visualization

Generate heatmaps from results:

```bash
cd scripts
python plot_deltastar_heatmaps.py        # For p1=0.1 results
python plot_deltastar_heatmaps_p1_0.01.py  # For p1=0.01 results
```

## Key Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| **n** | Network size | 2000 |
| **k** | Node degree | 8 |
| **p1** | Base adoption probability | 0.1 or 0.01 |
| **p2** | Reinforced adoption probability | 0.5 - 1.0 |
| **threshold** | Exposures needed for reinforcement | 2.0 - 5.0 |
| **beta** | Time of influence (memory) | 1.0 - 10.0 |
| **delta** | Temporal discount factor | 0.0 - 1.0 |

### Fractional Parameters

Both `threshold` and `beta` support fractional values:
- **threshold = 2.3**: Each node gets threshold=3 with probability 0.3, else threshold=2
- **beta = 5.7**: Each adopter stays influential for 6 timesteps with probability 0.7, else 5

## Delta-Star Metric

Delta-star (δ*) is the critical temporal discount factor where both networks achieve equal performance:

V(δ) = Σ c_t · δ^t

Where:
- c_t: New adopters at timestep t
- δ*: Value where V_Cpq(δ*) = V_Random(δ*)

**Interpretation**:
- High δ* (e.g., 0.7): C_{p,q/n} performs better for most time preferences
- Low δ* (e.g., 0.3): Random performs better for most time preferences
- No δ*: One network always dominates (Case A or B)

## Output Files

Results are saved to `output_data/`:
- `*_sweep_summary_*.csv`: Summary statistics per parameter combination
- `*_sweep_full_*.csv`: Individual network pair results

Case classification:
- **Case A**: C_{p,q/n} dominates (all δ values)
- **Case B**: Random dominates (all δ values)
- **Case C**: Intersection exists (δ* found)
- **Case Invalid**: No valid comparison region

## Citation

If you use this code, please cite:

```bibtex
@article{guilbeault2018complex,
  title={Complex contagions: A decade in review},
  author={Guilbeault, Douglas and Becker, Joshua and Centola, Damon},
  journal={Complex spreading phenomena in social systems},
  pages={3--25},
  year={2018},
  publisher={Springer}
}
```
