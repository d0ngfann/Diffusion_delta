# Beta=Infinity Delta-Star Analysis

Closed-form and simulation-based delta* computation for C_{p,q/n} vs Random Regular networks under **beta = infinity** (SI model: permanent infection).

## Networks

| Label | p | q | Lattice edges | Random edges | Total degree |
|-------|---|---|---------------|--------------|--------------|
| p1q6  | 1 | 6 | 2             | ~6           | ~8           |
| p2q4  | 2 | 4 | 4             | ~4           | ~8           |
| p3q2  | 3 | 2 | 6             | ~2           | ~8           |

Baseline: Random Regular (all nodes exactly degree 8).

## Folder Structure

```
beta_infinite_delta_star/
├── core/                          # Core modules
│   ├── theory.py                  # Closed-form delta* (MGF + integral)
│   ├── simulation.py              # SI network simulation
│   └── plotting.py                # Heatmap visualization
├── sweep_simulation/              # Simulation parameter sweeps
│   └── sweep_p2_vs_thrshld.py
├── sweep_theory/                  # Theoretical parameter sweeps
│   └── sweep_p2_vs_thrshld.py
├── sbatch/                        # HPC job scripts
│   ├── submit_all.sh              # Master submission script
│   ├── sim/                       # Simulation sbatch files (3)
│   └── theory/                    # Theory sbatch files (3)
├── results/                       # Output data
│   ├── simulation/{p1q6,p2q4,p3q2}/
│   ├── theory/{p1q6,p2q4,p3q2}/
│   └── plots/
└── README.md
```

## Quick Start

### Theoretical computation (local, fast)

```bash
cd scripts/beta_infinite_delta_star

# Single point
python core/theory.py --p 1 --p2 0.8 --threshold 3.0

# Full sweep (one network)
python sweep_theory/sweep_p2_vs_thrshld.py --p 1
python sweep_theory/sweep_p2_vs_thrshld.py --p 2
python sweep_theory/sweep_p2_vs_thrshld.py --p 3
```

### Simulation (local, small scale)

```bash
cd scripts/beta_infinite_delta_star

# Single point (10 pairs for quick test)
python core/simulation.py --p 1 --p2 0.8 --threshold 3.0 --n-pairs 10

# Full sweep (slow, use HPC for 500 pairs)
python sweep_simulation/sweep_p2_vs_thrshld.py --p 1 --n-pairs 10 --n-workers 4
```

### HPC (SLURM)

```bash
cd scripts/beta_infinite_delta_star/sbatch

# Submit all 6 jobs (3 theory + 3 simulation)
bash submit_all.sh

# Submit theory only (fast, ~minutes)
bash submit_all.sh theory

# Submit simulation only (slow, ~hours)
bash submit_all.sh sim
```

## Theory

### Adoption Rule (Step-Function Threshold)

```
Pi(m) = 0    if m = 0
Pi(m) = p1   if 1 <= m < i
Pi(m) = p2   if m >= i
```

### Speed Functions (beta=infinity)

**Random Regular** (m ~ Binomial(8, alpha)):
```
g_R(alpha) = sum_{m=1}^{8} C(8,m) alpha^m (1-alpha)^{8-m} Pi(m)
```

**C_{p,q/n}** (m = m_l + r, m_l ~ Bin(2p,alpha), r ~ Pois(q*alpha)):
```
g_C(alpha) = sum_{m_l=0}^{2p} Bin(2p,m_l,alpha) * sum_r Pois(r,q*alpha) * Pi(m_l+r)
```

### Growth Rate

```
r = 2 * g(1/2)
```

### Delta* (MGF Approach)

```
delta* = A^{-3*r_H / pi^2}
```

where:
- A = (K - alpha_0) / alpha_0, alpha_0 = ceil(i)/N, K ~ 1
- r_H = 2*r_C*r_R / (r_C + r_R) (harmonic mean)

### Delta* (Integral Approach)

```
delta* ~ (a / (a+b))^{1/t_x}
```

Based on speed-function crossover analysis. Primary method for p3q2.

## Sweep Parameters

- **p2**: 0.50 to 1.00 (step 0.01) — 51 values
- **threshold**: 2.0 to 5.0 (step 0.1) — 31 values
- **Total**: 51 x 31 = 1,581 parameter combinations per network
- **Simulation**: 500 network pairs per combination
- **Networks**: p1q6, p2q4, p3q2

## Dependencies

- Python 3.8+
- numpy, pandas, scipy, networkx, matplotlib
- Optional: tqdm (progress bars)
