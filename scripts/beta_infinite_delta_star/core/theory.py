"""
Closed-Form Delta-Star Calculator for beta = infinity (SI Model)

Implements the theoretical delta* formulas from:
  - closed_form_delta_star_beta_inf_p1q6_final.tex
  - closed_form_delta_star_beta_inf_p2q4_final.tex
  - closed_form_delta_star_beta_inf_p3q2_final.tex

Unified implementation: a single algorithm handles all three q>0 networks
(p1q6, p2q4, p3q2) with only p and q as parameters.

Two approaches:
  1. MGF/Logistic: delta* = A^{-3*r_H / pi^2}
  2. Speed-Function Integral: delta* ~ (a/(a+b))^{1/t_x}

Under beta=infinity (SI model), the system is memoryless:
  - All infected nodes remain infectious forever
  - State reduces to (s_t) only
  - Speed function: c_t = s_t * g(alpha_t)
  - Growth follows logistic curve (for networks with random edges)
"""

import numpy as np
from math import factorial, ceil, pi
from scipy import integrate

# Poisson truncation limit
R_MAX = 30


# ================================================================
# Combinatorial utilities
# ================================================================

def comb(n, k):
    """Binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))


def poisson_pmf(r, lam):
    """Poisson PMF: P(R=r | lambda)."""
    if lam <= 0:
        return 1.0 if r == 0 else 0.0
    if r < 0:
        return 0.0
    log_p = r * np.log(lam) - lam - sum(np.log(j) for j in range(1, r + 1))
    return np.exp(log_p)


def binomial_pmf(k, n, p):
    """Binomial PMF: P(X=k | n, p)."""
    if p <= 0:
        return 1.0 if k == 0 else 0.0
    if p >= 1:
        return 1.0 if k == n else 0.0
    return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))


# ================================================================
# Adoption rule Pi(m) -- Step-function threshold
# From tex: Pi(m) = 0 if m=0, p1 if 1<=m<i, p2 if m>=i
# ================================================================

def Pi(m, p1, p2, threshold):
    """
    Step-function adoption probability given m infectious neighbors.

    Pi(m) = 0   if m = 0
    Pi(m) = p1  if 1 <= m < threshold
    Pi(m) = p2  if m >= threshold
    """
    if m <= 0:
        return 0.0
    i = int(ceil(threshold))
    if m < i:
        return p1
    return p2


# ================================================================
# Speed function g_R(alpha) -- Random Regular, degree k=8
# g_R(alpha) = sum_{m=1}^{8} C(8,m) alpha^m (1-alpha)^{8-m} Pi(m)
# ================================================================

def g_random(alpha, p1, p2, threshold, k=8):
    """
    Speed function for Random Regular network with degree k.
    Each susceptible has k neighbors, each infectious with prob alpha.
    m ~ Binomial(k, alpha).
    """
    if alpha <= 0 or alpha >= 1:
        return 0.0
    g = 0.0
    for m in range(1, k + 1):
        binom_p = binomial_pmf(m, k, alpha)
        g += binom_p * Pi(m, p1, p2, threshold)
    return g


# ================================================================
# Speed function g_C(alpha) -- C_{p,q/n} network (bulk, beta=inf)
# M = m_l + r, m_l ~ Bin(2p, alpha), r ~ Pois(q*alpha)
# g_C(alpha) = sum_{m_l} Bin(2p,m_l,alpha) * sum_r Pois(r,q*alpha) * Pi(m_l+r)
# ================================================================

def g_cpq(alpha, p1, p2, threshold, p, q):
    """
    Speed function for C_{p,q/n} network under beta=infinity.

    Contact distribution: M = m_l + r
      m_l ~ Binomial(2p, alpha)  -- lattice contacts
      r   ~ Poisson(q * alpha)   -- random contacts

    Parameters:
        alpha: infected fraction
        p1, p2: adoption probabilities
        threshold: adoption threshold i
        p: cycle-power parameter (lattice: 2p edges)
        q: random component (expected q random edges)
    """
    if alpha <= 0 or alpha >= 1:
        return 0.0
    n_lat = 2 * p
    lam = q * alpha
    g = 0.0
    for m_l in range(n_lat + 1):
        bp = binomial_pmf(m_l, n_lat, alpha)
        if bp < 1e-18:
            continue
        for r in range(R_MAX + 1):
            pp = poisson_pmf(r, lam)
            if pp < 1e-15 and r > 0:
                break
            m_total = m_l + r
            if m_total > 0:
                g += bp * pp * Pi(m_total, p1, p2, threshold)
    return g


# ================================================================
# Growth rate: r = 2 * g(1/2)
# From logistic model at inflection point alpha = K/2 ~ 1/2
# ================================================================

def growth_rate(g_half):
    """Growth rate from speed function at alpha=1/2."""
    return 2.0 * g_half


def harmonic_mean(r1, r2):
    """Harmonic mean of two growth rates."""
    if r1 <= 0 or r2 <= 0:
        return 0.0
    return 2.0 * r1 * r2 / (r1 + r2)


# ================================================================
# Approach 1: MGF/Logistic closed form
# delta* = A^{-3*r_H / pi^2}
# A = (K - alpha_0) / alpha_0
# ================================================================

def delta_star_mgf(r_C, r_R, alpha_0, K=1.0):
    """
    MGF approach: delta* = A^{-3*r_H / pi^2}

    Parameters:
        r_C: growth rate of C_{p,q/n}
        r_R: growth rate of Random Regular
        alpha_0: initial infected fraction (seeds/N)
        K: final infected fraction (approx 1 for beta=inf)

    Returns:
        delta_star: float, or None if undefined
    """
    if r_C <= 0 or r_R <= 0:
        return None
    if alpha_0 <= 0 or alpha_0 >= K:
        return None

    A = (K - alpha_0) / alpha_0
    if A <= 0:
        return None

    r_H = harmonic_mean(r_C, r_R)
    exponent = -3.0 * r_H / (pi ** 2)
    delta_star = A ** exponent

    if delta_star < 0 or delta_star > 1:
        return None
    return delta_star


# ================================================================
# Approach 2: Speed-Function Integral Transform
# delta* ~ (a/(a+b))^{1/t_x}
# ================================================================

def find_speed_crossover(g_C_func, g_R_func, alpha_0, alpha_f, n_points=10000):
    """
    Find alpha_x where g_C(alpha_x) = g_R(alpha_x).

    Returns alpha_x or None if no crossover found.
    """
    alphas = np.linspace(alpha_0, alpha_f, n_points)
    diff_prev = g_C_func(alphas[0]) - g_R_func(alphas[0])

    for i in range(1, len(alphas)):
        diff_curr = g_C_func(alphas[i]) - g_R_func(alphas[i])
        if diff_prev * diff_curr < 0:
            # Linear interpolation
            a1, a2 = alphas[i - 1], alphas[i]
            alpha_x = a1 - diff_prev * (a2 - a1) / (diff_curr - diff_prev)
            return alpha_x
        diff_prev = diff_curr
    return None


def delta_star_integral(g_C_func, g_R_func, alpha_0, alpha_f, N,
                        n_points=5000):
    """
    Speed-function integral approach.

    delta* ~ (a/(a+b))^{1/t_x}

    Handles both directions:
      - C faster early, R faster late (typical for beta=inf)
      - R faster early, C faster late

    The formula uses the absolute magnitude of early and late advantages.

    Returns:
        dict with 'delta_star', 'alpha_x', 't_x', 'early_adv', 'late_adv',
        'T_total', 'early_leader'
    """
    alpha_x = find_speed_crossover(g_C_func, g_R_func, alpha_0, alpha_f)

    if alpha_x is None:
        g_C_mid = g_C_func(0.5)
        g_R_mid = g_R_func(0.5)
        if g_C_mid > g_R_mid:
            return {'delta_star': None, 'case': 'A', 'reason': 'C always faster'}
        else:
            return {'delta_star': None, 'case': 'B', 'reason': 'R always faster'}

    def integrand_time_avg(alpha):
        g_bar = (g_C_func(alpha) + g_R_func(alpha)) / 2.0
        if g_bar <= 0:
            return 0.0
        return 1.0 / ((1.0 - alpha) * g_bar)

    def integrand_diff(alpha):
        """(g_C - g_R) / g_bar — positive when C is faster."""
        g_bar = (g_C_func(alpha) + g_R_func(alpha)) / 2.0
        if g_bar <= 0:
            return 0.0
        return (g_C_func(alpha) - g_R_func(alpha)) / g_bar

    # t_x: time to crossover
    t_x, _ = integrate.quad(integrand_time_avg, alpha_0, alpha_x,
                            limit=200, epsabs=1e-10)

    # T: total time
    T_total, _ = integrate.quad(integrand_time_avg, alpha_0, alpha_f,
                                limit=200, epsabs=1e-10)

    # Early phase integral (alpha_0 to alpha_x)
    early_raw, _ = integrate.quad(integrand_diff, alpha_0, alpha_x,
                                  limit=200, epsabs=1e-10)
    early_adv = N * abs(early_raw)

    # Late phase integral (alpha_x to alpha_f)
    late_raw, _ = integrate.quad(integrand_diff, alpha_x, alpha_f,
                                 limit=200, epsabs=1e-10)
    late_adv = N * abs(late_raw)

    # Determine which network leads early
    early_leader = 'C' if early_raw > 0 else 'R'

    if t_x <= 0 or T_total <= t_x:
        return {'delta_star': None, 'case': 'Invalid',
                'reason': 'Time computation error'}

    a = early_adv / t_x
    b = late_adv / (T_total - t_x)

    if a + b <= 0 or a <= 0:
        return {'delta_star': None, 'case': 'Invalid',
                'reason': 'Deficit/surplus error'}

    ratio = a / (a + b)
    if ratio <= 0 or ratio >= 1:
        return {'delta_star': None, 'case': 'Invalid',
                'reason': f'Invalid ratio {ratio}'}

    delta_star = ratio ** (1.0 / t_x)

    return {
        'delta_star': delta_star,
        'case': 'C',
        'alpha_x': alpha_x,
        't_x': t_x,
        'T_total': T_total,
        'early_adv': early_adv,
        'late_adv': late_adv,
        'early_leader': early_leader,
    }


# ================================================================
# Main entry point: compute delta* for given parameters
# ================================================================

def compute_delta_star(p, q, p1, p2, threshold, N=2000, k=8,
                       method='both', alpha_f=0.999):
    """
    Compute theoretical delta* for C_{p,q/n} vs Random Regular
    under beta = infinity (SI model).

    Parameters:
        p: cycle-power parameter (1, 2, or 3)
        q: random component (6, 4, or 2)
        p1: base adoption probability
        p2: reinforced adoption probability
        threshold: adoption threshold i
        N: network size
        k: degree of Random Regular
        method: 'mgf', 'integral', or 'both'
        alpha_f: final infection fraction for integral method

    Returns:
        dict with keys:
          'delta_star_mgf': float or None
          'delta_star_integral': float or None
          'delta_star': float or None (best estimate)
          'case': 'A', 'B', 'C', or 'Invalid'
          'g_R_half', 'g_C_half': speed functions at alpha=1/2
          'r_R', 'r_C', 'r_H': growth rates
    """
    seeds = int(ceil(threshold))
    alpha_0 = seeds / N

    # Speed functions at alpha = 1/2
    g_R_half = g_random(0.5, p1, p2, threshold, k)
    g_C_half = g_cpq(0.5, p1, p2, threshold, p, q)

    # Growth rates
    r_R = growth_rate(g_R_half)
    r_C = growth_rate(g_C_half)
    r_H = harmonic_mean(r_C, r_R)

    result = {
        'g_R_half': g_R_half,
        'g_C_half': g_C_half,
        'r_R': r_R,
        'r_C': r_C,
        'r_H': r_H,
        'alpha_0': alpha_0,
        'seeds': seeds,
    }

    # MGF approach
    ds_mgf = None
    if method in ('mgf', 'both'):
        ds_mgf = delta_star_mgf(r_C, r_R, alpha_0)
        result['delta_star_mgf'] = ds_mgf

    # Integral approach
    ds_int = None
    if method in ('integral', 'both'):
        g_C_func = lambda a: g_cpq(a, p1, p2, threshold, p, q)
        g_R_func = lambda a: g_random(a, p1, p2, threshold, k)
        int_result = delta_star_integral(g_C_func, g_R_func,
                                         alpha_0, alpha_f, N)
        ds_int = int_result.get('delta_star')
        result['delta_star_integral'] = ds_int
        result['integral_details'] = int_result

    # Determine best estimate and case
    if ds_mgf is not None or ds_int is not None:
        result['case'] = 'C'
        # For beta=infinity, MGF is the primary method for all networks.
        # The integral approach has poor accuracy when the speed crossover
        # occurs very early (alpha_x < 0.15), which is typical for beta=inf.
        if ds_mgf is not None:
            result['delta_star'] = ds_mgf
        else:
            result['delta_star'] = ds_int
    else:
        # Check dominance direction
        if g_C_half > g_R_half:
            result['case'] = 'A'  # C always dominates
        elif g_C_half < g_R_half:
            result['case'] = 'B'  # Random always dominates
        else:
            result['case'] = 'Invalid'
        result['delta_star'] = None

    return result


# ================================================================
# Batch computation for parameter sweeps
# ================================================================

def compute_sweep_p2_vs_thrshld(p, q, p1=0.1, N=2000,
                                 p2_values=None, thrshld_values=None):
    """
    Compute theoretical delta* for a grid of (p2, threshold) values.

    Returns:
        list of dicts, one per parameter combination
    """
    if p2_values is None:
        p2_values = np.arange(0.50, 1.01, 0.01)
    if thrshld_values is None:
        thrshld_values = np.arange(2.0, 5.1, 0.1)

    results = []
    for thrshld in thrshld_values:
        for p2 in p2_values:
            if p1 > p2:
                continue
            r = compute_delta_star(p, q, p1, p2, thrshld, N)
            results.append({
                'p': p,
                'q': q,
                'p1': p1,
                'p2': round(float(p2), 4),
                'thrshld': round(float(thrshld), 2),
                'N': N,
                'seeds': r['seeds'],
                'g_R_half': r['g_R_half'],
                'g_C_half': r['g_C_half'],
                'r_R': r['r_R'],
                'r_C': r['r_C'],
                'r_H': r['r_H'],
                'case': r['case'],
                'delta_star': r.get('delta_star'),
                'delta_star_mgf': r.get('delta_star_mgf'),
                'delta_star_integral': r.get('delta_star_integral'),
            })
    return results


# ================================================================
# CLI
# ================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Theoretical delta* for C_{p,q/n} vs Random Regular (beta=inf)')
    parser.add_argument('--p', type=int, required=True,
                        help='Cycle-power parameter (1, 2, or 3)')
    parser.add_argument('--q', type=float, default=None,
                        help='Random component (default: 8-2p)')
    parser.add_argument('--p1', type=float, default=0.1)
    parser.add_argument('--p2', type=float, default=1.0)
    parser.add_argument('--threshold', type=float, default=2.0)
    parser.add_argument('--N', type=int, default=2000)
    args = parser.parse_args()

    q = args.q if args.q is not None else float(8 - 2 * args.p)

    print("=" * 60)
    print(f"Theoretical delta* for C_{{{args.p},{q:.0f}/n}} vs Random Regular")
    print(f"  beta = infinity (SI model)")
    print("=" * 60)
    print(f"  N={args.N}, p1={args.p1}, p2={args.p2}, threshold={args.threshold}")
    print(f"  seeds = {int(ceil(args.threshold))}")
    print("=" * 60)

    result = compute_delta_star(args.p, q, args.p1, args.p2,
                                args.threshold, args.N)

    print(f"\n  g_R(1/2) = {result['g_R_half']:.6f}")
    print(f"  g_C(1/2) = {result['g_C_half']:.6f}")
    print(f"  r_R      = {result['r_R']:.6f}")
    print(f"  r_C      = {result['r_C']:.6f}")
    print(f"  r_H      = {result['r_H']:.6f}")
    print(f"  alpha_0  = {result['alpha_0']:.6f}")
    print(f"\n  Case: {result['case']}")
    if result.get('delta_star_mgf') is not None:
        print(f"  delta* (MGF)      = {result['delta_star_mgf']:.6f}")
    if result.get('delta_star_integral') is not None:
        print(f"  delta* (Integral) = {result['delta_star_integral']:.6f}")
    if result.get('delta_star') is not None:
        print(f"  delta* (best)     = {result['delta_star']:.6f}")
    print("=" * 60)
