"""Two ways to estimate the failure probability P[g(X) < 0] of a limit state with independent basic variables:

form()        First-Order Reliability Method, Hasofer-Lind / Rackwitz-Fiessler iteration in standard normal space.
              Returns the reliability index beta, the design point and the sensitivity factors alpha.
importance()  Importance sampling centred on the FORM design point: checks the FORM result without the
              hundreds of millions of samples crude Monte Carlo would need at beta around 4 to 5.
crude_mc()    Plain Monte Carlo, used where the failure probability is large enough to count failures.
"""
import numpy as np
from scipy import stats


def _x(variables, u):
    return {name: d.x_from_u(u[..., i]) for i, (name, d) in enumerate(variables.items())}


def form(g, variables, tol=1e-5, max_iter=200, h=1e-5):
    names = list(variables)
    n = len(names)
    u = np.zeros(n)
    beta = 0.0
    for it in range(max_iter):
        g0 = float(g(_x(variables, u)))
        grad = np.empty(n)
        for i in range(n):                      # central differences in u-space
            du = np.zeros(n)
            du[i] = h
            grad[i] = (float(g(_x(variables, u + du))) - float(g(_x(variables, u - du)))) / (2 * h)
        norm = np.linalg.norm(grad)
        u_new = (grad @ u - g0) / norm ** 2 * grad   # HL-RF step
        if it >= 30:                                 # damp the step if plain HL-RF oscillates
            u_new = u + 0.5 * (u_new - u)
        beta_new = float(np.linalg.norm(u_new))
        if np.linalg.norm(u_new - u) < tol and abs(beta_new - beta) < tol:
            u = u_new
            beta = beta_new
            break
        u, beta = u_new, beta_new
    else:
        raise RuntimeError("FORM did not converge")
    alpha = -u / beta if beta > 0 else np.zeros(n)
    sign = 1.0 if float(g(_x(variables, np.zeros(n)))) > 0 else -1.0  # beta < 0 if the mean point already fails
    beta *= sign
    return {"beta": beta, "pf": float(stats.norm.cdf(-beta)), "u_star": u,
            "x_star": {k: float(v) for k, v in _x(variables, u).items()},
            "alpha": dict(zip(names, alpha.tolist())), "iterations": it + 1}


def importance(g, variables, u_star, n=200_000, seed=1):
    """Importance sampling with a unit-variance normal density centred on the design point."""
    rng = np.random.default_rng(seed)
    d = len(variables)
    v = rng.standard_normal((n, d)) + u_star
    fail = g(_x(variables, v)) < 0
    # weight = phi(v) / phi(v - u*)
    logw = -0.5 * np.sum(v ** 2, axis=1) + 0.5 * np.sum((v - u_star) ** 2, axis=1)
    w = np.where(fail, np.exp(logw), 0.0)
    pf = w.mean()
    se = w.std(ddof=1) / np.sqrt(n)
    return {"pf": float(pf), "beta": float(-stats.norm.ppf(pf)) if pf > 0 else np.inf,
            "cov_pf": float(se / pf) if pf > 0 else np.inf, "n": n}


def crude_mc(g, variables, n=2_000_000, seed=2, chunk=500_000):
    rng = np.random.default_rng(seed)
    fails, done = 0, 0
    while done < n:
        m = min(chunk, n - done)
        x = {name: d.sample(rng, m) for name, d in variables.items()}
        fails += int(np.sum(g(x) < 0))
        done += m
    pf = fails / n
    return {"pf": pf, "beta": float(-stats.norm.ppf(pf)) if pf > 0 else np.inf,
            "failures": fails, "n": n, "cov_pf": float(np.sqrt((1 - pf) / (n * pf))) if pf > 0 else np.inf}
