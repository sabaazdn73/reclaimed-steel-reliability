"""Marginal distributions defined by mean and coefficient of variation, with the two maps FORM needs:
u -> x (standard normal space to physical space) and x -> u."""
import numpy as np
from scipy import stats

EULER_GAMMA = 0.5772156649015329


class Normal:
    def __init__(self, mean, cov):
        self.mean, self.cov = mean, cov
        self.sd = abs(mean) * cov

    def x_from_u(self, u):
        return self.mean + self.sd * u

    def u_from_x(self, x):
        return (x - self.mean) / self.sd

    def sample(self, rng, n):
        return rng.normal(self.mean, self.sd, n)


class Lognormal:
    def __init__(self, mean, cov):
        self.mean, self.cov = mean, cov
        self.zeta = np.sqrt(np.log(1 + cov ** 2))
        self.lam = np.log(mean) - 0.5 * self.zeta ** 2

    def x_from_u(self, u):
        return np.exp(self.lam + self.zeta * u)

    def u_from_x(self, x):
        return (np.log(x) - self.lam) / self.zeta

    def sample(self, rng, n):
        return rng.lognormal(self.lam, self.zeta, n)


class Gumbel:
    """Largest-value type I (maxima), e.g. the 50-year maximum of a variable action."""

    def __init__(self, mean, cov):
        self.mean, self.cov = mean, cov
        sd = mean * cov
        self.beta = sd * np.sqrt(6) / np.pi          # scale
        self.mu = mean - EULER_GAMMA * self.beta     # location

    def x_from_u(self, u):
        # -log(Phi(u)) written as -log1p(-sf(u)) so the upper tail stays accurate for large u
        t = -np.log1p(-stats.norm.sf(u))
        return self.mu - self.beta * np.log(t)

    def u_from_x(self, x):
        sf = -np.expm1(-np.exp(-(x - self.mu) / self.beta))   # 1 - F(x), accurate in the upper tail
        return stats.norm.isf(sf)

    def sample(self, rng, n):
        return rng.gumbel(self.mu, self.beta, n)


def characteristic(dist, fractile):
    """Value of the distribution at a given fractile (e.g. 0.05 for a characteristic strength)."""
    return float(dist.x_from_u(stats.norm.ppf(fractile)))
