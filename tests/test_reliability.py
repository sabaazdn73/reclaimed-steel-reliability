"""Checks against cases with an exact answer.   python3 -m unittest discover -s tests -v"""
import unittest

import numpy as np

from reliability.distributions import Gumbel, Lognormal, Normal, characteristic
from reliability.methods import crude_mc, form, importance
from reliability.models import BENDING, MODERN, NET_SECTION, RECLAIMED_R2


class ExactCases(unittest.TestCase):
    def test_linear_normal(self):
        # g = R - S, R ~ N(10, 1), S ~ N(5, 1): beta = 5 / sqrt(2) exactly
        v = {"R": Normal(10, 0.1), "S": Normal(5, 0.2)}
        r = form(lambda x: x["R"] - x["S"], v)
        self.assertAlmostEqual(r["beta"], 5 / np.sqrt(2), places=4)

    def test_lognormal_ratio(self):
        # g = ln R - ln S with lognormal R and S is linear in u: beta = (lam_R - lam_S) / sqrt(zeta_R^2 + zeta_S^2)
        R, S = Lognormal(2.0, 0.1), Lognormal(1.0, 0.2)
        r = form(lambda x: np.log(x["R"]) - np.log(x["S"]), {"R": R, "S": S})
        exact = (R.lam - S.lam) / np.hypot(R.zeta, S.zeta)
        self.assertAlmostEqual(r["beta"], exact, places=4)

    def test_alpha_is_a_unit_vector_with_signs(self):
        r = form(lambda x: x["R"] - x["S"], {"R": Normal(10, 0.1), "S": Normal(5, 0.2)})
        self.assertAlmostEqual(np.hypot(*r["alpha"].values()), 1.0, places=6)
        self.assertGreater(r["alpha"]["R"], 0)   # resistance: positive alpha
        self.assertLess(r["alpha"]["S"], 0)      # load: negative alpha


class Distributions(unittest.TestCase):
    def test_mean_and_cov_by_sampling(self):
        rng = np.random.default_rng(0)
        for d in (Normal(3, 0.1), Lognormal(1.15, 0.07), Gumbel(0.6, 0.35)):
            x = d.sample(rng, 400_000)
            self.assertAlmostEqual(x.mean() / d.mean, 1, places=2)
            self.assertAlmostEqual(x.std() / x.mean() / d.cov, 1, places=1)

    def test_u_x_round_trip(self):
        u = np.linspace(-5, 8, 27)
        for d in (Normal(3, 0.1), Lognormal(1.15, 0.07), Gumbel(0.6, 0.35)):
            np.testing.assert_allclose(d.u_from_x(d.x_from_u(u)), u, atol=1e-6)

    def test_characteristic_yield_is_near_the_lower_tail(self):
        # with mean fy/fyk = 1.15 and CoV 0.07, fyk lies close to the 2 to 5% fractile
        fy = Lognormal(MODERN.fy_mean, MODERN.fy_cov)
        self.assertTrue(characteristic(fy, 0.02) < 1.0 < characteristic(fy, 0.06))


class FormAgainstSampling(unittest.TestCase):
    def test_form_matches_importance_sampling(self):
        for ls, steel, chi in ((BENDING, MODERN, 0.4), (NET_SECTION, RECLAIMED_R2, 0.3)):
            g, v = ls.g(chi, ls.gamma_m_default), ls.variables(steel, chi)
            f = form(g, v)
            s = importance(g, v, f["u_star"], n=100_000)
            self.assertLess(abs(f["beta"] - s["beta"]), 0.08)

    def test_form_matches_crude_monte_carlo_where_failures_are_common(self):
        g, v = BENDING.g(0.1, 1.0), BENDING.variables(RECLAIMED_R2, 0.1)
        f, mc = form(g, v), crude_mc(g, v, n=1_000_000)
        self.assertGreater(mc["failures"], 5_000)
        self.assertLess(abs(f["beta"] - mc["beta"]), 0.05)


class Design(unittest.TestCase):
    def test_designed_at_the_limit(self):
        # with the resistance at its design value (fk * factor / gamma_M) and the design loads, g equals zero
        from reliability.models import GAMMA_G, GAMMA_Q
        for ls in (BENDING, NET_SECTION):
            chi = 0.4
            g = ls.g(chi, ls.gamma_m_default)
            f_design = ls.design_factor / ls.gamma_m_default
            x = {"KR": 1.0, "geo": 1.0, "f": f_design,
                 "KE": 1.0, "G": GAMMA_G * (1 - chi), "Q": GAMMA_Q * chi}
            self.assertAlmostEqual(float(g(x)), 0.0, places=9)


if __name__ == "__main__":
    unittest.main()
