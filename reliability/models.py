"""Probabilistic models, steel scenarios and limit states.

Everything is normalised: Gk + Qk = 1 and the load ratio chi = Qk / (Gk + Qk). Strengths are expressed as
ratios to their characteristic values (fy / fyk, fu / fuk), so one model serves any steel grade.

Each member is designed exactly to the Eurocode, at the limit, with the 6.10 combination:
    Ed = gamma_G * Gk + gamma_Q * Qk
    bending, Class 1 or 2 section (EN 1993-1-1, 6.2.5):   Wpl * fyk / gamma_M0            = Ed
    net section of a bolted tension member (6.2.3):         0.9 * A_net * fuk / gamma_M2   = Ed
and its reliability is then computed with the probabilistic models below.

The models for modern steel and the actions follow values commonly used in Eurocode calibration studies
(JCSS Probabilistic Model Code; Nadolski, Rozsas and Sykora 2019, Table 3). The reclaimed-steel scenarios
are ILLUSTRATIVE assumptions, not measured data: they show how the method responds to a lower mean and a
larger scatter of strength, and are meant to be replaced by test statistics of real reclaimed steel.
"""
from dataclasses import dataclass, field

import numpy as np

from .distributions import Gumbel, Lognormal, Normal

GAMMA_G, GAMMA_Q = 1.35, 1.5          # EN 1990, 6.10 (recommended values)
GAMMA_M0, GAMMA_M2 = 1.0, 1.25        # EN 1993-1-1, recommended values
BETA_TARGET_50Y = 3.8                 # EN 1990 Annex C, RC2, 50-year reference period


@dataclass
class Steel:
    name: str
    fy_mean: float      # mean of fy / fyk
    fy_cov: float
    fu_mean: float      # mean of fu / fuk
    fu_cov: float
    geo_cov: float = 0.02   # section property (Wpl or A_net) relative to nominal, mean 1.0
    note: str = ""


MODERN = Steel("modern S355", 1.15, 0.07, 1.10, 0.05, 0.02,
               "JCSS-type model: lognormal, fyk about two standard deviations below the mean")
RECLAIMED_R1 = Steel("reclaimed, grade confirmed by tests (R1)", 1.15, 0.10, 1.10, 0.07, 0.02,
                     "illustrative: same mean as modern steel, larger scatter")
RECLAIMED_R2 = Steel("reclaimed, older steel of uncertain grade (R2)", 1.05, 0.12, 1.05, 0.10, 0.04,
                     "illustrative: lower mean, larger scatter, more geometric variation (wear, tolerances)")
SCENARIOS = [MODERN, RECLAIMED_R1, RECLAIMED_R2]


def actions(chi):
    """Permanent, variable (50-year maximum, imposed load) and load-effect model uncertainty."""
    gk, qk = 1.0 - chi, chi
    return {
        "G": Normal(gk, 0.10) if gk > 0 else Normal(1e-9, 0.10),
        "Q": Gumbel(0.6 * qk, 0.35) if qk > 0 else Normal(1e-9, 0.35),
        "KE": Lognormal(1.0, 0.10),
    }


def design_load(chi):
    return GAMMA_G * (1.0 - chi) + GAMMA_Q * chi


@dataclass
class LimitState:
    key: str
    title: str
    strength: str           # "fy" or "fu"
    gamma_m_default: float
    design_factor: float    # resistance at the design point = design_factor * section * f_k / gamma_M

    def variables(self, steel, chi):
        s_mean, s_cov = ((steel.fy_mean, steel.fy_cov) if self.strength == "fy" else (steel.fu_mean, steel.fu_cov))
        v = {"KR": Lognormal(1.0, 0.05), "geo": Normal(1.0, steel.geo_cov), "f": Lognormal(s_mean, s_cov)}
        v.update(actions(chi))
        return v

    def g(self, chi, gamma_m):
        # nominal section chosen so that design resistance equals design load: section * f_k = gamma_M * Ed / factor
        r_nom = gamma_m * design_load(chi) / self.design_factor

        def g(x):
            return x["KR"] * x["geo"] * r_nom * x["f"] - x["KE"] * (x["G"] + x["Q"])
        return g


BENDING = LimitState("bending", "Bending resistance of a Class 1 beam (EN 1993-1-1, 6.2.5)", "fy", GAMMA_M0, 1.0)
NET_SECTION = LimitState("net_section", "Net-section failure of a bolted tension member (EN 1993-1-1, 6.2.3)",
                         "fu", GAMMA_M2, 0.9)
LIMIT_STATES = [BENDING, NET_SECTION]

CHI = np.round(np.arange(0.1, 0.81, 0.1), 2)
CHI_TYPICAL = (0.2, 0.6)   # load ratios typical of steel buildings, used for calibration
