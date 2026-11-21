"""Run the whole study: reliability of Eurocode-designed members in modern and reclaimed steel, checks of FORM
against sampling, and the partial factor reclaimed steel would need.

    python3 run.py          writes results/*.csv, results/summary.json and figures/*.png
"""
import csv
import json
import os
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from reliability.methods import crude_mc, form, importance
from reliability.models import (BETA_TARGET_50Y, CHI, CHI_TYPICAL, LIMIT_STATES, MODERN, SCENARIOS)

HERE = os.path.dirname(os.path.abspath(__file__))
RES, FIG = os.path.join(HERE, "results"), os.path.join(HERE, "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
CHI_CAL = [c for c in CHI if CHI_TYPICAL[0] - 1e-9 <= c <= CHI_TYPICAL[1] + 1e-9]


def beta(ls, steel, chi, gamma_m=None):
    gm = ls.gamma_m_default if gamma_m is None else gamma_m
    return form(ls.g(chi, gm), ls.variables(steel, chi))


def mean_beta(ls, steel, gamma_m=None):
    return float(np.mean([beta(ls, steel, c, gamma_m)["beta"] for c in CHI_CAL]))


def min_beta(ls, steel, gamma_m=None):
    return float(min(beta(ls, steel, c, gamma_m)["beta"] for c in CHI_CAL))


def main():
    rows, checks, summary = [], [], {"beta_target_50y": BETA_TARGET_50Y, "chi_calibration": CHI_CAL}

    # 1. Reliability index against load ratio, FORM checked by importance sampling at every point
    for ls in LIMIT_STATES:
        for st in SCENARIOS:
            for chi in CHI:
                f = beta(ls, st, chi)
                ims = importance(ls.g(chi, ls.gamma_m_default), ls.variables(st, chi), f["u_star"], n=100_000)
                rows.append({"limit_state": ls.key, "steel": st.name, "chi": float(chi),
                             "gamma_M": ls.gamma_m_default, "beta_FORM": round(f["beta"], 3),
                             "beta_IS": round(ims["beta"], 3), "pf_FORM": f"{f['pf']:.3e}",
                             "pf_IS": f"{ims['pf']:.3e}", "cov_IS": round(ims["cov_pf"], 3)})
    with open(os.path.join(RES, "beta_vs_chi.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    diff = max(abs(r["beta_FORM"] - r["beta_IS"]) for r in rows)
    summary["max_abs_beta_difference_FORM_vs_IS"] = round(diff, 3)

    # crude Monte Carlo where failures are common enough to count (the least reliable point)
    worst = min(rows, key=lambda r: r["beta_FORM"])
    ls = next(l for l in LIMIT_STATES if l.key == worst["limit_state"])
    st = next(s for s in SCENARIOS if s.name == worst["steel"])
    mc = crude_mc(ls.g(worst["chi"], ls.gamma_m_default), ls.variables(st, worst["chi"]), n=4_000_000)
    summary["crude_mc_check"] = {"limit_state": ls.key, "steel": st.name, "chi": worst["chi"],
                                 "beta_FORM": worst["beta_FORM"], "beta_IS": worst["beta_IS"],
                                 "beta_MC": round(mc["beta"], 3), "failures": mc["failures"], "n": mc["n"]}

    # 2. Mean and minimum beta over typical load ratios, and the partial factor reclaimed steel would need
    cal = []
    for ls in LIMIT_STATES:
        ref = mean_beta(ls, MODERN)
        for st in SCENARIOS:
            entry = {"limit_state": ls.key, "steel": st.name, "gamma_M_recommended": ls.gamma_m_default,
                     "mean_beta": round(mean_beta(ls, st), 3), "min_beta": round(min_beta(ls, st), 3)}
            if st is not MODERN:
                entry["gamma_M_same_reliability_as_modern"] = round(
                    brentq(lambda gm: mean_beta(ls, st, gm) - ref, 0.9, 2.2, xtol=1e-4), 3)
            entry["gamma_M_for_min_beta_3.8"] = round(
                brentq(lambda gm: min_beta(ls, st, gm) - BETA_TARGET_50Y, 0.9, 2.2, xtol=1e-4), 3)
            cal.append(entry)
    with open(os.path.join(RES, "calibration.csv"), "w", newline="") as fh:
        keys = sorted({k for e in cal for k in e}, key=lambda k: list(cal[1]).index(k) if k in cal[1] else 99)
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(cal)
    summary["calibration"] = cal

    # 3. Which variables drive failure: FORM sensitivity factors at chi = 0.4, bending
    ls = LIMIT_STATES[0]
    sens = {st.name: {k: round(v, 3) for k, v in beta(ls, st, 0.4)["alpha"].items()} for st in SCENARIOS}
    summary["alpha_bending_chi_0.4"] = sens

    # 4. Required gamma_M0 for equal reliability, against the scatter of fy (bending), for two mean ratios
    covs = np.round(np.arange(0.05, 0.151, 0.01), 3)
    ref = mean_beta(ls, MODERN)
    sweep = {}
    for mean in (1.15, 1.05):
        sweep[mean] = [brentq(lambda gm: mean_beta(ls, replace(MODERN, fy_mean=mean, fy_cov=c), gm) - ref,
                              0.8, 2.2, xtol=1e-4) for c in covs]
    summary["gamma_M0_vs_fy_cov"] = {str(m): dict(zip(map(float, covs), map(lambda v: round(v, 3), vals)))
                                      for m, vals in sweep.items()}

    with open(os.path.join(RES, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    # figures
    colours = {MODERN.name: "#1f4e79", SCENARIOS[1].name: "#d08c00", SCENARIOS[2].name: "#b0302a"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, ls in zip(axes, LIMIT_STATES):
        for st in SCENARIOS:
            b = [r["beta_FORM"] for r in rows if r["limit_state"] == ls.key and r["steel"] == st.name]
            ax.plot(CHI, b, marker="o", ms=4, lw=2, color=colours[st.name], label=st.name)
        ax.axhline(BETA_TARGET_50Y, color="0.4", ls="--", lw=1)
        ax.text(CHI[0], BETA_TARGET_50Y + 0.06, "EN 1990 target, 50 years (3.8)", fontsize=8, color="0.3")
        ax.set_title(ls.title.split(" (")[0] + f"\n(gamma_M = {ls.gamma_m_default})", fontsize=10)
        ax.set_xlabel("load ratio  chi = Qk / (Gk + Qk)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("reliability index beta (FORM)")
    axes[1].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "1_beta_vs_load_ratio.png"), dpi=160)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for mean, vals in sweep.items():
        ax.plot(covs, vals, marker="o", ms=4, lw=2, label=f"mean fy / fyk = {mean}")
    ax.axhline(1.0, color="0.4", ls="--", lw=1)
    ax.set_xlabel("coefficient of variation of fy")
    ax.set_ylabel("gamma_M0 for the same reliability as modern steel")
    ax.set_title("Bending: partial factor needed for equal reliability", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "2_gamma_M0_vs_scatter.png"), dpi=160)

    fig, ax = plt.subplots(figsize=(6, 3.6))
    names = list(next(iter(sens.values())))
    width = 0.27
    for i, st in enumerate(SCENARIOS):
        ax.bar(np.arange(len(names)) + (i - 1) * width, [sens[st.name][n] ** 2 for n in names], width,
               color=colours[st.name], label=st.name)
    ax.set_xticks(np.arange(len(names)), names)
    ax.set_ylabel("share of variance at the design point (alpha^2)")
    ax.set_title("Bending, chi = 0.4: what drives failure", fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "3_sensitivity.png"), dpi=160)

    print(json.dumps({k: v for k, v in summary.items() if k != "gamma_M0_vs_fy_cov"}, indent=1))


if __name__ == "__main__":
    main()
