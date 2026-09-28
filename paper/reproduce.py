#!/usr/bin/env python3
"""Reproduce the figures and numerical checks in the accompanying paper.

Run from this directory: python3 reproduce.py
Requires numpy, scipy and matplotlib. All three figures are vector PDFs.
"""

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


OUT = Path(__file__).resolve().parent
FIGURES = OUT / "figures"
G = 6.67430e-11
C = 299792458.0
H = 6.62607015e-34
HBAR = H / (2.0 * math.pi)
KB = 1.380649e-23
M0 = 1.0e7
RHO = 1.0e3
CS = 1500.0
LAM = 0.25
YEAR = 365.25 * 86400.0

B = HBAR * C**4 / (15360.0 * math.pi * G**2)
K = B * C**2
A = 4.0 * math.pi * LAM * RHO * G**2 / CS**3
MC = (B / A)**0.25
TSTAR = 1.0 / (A * MC)
TAU_BB = M0**3 / (3.0 * B)


def primitive(x):
    """Antiderivative of x**2/(x**4-1); x must be positive and != 1."""
    return 0.25 * np.log(np.abs((x - 1.0) / (x + 1.0))) + 0.5 * np.arctan(x)


def plot_style():
    plt.rcParams.update({
        "font.size": 11.5,
        "axes.labelsize": 12.0,
        "xtick.labelsize": 11.0,
        "ytick.labelsize": 11.0,
        "legend.fontsize": 10.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
    })


def figure_one():
    x = np.linspace(0.5, 1.75, 501)
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.axvspan(0.5, 1.0, color="#c74436", alpha=0.055)
    ax.axvspan(1.0, 1.75, color="#1675ae", alpha=0.06)
    ax.plot(x, x**2 - x**-2, color="#263846", lw=2.3)
    ax.axhline(0, color="#71808a", lw=0.9)
    ax.axvline(1, color="#71808a", lw=0.9, ls=":")
    ax.annotate("", xy=(0.62, 0), xytext=(0.91, 0),
                arrowprops={"arrowstyle": "-|>", "color": "#c74436", "lw": 2})
    ax.annotate("", xy=(1.52, 0), xytext=(1.09, 0),
                arrowprops={"arrowstyle": "-|>", "color": "#1675ae", "lw": 2})
    ax.plot([1], [0], "o", markersize=7, markerfacecolor="white",
            markeredgecolor="#263846", markeredgewidth=1.6, zorder=5)
    ax.annotate("Unstable balance", xy=(1, 0), xytext=(0.55, 1.38),
                arrowprops={"arrowstyle": "-", "color": "#71808a"}, fontsize=11)
    ax.text(0.70, -3.2, "Mass decreases", color="#c74436")
    ax.text(1.18, 2.47, "Mass increases", color="#1675ae")
    ax.set(xlim=(0.5, 1.75), ylim=(-4.25, 3.0),
           xlabel=r"Mass ratio $x=M/M_{\rm crit}$", ylabel=r"Net mass rate $dx/ds$")
    ax.set_xticks([0.5, 0.75, 1, 1.25, 1.5, 1.75])
    fig.savefig(FIGURES / "figure1.pdf")
    plt.close(fig)


def trajectory(x0):
    """Integrate nonsingular transformed equations up to a visible cutoff."""
    below = x0 < 1.0
    cutoff = 0.05 if below else 2.2
    if below:
        y0 = [x0**3]
        y_stop = cutoff**3

        def rhs(_, y):
            return [3.0 * (max(y[0], 0.0)**(4.0 / 3.0) - 1.0)]

        to_x = lambda y: np.cbrt(y[0])
    else:
        y0 = [1.0 / x0]
        y_stop = 1.0 / cutoff

        def rhs(_, y):
            return [y[0]**4 - 1.0]

        to_x = lambda y: 1.0 / y[0]

    def event(_, y):
        return y[0] - y_stop

    event.terminal = True
    event.direction = -1
    sol = solve_ivp(rhs, (0, 2), y0, method="DOP853", events=event,
                    rtol=1e-10, atol=1e-12, max_step=0.02, dense_output=True)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError(f"Integration failed at x0={x0}: {sol.message}")
    end = float(sol.t_events[0][0])
    times = np.linspace(0.0, end, 501)
    masses = to_x(sol.sol(times))
    residual = float(np.max(np.abs(times - (primitive(masses) - primitive(x0)))))
    event_error = abs(end - float(primitive(cutoff) - primitive(x0)))
    if residual > 1e-8 or event_error > 1e-8:
        raise AssertionError(f"Analytic check failed at x0={x0}")
    return times, masses, residual, event_error


def figure_two():
    fig, ax = plt.subplots(figsize=(6.5, 3.65))
    checks = {}
    for x0, color, style in [(.80, "#c74436", "-"), (.95, "#c74436", "--"),
                             (1.05, "#1675ae", "--"), (1.20, "#1675ae", "-")]:
        times, masses, residual, event_error = trajectory(x0)
        ax.plot(times, masses, color=color, ls=style, lw=2, label=f"$x_0={x0:.2f}$")
        ax.plot([times[-1]], [masses[-1]], "o", color=color, markersize=4.5)
        checks[f"{x0:.2f}"] = {"event_time_s": float(times[-1]),
                              "maximum_implicit_residual": residual,
                              "event_time_error": event_error}
    ax.axhline(1, color="#71808a", linestyle=":", lw=1.3)
    ax.text(.41, 1.05, r"$x_0=1$: constant mass", color="#556571", fontsize=10.5)
    ax.set(xlim=(0, .9), ylim=(0, 2.34), xlabel=r"Time ratio $s=t/t_*$",
           ylabel=r"Mass ratio $x=M/M_{\rm crit}$")
    ax.legend(ncol=4, frameon=False, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), borderaxespad=0.0,
              columnspacing=0.9, handlelength=1.6)
    fig.savefig(FIGURES / "figure2.pdf")
    plt.close(fig)
    return checks


def figure_three():
    q = np.logspace(-1, 16, 700)
    m_ratio = MC / M0 * q**-0.25
    q_balance = (MC / M0)**4
    fig, ax = plt.subplots(figsize=(6.5, 3.1))
    ax.loglog(q, m_ratio, color="#263846", lw=2.3)
    ax.axhline(1, ls="--", lw=1, color="#71808a")
    ax.scatter([1, q_balance], [MC / M0, 1], c=["#c74436", "#1675ae"], zorder=3)
    ax.annotate("Reference parameters", xy=(1, MC/M0), xytext=(20, 7600),
                arrowprops={"arrowstyle": "-", "color": "#71808a"}, color="#c74436")
    ax.annotate(r"Balance at $M_0$", xy=(q_balance, 1), xytext=(3e7, 3.5),
                arrowprops={"arrowstyle": "-", "color": "#71808a"}, color="#1675ae")
    ax.set(xlim=(.1, 1e16), ylim=(.35, 1e4),
           xlabel=r"Accretion coefficient ratio $Q=A/A_{\rm ref}$",
           ylabel=r"Critical mass ratio $M_{\rm crit}/M_0$")
    ax.set_xticks([1, 1e4, 1e8, 1e12, 1e16])
    ax.grid(axis="y", alpha=.3)
    fig.savefig(FIGURES / "figure3.pdf")
    plt.close(fig)
    return q_balance


def main():
    FIGURES.mkdir(exist_ok=True)
    plot_style()
    figure_one()
    checks = figure_two()
    q_balance = figure_three()
    values = {
        "A_kg^-1_s^-1": A, "B_kg^3_s^-1": B, "K_W_kg^2": K,
        "Mcrit_kg": MC, "M0_over_Mcrit": M0 / MC,
        "tstar_years": TSTAR / YEAR, "efold_years": TSTAR / (4*YEAR),
        "tau_bb_hours": TAU_BB / 3600.0, "Q_balance": q_balance,
        "blackbody_power_M0_W": K / M0**2,
        "hawking_temperature_M0_K": HBAR*C**3/(8*math.pi*G*M0*KB),
        "capture_radius_M0_m": G*M0/CS**2,
        "evap_time_x0_0.95_tstar": -float(primitive(.95)),
        "double_time_x0_1.05_tstar": float(primitive(2.1)-primitive(1.05)),
        "numerical_checks": checks,
        "software": {"numpy": np.__version__, "scipy": __import__("scipy").__version__,
                     "matplotlib": matplotlib.__version__},
    }
    (OUT / "results.json").write_text(json.dumps(values, indent=2) + "\n")
    print(json.dumps(values, indent=2))


if __name__ == "__main__":
    main()
