#!/usr/bin/env python3
"""Reproduce the numbers and figures of Paper II (coupled Hawking-Bondi thermal feedback).

Run from this directory: python3 reproduce.py
Requires numpy, scipy and matplotlib. Writes results.json and figures/figure1-3.pdf.

Dimensionless model (Eq. 2 of the paper):
    dx/ds     = theta**(-3/2) * x**2 - x**(-2)
    dtheta/ds = eta * x**(-2) - zeta * (theta - 1),   eta = chi * zeta
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent
FIGURES = OUT / "figures"
SOLVER = {"method": "DOP853", "rtol": 1e-10, "atol": 1e-12, "max_step": 0.01}


def rhs(_s, y, chi, zeta):
    x, th = y
    return [th**-1.5 * x**2 - x**-2, chi * zeta * x**-2 - zeta * (th - 1.0)]


def equilibrium(chi):
    """Unique positive root of theta - 1 - chi*theta**(-3/4) = 0 (Eq. 3)."""
    h = lambda th: th - 1.0 - chi * th**-0.75
    th = brentq(h, 1.0, 2.0 + chi, xtol=1e-15, rtol=1e-15)
    return th**0.375, th


def jacobian(chi, zeta):
    x, th = equilibrium(chi)
    eta = chi * zeta
    a = 4.0 * x**-3
    b = -1.5 * th**-2.5 * x**2
    c = -2.0 * eta * x**-3
    return np.array([[a, b], [c, -zeta]]), x, th


def rates(chi, zeta):
    J, x, th = jacobian(chi, zeta)
    lam = np.sort(np.linalg.eigvals(J).real)
    return {"chi": chi, "zeta": zeta, "x_star": x, "theta_star": th,
            "a": J[0, 0], "lambda_minus": lam[0], "lambda_plus": lam[1],
            "det": float(np.linalg.det(J))}


def plot_style():
    plt.rcParams.update({
        "font.size": 11.0, "axes.labelsize": 11.5, "legend.fontsize": 9.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
    })


def figure_one():
    chis = np.logspace(-3, 3, 241)
    eq = np.array([equilibrium(c) for c in chis])
    r = [rates(c, 1.0) for c in chis]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.4))
    ax1.loglog(chis, eq[:, 1], color="#d0632a", lw=2, label=r"$\theta_*$")
    ax1.loglog(chis, eq[:, 0], color="#1f5a99", lw=2, label=r"$x_*$")
    ax1.set(xlabel=r"Feedback strength $\chi$", ylabel="Equilibrium ratio")
    ax1.legend(frameon=False)
    ax2.loglog(chis, [q["lambda_plus"] for q in r], color="#d0632a", lw=2,
               label=r"Coupled $\lambda_+$")
    ax2.loglog(chis, [q["a"] for q in r], color="#1f5a99", lw=1.8, ls="--",
               label=r"Frozen $a$")
    ax2.set(xlabel=r"Feedback strength $\chi$", ylabel="Dimensionless growth rate")
    ax2.legend(frameon=False)
    for ax in (ax1, ax2):
        ax.grid(which="both", alpha=0.15)
    fig.savefig(FIGURES / "figure1.pdf")
    plt.close(fig)


def figure_two(chi=1.0, zeta=1.0):
    J, xs, ts = jacobian(chi, zeta)
    w, v = np.linalg.eig(J)
    order = np.argsort(w.real)
    stable, unstable = v[:, order[0]].real, v[:, order[1]].real
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    X, T = np.meshgrid(np.linspace(0.75, 1.95, 18), np.linspace(0.9, 2.5, 15))
    U, V = rhs(0, [X, T], chi, zeta)
    n = np.hypot(U, V)
    ax.quiver(X, T, U / n, V / n, color="#c9cdd1", width=0.0025, scale=40)
    th = np.linspace(0.85, 2.5, 300)
    ax.plot(th**0.375, th, color="#71808a", ls="--", label=r"Mass nullcline $F=0$")
    x = np.linspace(0.75, 1.95, 300)
    ax.plot(x, 1 + chi * x**-2, color="#71808a", ls=":", label=r"Thermal nullcline $G=0$")
    stop = lambda s, y, *_: min(y[0] - 0.75, 1.95 - y[0], y[1] - 0.85, 2.5 - y[1])
    stop.terminal = True
    for sign in (1, -1):
        y0 = np.array([xs, ts]) + sign * 1e-6 * stable
        sol = solve_ivp(rhs, (0, -40), y0, args=(chi, zeta), events=stop, **SOLVER)
        ax.plot(sol.y[0], sol.y[1], color="#1f5a99", lw=2,
                label="Stable manifold" if sign == 1 else None)
        y0 = np.array([xs, ts]) + sign * 1e-6 * unstable
        sol = solve_ivp(rhs, (0, 40), y0, args=(chi, zeta), events=stop, **SOLVER)
        ax.plot(sol.y[0], sol.y[1], color="#d0632a", lw=2,
                label="Unstable manifold" if sign == 1 else None)
    ax.plot([xs], [ts], "o", ms=7, mfc="white", mec="#263846", mew=1.5, zorder=5)
    ax.set(xlim=(0.75, 1.95), ylim=(0.85, 2.5), xlabel=r"Mass ratio $x=M/M_{c0}$",
           ylabel=r"Temperature ratio $\theta=T/T_0$")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.savefig(FIGURES / "figure2.pdf")
    plt.close(fig)


def figure_three(chi=1.0, zeta=1.0, eps=0.02):
    xs, ts = equilibrium(chi)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))

    def run(f, y0, scale):
        lo = lambda s, y, *_: y[0] / scale - 0.65
        hi = lambda s, y, *_: y[0] / scale - 1.65
        lo.terminal = hi.terminal = True
        sol = solve_ivp(f, (0, 10), y0, events=[lo, hi], **SOLVER)
        return sol.t, sol.y[0] / scale

    fixed = lambda s, y: [y[0]**2 - y[0]**-2]
    coupled = lambda s, y: rhs(s, y, chi, zeta)
    for sign, ls in ((1, "-"), (-1, "--")):
        lab = f"{sign * 2:+d}%"
        t, m = run(fixed, [1 + sign * eps], 1.0)
        ax.plot(t, m, color="#1f5a99", ls=ls, lw=2, label=f"Fixed medium, {lab}")
        t, m = run(coupled, [xs * (1 + sign * eps), ts], xs)
        ax.plot(t, m, color="#d0632a", ls=ls, lw=2, label=f"Thermal feedback, {lab}")
    ax.axhline(1, color="#71808a", lw=0.8, ls=":")
    ax.set(xlabel=r"Dimensionless time $s=t/t_*$", ylabel="Mass / its equilibrium mass")
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    fig.savefig(FIGURES / "figure3.pdf")
    plt.close(fig)


def main():
    FIGURES.mkdir(exist_ok=True)
    plot_style()
    table = [rates(c, 1.0) for c in (0.01, 0.1, 1.0, 10.0, 100.0)]
    grid = [rates(c, z) for c in np.logspace(-2, 2, 41) for z in np.logspace(-2, 2, 41)]
    checks = {
        "max_equilibrium_residual": max(
            abs(q["theta_star"] - 1 - q["chi"] * q["theta_star"]**-0.75) for q in grid),
        "all_det_negative": all(q["det"] < 0 for q in grid),
        "all_lambda_plus_positive": all(q["lambda_plus"] > 0 for q in grid),
        "all_lambda_plus_above_a": all(q["lambda_plus"] > q["a"] for q in grid),
        "grid_points": len(grid),
    }
    figure_one()
    figure_two()
    figure_three()
    results = {"table1_zeta_1": table, "reference_point": rates(1.0, 1.0), "checks": checks}
    (OUT / "results.json").write_text(json.dumps(results, indent=2))
    print(f"{'chi':>6} {'theta*':>10} {'x*':>9} {'a':>9} {'lam-':>10} {'lam+':>9}")
    for q in table:
        print(f"{q['chi']:6g} {q['theta_star']:10.5f} {q['x_star']:9.5f} {q['a']:9.5f} "
              f"{q['lambda_minus']:10.5f} {q['lambda_plus']:9.5f}")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
