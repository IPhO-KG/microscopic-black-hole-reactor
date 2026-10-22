#!/usr/bin/env python3
"""Reproduce the numbers, tables and figures of the thermal-feedback paper.

Run from this directory: python3 reproduce.py
Requires numpy, scipy and matplotlib. Writes results.json and four vector
figures in figures/. Figure 1 (the feedback loop) is drawn in main.tex.

Dimensionless model (Eq. 4 of the paper):
    dx/ds     = theta**(-3/2) * x**2 - x**(-2)
    dtheta/ds = eta * x**(-2) - zeta * (theta - 1),   eta = chi * zeta
"""

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent
FIGURES = OUT / "figures"
SOLVER = {"method": "DOP853", "rtol": 1e-10, "atol": 1e-12, "max_step": 0.01}

BLUE, ORANGE, GRAY, INK = "#1f5a99", "#d0632a", "#71808a", "#263846"

# SI constants (CODATA 2018 exact or recommended values)
G = 6.67430e-11
C = 299792458.0
HBAR = 6.62607015e-34 / (2.0 * math.pi)
KB = 1.380649e-23
EV = 1.602176634e-19
MP = 1.67262192369e-27
SIGMA_T = 6.6524587321e-29
K_BB = HBAR * C**6 / (15360.0 * math.pi * G**2)   # P_bb = K / M^2
B_EV = K_BB / C**2                                 # dM/dt = -B / M^2


# ---------------------------------------------------------------- model
def rhs(_s, y, chi, zeta):
    x, th = y
    return [th**-1.5 * x**2 - x**-2, chi * zeta * x**-2 - zeta * (th - 1.0)]


def equilibrium(chi):
    """Unique positive root of theta - 1 - chi*theta**(-3/4) = 0."""
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


# ------------------------------------------------- physical diagnostics
def plasma_row(T0, n):
    """Order-of-magnitude scales for fully ionized hydrogen (Sec. 6, Table 2).

    gamma = 5/3, mu = 1/2, lambda = 1/4. Free-free cooling and the ion
    collision rate follow the NRL Plasma Formulary (cgs/eV inputs converted).
    """
    cs = math.sqrt(5.0 / 3.0 * KB * T0 / (0.5 * MP))
    rho = n * MP
    A0 = 4.0 * math.pi * 0.25 * rho * G**2 / cs**3
    Mc0 = (B_EV / A0) ** 0.25
    rB = G * Mc0 / cs**2
    t_flow = rB / cs
    tau_T = n * SIGMA_T * rB
    P_bb = K_BB / Mc0**2
    E_B = 3.0 * n * KB * T0 * (4.0 * math.pi * rB**3 / 3.0)
    t_heat = E_B / (tau_T * P_bb)
    T_eV = KB * T0 / EV
    q_ff = 1.69e-38 * n**2 * math.sqrt(T_eV)          # W m^-3, n in m^-3
    t_ff = 3.0 * n * KB * T0 / q_ff
    nu_ii = 4.80e-8 * (n * 1e-6) * T_eV**-1.5 * 10.0   # s^-1, ln(Lambda) = 10
    v_ti = 9.79e3 * math.sqrt(T_eV)                    # m s^-1, (k T / m_p)^1/2
    return {"T0_K": T0, "n_m3": n, "cs_m_s": cs, "Mc0_kg": Mc0, "rB_m": rB,
            "t_flow_s": t_flow, "tau_T": tau_T, "tauT_Pbb_W": tau_T * P_bb,
            "t_heat_T_s": t_heat, "t_ff_s": t_ff, "mfp_ii_over_rB": v_ti / nu_ii / rB,
            "t_heat_over_t_flow": t_heat / t_flow, "t_ff_over_t_flow": t_ff / t_flow,
            "chi_if_tau_c_eq_t_flow": t_flow / t_heat}


# -------------------------------------------------------------- figures
def plot_style():
    plt.rcParams.update({
        "font.size": 10.5, "axes.labelsize": 11.0, "legend.fontsize": 9.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
    })


def figure_equilibrium():
    chis = np.logspace(-3, 3, 241)
    eq = np.array([equilibrium(c) for c in chis])
    r = [rates(c, 1.0) for c in chis]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.3))
    ax1.loglog(chis, eq[:, 1], color=ORANGE, lw=2, label=r"$\theta_*$")
    ax1.loglog(chis, eq[:, 0], color=BLUE, lw=2, label=r"$x_*$")
    ax1.set(xlabel=r"Feedback strength $\chi$", ylabel="Equilibrium value")
    ax1.legend(frameon=False)
    ax2.loglog(chis, [q["lambda_plus"] for q in r], color=ORANGE, lw=2,
               label=r"$\lambda_+$ (temperature responds)")
    ax2.loglog(chis, [q["a"] for q in r], color=BLUE, lw=1.8, ls="--",
               label=r"$a$ (temperature frozen at $\theta_*$)")
    ax2.set(xlabel=r"Feedback strength $\chi$", ylabel=r"Growth rate ($\zeta=1$)")
    ax2.legend(frameon=False, loc="lower left")
    for ax, tag in ((ax1, "(a)"), (ax2, "(b)")):
        ax.grid(which="major", alpha=0.18)
        ax.text(-0.16, 1.02, tag, transform=ax.transAxes, fontweight="bold")
    fig.tight_layout(w_pad=2.5)
    fig.savefig(FIGURES / "equilibrium.pdf")
    plt.close(fig)


def _arrow_on(ax, xs, ys, dist, color, center):
    """Arrowhead along (xs, ys), in forward time, at distance `dist` from `center`."""
    d = np.hypot(xs - center[0], ys - center[1])
    i = int(np.argmin(np.abs(d - dist)))
    i = min(max(i, 1), len(xs) - 2)
    ax.annotate("", xy=(xs[i + 1], ys[i + 1]), xytext=(xs[i - 1], ys[i - 1]),
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.6,
                            "mutation_scale": 14})


def figure_phase(chi=1.0, zeta=1.0):
    J, xs, ts = jacobian(chi, zeta)
    w, v = np.linalg.eig(J)
    order = np.argsort(w.real)
    stable, unstable = v[:, order[0]].real, v[:, order[1]].real
    xlim, ylim = (0.75, 1.95), (0.85, 2.5)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    X, T = np.meshgrid(np.linspace(*xlim, 19), np.linspace(*ylim, 15))
    U, V = rhs(0, [X, T], chi, zeta)
    n = np.hypot(U, V)
    ax.quiver(X, T, U / n, V / n, color="#c9cdd1", width=0.0025, scale=42)
    th = np.linspace(*ylim, 300)
    ax.plot(th**0.375, th, color=GRAY, ls="--", lw=1.4, label=r"Mass nullcline $\dot x=0$")
    x = np.linspace(*xlim, 300)
    ax.plot(x, 1 + chi * x**-2, color=GRAY, ls=":", lw=1.6,
            label=r"Thermal nullcline $\dot\theta=0$")

    def stop(_s, y, *_):
        return min(y[0] - xlim[0], xlim[1] - y[0], y[1] - ylim[0], ylim[1] - y[1])
    stop.terminal = True
    for sign in (1, -1):
        sol = solve_ivp(rhs, (0, -40), np.array([xs, ts]) + sign * 1e-6 * stable,
                        args=(chi, zeta), events=stop, **SOLVER)
        xr, yr = sol.y[0][::-1], sol.y[1][::-1]           # forward-time order
        ax.plot(xr, yr, color=BLUE, lw=2, label="Stable manifold" if sign == 1 else None)
        _arrow_on(ax, xr, yr, 0.38, BLUE, (xs, ts))
        sol = solve_ivp(rhs, (0, 40), np.array([xs, ts]) + sign * 1e-6 * unstable,
                        args=(chi, zeta), events=stop, **SOLVER)
        ax.plot(sol.y[0], sol.y[1], color=ORANGE, lw=2,
                label="Unstable manifold" if sign == 1 else None)
        _arrow_on(ax, sol.y[0], sol.y[1], 0.38, ORANGE, (xs, ts))
    ax.plot([xs], [ts], "o", ms=7, mfc="white", mec=INK, mew=1.5, zorder=5)
    ax.text(1.62, 2.30, r"$\dot x>0$", color=INK)
    ax.text(0.80, 1.02, r"$\dot x<0$", color=INK)
    ax.set(xlim=xlim, ylim=ylim, xlabel=r"Mass ratio $x=M/M_{c0}$",
           ylabel=r"Temperature ratio $\theta=T/T_0$")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    fig.savefig(FIGURES / "phase_portrait.pdf")
    plt.close(fig)


def figure_scan(n=201):
    chis = np.logspace(-2, 2, n)
    zetas = np.logspace(-2, 2, n)
    lam = np.array([[rates(c, z)["lambda_plus"] for c in chis] for z in zetas])
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    mesh = ax.pcolormesh(chis, zetas, lam, shading="gouraud", cmap="viridis",
                         norm=LogNorm(vmin=lam.min(), vmax=lam.max()), rasterized=True)
    cs = ax.contour(chis, zetas, lam, levels=[0.3, 0.5, 1, 2, 3], colors="white",
                    linewidths=0.9)
    for label in ax.clabel(cs, fmt="%g", fontsize=8.5, inline_spacing=6,
                           manual=[(60, 1), (25, 1), (8, 1), (1.5, 1), (0.4, 1)]):
        label.set_rotation(0)
    ax.set(xscale="log", yscale="log", xlabel=r"Feedback strength $\chi$",
           ylabel=r"Cooling rate $\zeta=t_*/\tau_c$")
    ax.spines[["top", "right"]].set_visible(True)
    cb = fig.colorbar(mesh, ax=ax, ticks=[0.3, 0.5, 1, 2, 3])
    cb.ax.yaxis.set_minor_formatter(plt.NullFormatter())
    cb.ax.set_yticklabels(["0.3", "0.5", "1", "2", "3"])
    cb.set_label(r"Unstable eigenvalue $\lambda_+$")
    fig.savefig(FIGURES / "parameter_scan.pdf", dpi=300)
    plt.close(fig)
    return float(lam.min()), float(lam.max())


def figure_departures(chi=1.0, zeta=1.0, eps=0.02):
    xs, ts = equilibrium(chi)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))

    def run(f, y0, scale):
        lo = lambda s, y: y[0] / scale - 0.65
        hi = lambda s, y: y[0] / scale - 1.65
        lo.terminal = hi.terminal = True
        sol = solve_ivp(f, (0, 10), y0, events=[lo, hi], **SOLVER)
        return sol.t, sol.y[0] / scale

    fixed = lambda s, y: [y[0]**2 - y[0]**-2]
    coupled = lambda s, y: rhs(s, y, chi, zeta)
    for sign, ls in ((1, "-"), (-1, "--")):
        lab = f"{sign * 2:+d}%"
        t, m = run(fixed, [1 + sign * eps], 1.0)
        ax.plot(t, m, color=BLUE, ls=ls, lw=2, label=f"Fixed medium, {lab}")
        t, m = run(coupled, [xs * (1 + sign * eps), ts], xs)
        ax.plot(t, m, color=ORANGE, ls=ls, lw=2, label=f"Thermal feedback, {lab}")
    ax.axhline(1, color=GRAY, lw=0.8, ls=":")
    ax.set(xlabel=r"Dimensionless time $s=t/t_*$", ylabel="Mass / equilibrium mass")
    ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    fig.savefig(FIGURES / "departures.pdf")
    plt.close(fig)


# ----------------------------------------------------------------- main
def main():
    FIGURES.mkdir(exist_ok=True)
    plot_style()
    table1 = [rates(c, 1.0) for c in (0.01, 0.1, 1.0, 10.0, 100.0)]
    grid = [rates(c, z) for c in np.logspace(-2, 2, 41) for z in np.logspace(-2, 2, 41)]
    table2 = {"hot_plasma": plasma_row(1.0e6, 1.0e20),
              "core_proxy": plasma_row(1.56e7, 1.5e5 / MP)}
    figure_equilibrium()
    figure_phase()
    lam_min, lam_max = figure_scan()
    figure_departures()
    checks = {
        "grid_points": len(grid),
        "max_equilibrium_residual": max(
            abs(q["theta_star"] - 1 - q["chi"] * q["theta_star"]**-0.75) for q in grid),
        "all_det_negative": all(q["det"] < 0 for q in grid),
        "all_lambda_plus_above_a": all(q["lambda_plus"] > q["a"] > 0 for q in grid),
        "all_lambda_minus_below_minus_zeta": all(q["lambda_minus"] < -q["zeta"] for q in grid),
        "fine_scan_lambda_plus_range": [lam_min, lam_max],
    }
    results = {"reference_point": rates(1.0, 1.0), "table1_zeta_1": table1,
               "table2_plasma_diagnostics": table2, "checks": checks}
    (OUT / "results.json").write_text(json.dumps(results, indent=2))

    print(f"{'chi':>6} {'theta*':>10} {'x*':>9} {'a':>9} {'lam-':>10} {'lam+':>9} {'det':>9}")
    for q in table1:
        print(f"{q['chi']:6g} {q['theta_star']:10.5f} {q['x_star']:9.5f} {q['a']:9.5f} "
              f"{q['lambda_minus']:10.5f} {q['lambda_plus']:9.5f} {q['det']:9.5f}")
    for name, row in table2.items():
        print(name, {k: f"{v:.3g}" for k, v in row.items()})
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
