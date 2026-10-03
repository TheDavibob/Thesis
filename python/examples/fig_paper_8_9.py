"""Corrected paper Figs 8-9: wall pressure and far-field noise for the base and
sheared boundary layers of paper Fig. 5 (inspired by Afshari et al. 2016).

A unit vortex sheet at y2 in {0.1, ..., 0.5} delta, hard plate, c0 = 35 U_inf.
Frequency f = 10 omega / (2 pi 0.02) Hz (U_inf = 10 m/s, delta = 20 mm), so the
paper's 10^2..10^4 Hz is omega ~ 1.3..126.

* p_w: wall pressure of the regularised Green's function, eq. (6.4.4).
* p_ff: far field directly above the plate (theta = pi/2): p_w times the
  physical transfer amplitude lim r^(1/2)|p'| per unit wall pressure, with
  kappa1 = omega / U(y2) (corrected factorisation: winding cancelled, tanh
  contours above all wake modes).

Two routes:

* "fit" (primary): degree-8 analytic fits to the two-piece profiles
  (``profiles.fit_analytic``; max |dU| 0.0012 / 0.0032).
* "two-piece" (cross-check at a subset of frequencies): the original piecewise
  profiles, integrated on real-axis paths lifted locally over the critical
  point.  Differences between the routes reflect the profiles' different U''
  at the critical level (critical-layer physics), not numerical error.

Usage::

    uv run --extra plot python examples/fig_paper_8_9.py [--fast] [--replot] [--tikz]

Plots stop at ``F_MAX`` = 5 kHz: above it the levels are 150-250 dB below the
low-frequency ones and the sheared/base ratios are not meaningful.
``--tikz`` also writes one matlab2tikz-style ``.tex`` per panel.
"""
import logging
import os
import sys
import time

import numpy as np

fast = "--fast" in sys.argv
stem = "fig_paper_8_9" + ("_fast" if fast else "")
Y2 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
C0 = 35.0
F_HZ = np.logspace(2, 4, 5 if fast else 14)
OMEGA = F_HZ * 2 * np.pi * 0.02 / 10
TWO_PIECE_IDX = [0, 2, 4] if fast else [0, 4, 8, 13]
ROUTES = ("fit", "two-piece")
F_MAX = 5000.0   # plotted range (Hz)
CASES = ("base", "sheared")


def n_bl(omega, U0):
    """Boundary-layer RK steps: resolve |k1| ~ omega/U0 (the original used 3000)."""
    return int(max(1500 if fast else 3000, 8 * omega / U0))


def transfer_at_90(ff, omega, U):
    """Physical far-field amplitude per unit wall pressure at theta = pi/2."""
    from shear_flow_scattering.steepest_descent import far_field_saddle

    sd = far_field_saddle([np.pi / 2], omega, U.f(1.0), C0)
    # (cos^2 + beta^2 sin^2)^(-1/4) = beta^(-1/2) at theta = pi/2
    return np.abs(np.ravel(ff.FF1p.p)) * np.sqrt(sd.K0 / (2 * np.pi)) / np.sqrt(sd.beta)


def compute():
    from shear_flow_scattering.profiles import afshari_profiles, c_const, fit_analytic
    from shear_flow_scattering.scripts.common import wall_operators
    from shear_flow_scattering.scripts.trailing_edge import (far_field_trailing_edge_lined,
                                                             vortex_wall_pressure)

    c = c_const(C0)
    base, sheared = afshari_profiles()
    profs = {("fit", "base"): fit_analytic(base), ("fit", "sheared"): fit_analytic(sheared),
             ("two-piece", "base"): base, ("two-piece", "sheared"): sheared}
    res = {}
    for (route, name), U in profs.items():
        idx = range(OMEGA.size) if route == "fit" else TWO_PIECE_IDX
        pw = np.full((OMEGA.size, Y2.size), np.nan + 0j)
        T = np.full((OMEGA.size, Y2.size), np.nan)
        for i in idx:
            w, t0 = OMEGA[i], time.time()
            N = n_bl(w, float(np.real(U.f(0))))
            L = wall_operators(w, U).hard
            pw[i] = vortex_wall_pressure(Y2, w, 0.0, U, c, 1.0, N, L, 1.0, "regularised")
            ff = far_field_trailing_edge_lined(Up=U, Um=U, cp=c, cm=c, omega=w, Z=np.inf,
                                               theta1=[np.pi / 2], theta2=[-np.pi / 2],
                                               y2=Y2, kappa1=w / U.f(Y2), pressure_upstream=1.0,
                                               N_bl=N, prec_N=1500 if fast else 3000)
            T[i] = transfer_at_90(ff, w, U)
            logging.info("%s %s f=%.0f Hz (omega=%.2f, N_bl=%d): %.0fs",
                         route, name, F_HZ[i], w, N, time.time() - t0)
        res[(route, name)] = (pw, T)
        save(res)  # checkpoint after each profile
    return res


def save(res):
    keys = list(res)
    np.savez(stem + ".npz", f=F_HZ, omega=OMEGA, y2=Y2,
             route=np.array([k[0] for k in keys]), case=np.array([k[1] for k in keys]),
             **{f"pw_{i}": res[k][0] for i, k in enumerate(keys)},
             **{f"T_{i}": res[k][1] for i, k in enumerate(keys)})


def load():
    d = np.load(stem + ".npz")
    return {(str(r), str(n)): (d[f"pw_{i}"], d[f"T_{i}"])
            for i, (r, n) in enumerate(zip(d["route"], d["case"]))}


def db(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20 * np.log10(np.abs(x))


STYLE = {"fit": dict(ls="-", lw=1.8), "two-piece": dict(ls="none", marker="o", ms=5)}

#: Two-piece points flagged unreliable by the solver (see fig89_full.log):
#: a source at the profile join puts the critical point on the U'' jump, and
#: for the sheared case at higher frequency the kernel's critical points could
#: not be passed within one analytic piece.  (fit route: no warnings.)
JOIN = {"base": 0.1, "sheared": 0.3}
FF_FLAGGED = {"sheared": [412, 1701, 10000]}   # Hz, far field only


def masked(route, name, what, val):
    """Apply the reliability mask to a two-piece series."""
    if route != "two-piece":
        return val
    val = val.copy()
    for f in FF_FLAGGED.get(name, []) if what == "pff" else []:
        val[np.isclose(F_HZ, f, rtol=0.01)] = np.nan
    return val


def tikz(fig, stem_, names):
    if "--tikz" in sys.argv:
        from shear_flow_scattering.tikz import save_figure_tikz
        print("wrote", *save_figure_tikz(fig, stem_, names))


def plot(res):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True, sharex=True)
    for col, name in enumerate(CASES):
        for row, what in enumerate(("pw", "pff")):
            ax = axes[row, col]
            for j, y in enumerate(Y2):
                for route in ROUTES:
                    if (route, name) not in res:
                        continue
                    pw, T = res[(route, name)]
                    if route == "two-piece" and np.isclose(y, JOIN[name]):
                        continue
                    val = masked(route, name, what, pw[:, j] if what == "pw" else pw[:, j] * T[:, j])
                    ax.plot(F_HZ, db(val), color=f"C{j}",
                            label=(f"$y_2 = {y}\\delta$" if route == "fit" else None), **STYLE[route])
            ax.set_xscale("log")
            ax.set_xlim(F_HZ[0], F_MAX)
            ax.set_ylim(-250, 10)
            label = "Wall" if what == "pw" else "Far-field"
            ax.set_title(f"({'abcd'[2 * row + col]}) {label} pressure, {name} case")
            ax.set_ylabel(r"$20\log_{10}|p_w|$" if what == "pw" else r"$20\log_{10}|p_{ff}|$")
            if row == 1:
                ax.set_xlabel("Frequency, f (Hz)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Corrected Fig. 8 -- lines: analytic fits; circles: two-piece profiles (indented paths; "
                 "points flagged unreliable omitted)", fontsize=10)
    fig.savefig(stem + "_fig8.png", dpi=130)
    tikz(fig, stem + "_fig8", ["WallPressure_base", "WallPressure_sheared", "FFPressure_base", "FFPressure_sheared"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True, sharey=True)
    for ax, what in zip(axes, ("pw", "pff")):
        for j, y in enumerate(Y2):
            for route in ROUTES:
                if (route, "base") not in res or (route, "sheared") not in res:
                    continue
                if route == "two-piece" and any(np.isclose(y, v) for v in JOIN.values()):
                    continue
                pb, Tb = res[(route, "base")]
                ps, Ts = res[(route, "sheared")]
                r = ps[:, j] / pb[:, j] if what == "pw" else ps[:, j] * Ts[:, j] / (pb[:, j] * Tb[:, j])
                r = masked(route, "sheared", what, r)
                ax.plot(F_HZ, db(r), color=f"C{j}",
                        label=(f"$y_2 = {y}\\delta$" if route == "fit" else None), **STYLE[route])
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_xlim(F_HZ[0], F_MAX)
        ax.set_ylim(-40, 20)
        ax.set_xlabel("Frequency, f (Hz)")
        ax.set_title("(a) wall pressure, sheared / base" if what == "pw" else "(b) far field, sheared / base")
    axes[0].set_ylabel("dB")
    axes[0].legend(fontsize=8)
    fig.suptitle("Corrected Fig. 9 -- lines: analytic fits; circles: two-piece profiles", fontsize=10)
    fig.savefig(stem + "_fig9.png", dpi=130)
    tikz(fig, stem + "_fig9", ["WallPressure_reduction", "FFPressure_reduction"])
    print("saved", stem + "_fig8.png", stem + "_fig9.png")


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    plot(load() if "--replot" in sys.argv else compute())
