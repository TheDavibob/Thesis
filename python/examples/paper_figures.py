"""All computed figures of the JFM trailing-edge noise draft, recomputed with the
corrected methods only, as PNG and as matlab2tikz-style TikZ (one file per panel).

Outputs go to ``examples/figures/``, named ``figNN<panel>_<description>`` after
the current paper figure numbers (``figNew42`` is the proposed section 4.2
figure).  Figs 1, 2, 12 and 13 are schematics and are not generated.

=====  ==========================================================================
Fig    What, and what changed relative to the draft
=====  ==========================================================================
3      Vortex sheet on a hard trailing edge, Re p' total (thesis Fig. 6.2.2).
       Regularised incident field, so the amplitude is that of a unit sheet.
4      Distributed vorticity omega_c = 1, Re p' total.  Corrected weighting
       (``weighting="jacobian"``, the literal y2 -> kappa1 change of variables).
5      The two-piece base/sheared profiles and the degree-8 analytic fits used
       for Figs 8-9.
6      Wall-pressure-to-far-field transfer, corrected only (winding cancelled,
       tanh contours above all wake modes); physical amplitude lim r^1/2 |p'|.
7      Wall pressure of a unit vortex sheet over the parabolic profile U0 = 0.3,
       Ma = 0.2, on the causal branch.  Points more than 150 dB down at the
       previous frequency are not computed (they are far below the colour scale).
8, 9   Corrected wall and far-field pressure (analytic fits; two-piece cross-
       checks as circles), truncated at 5 kHz where the results are converged.
10     How the semicircular path C+- moves the critical layer, U0 = 0.2, omega = 1.
11     Phase of D_h = v_d(0) via the raw semicircle C+ (no path safety step), with
       a perceptually uniform cyclic colour map (CET-C6) and log10|D_h| contours.
New42  Section 4.2 parameter studies: asymmetry, Mach number, profile shape.
=====  ==========================================================================

Usage (from ``examples/``)::

    uv run --extra plot python paper_figures.py [fig03 fig07 ...] [--replot] [--workers N]

With no figure names, all are made.  Results are cached in
``figures/cache/<name>.npz``; ``--replot`` replots from the cache.
"""
import logging
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
TIKZ = os.path.join(OUT, "tikz")
CACHE = os.path.join(OUT, "cache")
REPLOT = "--replot" in sys.argv
WORKERS = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def cached(name, compute):
    path = os.path.join(CACHE, name + ".npz")
    if REPLOT:
        d = np.load(path, allow_pickle=False)
        return {k: d[k] for k in d.files}
    t0 = time.time()
    data = compute()
    np.savez(path, **data)
    logging.info("%s computed in %.0f s", name, time.time() - t0)
    return data


def save(fig, stem, panels):
    """PNG of the whole figure plus one TikZ file per panel."""
    from shear_flow_scattering.tikz import save_tikz
    fig.savefig(os.path.join(OUT, stem + ".png"), dpi=150)
    num = stem.split("_")[0]
    axes = [a for a in fig.axes if not hasattr(a, "_colorbar")]
    assert len(axes) == len(panels), (stem, len(axes))
    letters = "abcdefgh" if len(axes) > 1 else [""]
    paths = [save_tikz(a, os.path.join(TIKZ, f"{num}{l}_{n}.tex")) for a, l, n in zip(axes, letters, panels)]
    logging.info("saved %s.png and %s", stem, ", ".join(os.path.basename(p) for p in paths))


def field_axes(ax, x1, x2, data, title, cmap="viridis"):
    lo, hi = np.percentile(data, [1, 99])   # clip the few extreme pixels near the critical layer
    m = ax.pcolormesh(x1, x2, data, shading="auto", cmap=cmap, vmin=lo, vmax=hi)
    ax.figure.colorbar(m, ax=ax, pad=0.02)
    ax.plot([x1.min(), 0], [0, 0], "k", lw=4)
    ax.set(xlabel="$x_1$", ylabel="$x_2$", title=title)


def db(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 20 * np.log10(np.abs(x))


# --------------------------------------------------------------------------
# Fig. 3 -- vortex sheet on a hard trailing edge
# --------------------------------------------------------------------------
def fig03():
    from shear_flow_scattering.profiles import c_const, u_parabolic
    from shear_flow_scattering.scripts.trailing_edge import trailing_edge_lined

    def compute():
        res = trailing_edge_lined(
            Up=u_parabolic(0.3, 1.0, 1.0), Um=u_parabolic(0.3, 1.5, 1.0),
            cp=c_const(5.0), cm=c_const(5.0), omega=1.0, Z=np.inf, y2=0.5, incident="regularised",
            x1=np.arange(-4, 10 + 1e-9, 0.05), x2=np.arange(-2, 2 + 1e-9, 0.02),
            N_bl=2000, N_C=500, N_CL=1000, prec_N=2000)
        return dict(x1=res.x1, x2=res.x2, p=res.I_total.p)

    d = cached("fig03", compute)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 3.2), constrained_layout=True)
    P = d["p"].real.T.copy()
    on = np.flatnonzero(np.isclose(d["x2"], 0.5))      # grid row on the sheet: one-sided limit,
    P[on] = 0.5 * (P[on - 1] + P[on + 1])               # p is continuous, so use the neighbours
    field_axes(ax, d["x1"], d["x2"], P, r"Re $p'_0$, unit vortex sheet at $y_2 = 0.5$")
    save(fig, "fig03_vortex_sheet_total_pressure", ["vortex_sheet_total_pressure"])


# --------------------------------------------------------------------------
# Fig. 4 -- distributed vorticity
# --------------------------------------------------------------------------
def fig04():
    from shear_flow_scattering.scripts.trailing_edge import distributed_trailing_edge

    def compute():
        res = distributed_trailing_edge(
            x1=np.arange(-5, 5 + 1e-9, 0.05), x2=np.arange(-2, 2 + 1e-9, 0.02),
            N_bl=2000, N_C=500, N_CL=1000, N_distribution=80, prec_N=2000, weighting="jacobian")
        return dict(x1=res.x1, x2=res.x2, p=res.I_total.p)

    d = cached("fig04", compute)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.2), constrained_layout=True)
    field_axes(ax, d["x1"], d["x2"], d["p"].real.T, r"Re $p'$, $\omega_c = 1$ on $0 < x_2 < 1$")
    save(fig, "fig04_distributed_vorticity_total_pressure", ["distributed_vorticity_total_pressure"])


# --------------------------------------------------------------------------
# Fig. 5 -- profiles
# --------------------------------------------------------------------------
def fig05():
    import matplotlib.pyplot as plt
    from shear_flow_scattering.profiles import afshari_profiles, fit_analytic

    x = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(4.5, 4), constrained_layout=True)
    for j, (name, U) in enumerate(zip(("Base", "Sheared"), afshari_profiles())):
        F = fit_analytic(U)
        ax.plot(np.real(U.f(x)), x, color=f"C{j}", lw=2, label=f"{name} profile")
        ax.plot(np.real(F.f(x)), x, color="k", lw=0.8, ls="--",
                label="degree-8 analytic fits" if j else None)
    ax.grid(True, color="0.85", lw=0.5)
    ax.axvline(1, color="0.35", lw=1.2)              # U = U_inf
    ax.axhline(1, color="0.35", lw=1.2)              # x2 = delta
    ax.set(xlabel=r"$U/U_\infty$", ylabel=r"$x_2/\delta$", xlim=(0, 1.1), ylim=(0, 1.1))
    ax.legend(loc="upper left")
    save(fig, "fig05_boundary_layer_profiles", ["boundary_layer_profiles"])


# --------------------------------------------------------------------------
# Fig. 6 -- wall-pressure-to-far-field transfer (corrected only)
# --------------------------------------------------------------------------
U0S_6, C0S_6 = (0.2, 0.6, 1.0), (50.0, 5.0)


def fig06():
    def compute():
        from sec42_studies import amplitude, theta_grids
        from shear_flow_scattering.profiles import c_const, u_parabolic
        from shear_flow_scattering.scripts.trailing_edge import far_field_trailing_edge_lined

        th1, th2 = theta_grids(181)
        out = {}
        for c0 in C0S_6:
            for U0 in U0S_6:
                logging.info("fig06 c0=%g U0=%g", c0, U0)
                U, c = u_parabolic(U0, 1, 1), c_const(c0)
                ff = far_field_trailing_edge_lined(Up=U, Um=U, cp=c, cm=c, omega=5.0, Z=np.inf, theta1=th1,
                                                   theta2=th2, pressure_upstream=1.0, kappa1=5.0 / 0.8,
                                                   factor_contour="tanh", winding="auto",
                                                   N_bl=1500, prec_N=4000)
                th, A = amplitude(th1, th2, ff.FF1p.p, ff.FF2p.p, U, U, c0)
                _, Am = amplitude(th1, th2, ff.FF1m.p, ff.FF2m.p, U, U, c0)
                down = (th < np.pi / 2) | (th > 3 * np.pi / 2)
                A[down] = Am[down]
                out[f"A_{c0:g}_{U0:g}"] = A
        out["theta"] = th
        return out

    d = cached("fig06", compute)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, subplot_kw={"projection": "polar"}, figsize=(10, 5), constrained_layout=True)
    for ax, c0, title in zip(axes, C0S_6, ("(a) Ma = 0.02, acoustically compact", "(b) Ma = 0.2")):
        for j, U0 in enumerate(U0S_6):
            ax.plot(d["theta"], d[f"A_{c0:g}_{U0:g}"], color=f"C{j}", lw=2,
                    label=(f"$U_0 = {U0:g}U_\\infty$" if U0 < 1 else "$U_0 = U_\\infty$"))
        ax.set_title(title, pad=12)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3, fontsize=8, frameon=False)
    save(fig, "fig06_wall_pressure_far_field_transfer",
         ["transfer_compact_Ma002", "transfer_noncompact_Ma02"])


# --------------------------------------------------------------------------
# Fig. 7 -- wall pressure vs (omega, y2), parabolic U0 = 0.3, Ma = 0.2
# --------------------------------------------------------------------------
OMEGA_7 = np.logspace(-1, 3, 25)
Y2_7 = np.linspace(0.05, 0.95, 19)
SKIP_DB = -150.0


def _fig07_row(args):
    """Wall pressure at one frequency for the source heights in ``ys``."""
    w, ys = args
    from shear_flow_scattering.profiles import c_const, u_parabolic
    from shear_flow_scattering.scripts.common import wall_operators
    from shear_flow_scattering.scripts.trailing_edge import vortex_wall_pressure
    U = u_parabolic(0.3, 1, 1)
    N = int(max(1500, 8 * w / 0.3))
    t0 = time.time()
    p = vortex_wall_pressure(np.asarray(ys), w, 0.0, U, c_const(5.0), 1.0, N, wall_operators(w, U).hard)
    return p, time.time() - t0


def fig07():
    def compute():
        from concurrent.futures import ProcessPoolExecutor
        pw = np.full((OMEGA_7.size, Y2_7.size), np.nan + 0j)
        live = np.ones(Y2_7.size, bool)
        with ProcessPoolExecutor(WORKERS) as pool:
            for i, w in enumerate(OMEGA_7):
                idx = np.flatnonzero(live)
                if idx.size == 0:
                    break
                chunks = [c for c in np.array_split(idx, WORKERS) if c.size]
                for c, (p, dt) in zip(chunks, pool.map(_fig07_row, [(w, Y2_7[c]) for c in chunks])):
                    pw[i, c] = p
                live &= ~(db(pw[i]) < SKIP_DB)
                logging.info("fig07 omega=%.3g: %d heights, %.0f s", w, idx.size, dt)
        return dict(omega=OMEGA_7, y2=Y2_7, pw=pw)

    d = cached("fig07", compute)
    import matplotlib.pyplot as plt
    L = db(d["pw"]).T
    L[~np.isfinite(L)] = -1e3              # not computed: far below the colour scale
    L = np.maximum(L, -120)                # clamp before interpolating
    # display on a fine grid: linear interpolation in (log10 omega, y2)
    lw, y = np.log10(d["omega"]), d["y2"]
    lwf, yf = np.linspace(lw[0], lw[-1], 400), np.linspace(y[0], y[-1], 200)
    Lf = np.array([np.interp(lwf, lw, row) for row in L])
    Lf = np.array([np.interp(yf, y, col) for col in Lf.T]).T
    dl, dy = lwf[1] - lwf[0], yf[1] - yf[0]
    we = 10 ** np.append(lwf - dl / 2, lwf[-1] + dl / 2)    # cell edges, uniform in log omega
    ye = np.append(yf - dy / 2, yf[-1] + dy / 2)
    fig, ax = plt.subplots(figsize=(6, 3.6), constrained_layout=True)
    m = ax.pcolormesh(we, ye, Lf, cmap="viridis", vmin=-100, vmax=0)
    fig.colorbar(m, ax=ax, pad=0.02, label=r"$20\log_{10}|p_w|$")
    ax.set(xscale="log", xlabel=r"$\omega\delta/U_\infty$", ylabel=r"$y_2/\delta$")
    save(fig, "fig07_wall_pressure_vortex_sheet", ["wall_pressure_vortex_sheet"])


# --------------------------------------------------------------------------
# Figs 8 and 9 -- base vs sheared profiles
# --------------------------------------------------------------------------
def fig08_09():
    import fig_paper_8_9 as f89

    def compute():
        res = f89.compute()
        out = dict(f=f89.F_HZ, y2=f89.Y2)
        for (route, name), (pw, T) in res.items():
            out[f"pw_{route}_{name}"], out[f"T_{route}_{name}"] = pw, T
        return out

    d = cached("fig08_09", compute)
    import matplotlib.pyplot as plt
    F, Y2 = d["f"], d["y2"]
    get = lambda route, name: (d[f"pw_{route}_{name}"], d[f"T_{route}_{name}"])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True, sharex=True)
    for col, name in enumerate(f89.CASES):
        for row, what in enumerate(("pw", "pff")):
            ax = axes[row, col]
            for j, y in enumerate(Y2):
                for route in f89.ROUTES:
                    if route == "two-piece" and np.isclose(y, f89.JOIN[name]):
                        continue
                    pw, T = get(route, name)
                    val = f89.masked(route, name, what, pw[:, j] if what == "pw" else pw[:, j] * T[:, j])
                    ax.plot(F, db(val), color=f"C{j}", label=(f"$y_2 = {y:g}\\delta$" if route == "fit" else None),
                            **f89.STYLE[route])
            ax.set(xscale="log", xlim=(F[0], f89.F_MAX), ylim=(-250, 10),
                   title=f"({'abcd'[2 * row + col]}) {'Wall' if what == 'pw' else 'Far-field'} pressure, {name} case",
                   ylabel=r"$20\log_{10}|p_w|$" if what == "pw" else r"$20\log_{10}|p_{ff}|$")
            if row == 1:
                ax.set_xlabel("Frequency, f (Hz)")
    axes[0, 0].legend(fontsize=8)
    save(fig, "fig08_wall_and_far_field_pressure",
         ["wall_pressure_base", "wall_pressure_sheared", "far_field_pressure_base", "far_field_pressure_sheared"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True, sharey=True)
    for ax, what in zip(axes, ("pw", "pff")):
        for j, y in enumerate(Y2):
            for route in f89.ROUTES:
                if route == "two-piece" and any(np.isclose(y, v) for v in f89.JOIN.values()):
                    continue
                (pb, Tb), (ps, Ts) = get(route, "base"), get(route, "sheared")
                r = ps[:, j] / pb[:, j] if what == "pw" else ps[:, j] * Ts[:, j] / (pb[:, j] * Tb[:, j])
                r = f89.masked(route, "sheared", what, r)
                ax.plot(F, db(r), color=f"C{j}", label=(f"$y_2 = {y:g}\\delta$" if route == "fit" else None),
                        **f89.STYLE[route])
        ax.axhline(0, color="k", lw=0.6)
        ax.set(xscale="log", xlim=(F[0], f89.F_MAX), ylim=(-40, 20), xlabel="Frequency, f (Hz)",
               title="(a) Wall pressure, sheared / base" if what == "pw" else "(b) Far field, sheared / base")
    axes[0].set_ylabel("dB")
    axes[0].legend(fontsize=8)
    save(fig, "fig09_change_sheared_vs_base", ["wall_pressure_change", "far_field_change"])


# --------------------------------------------------------------------------
# Fig. 10 -- critical-layer deformation by the semicircular paths
# --------------------------------------------------------------------------
def fig10():
    import matplotlib.pyplot as plt
    from shear_flow_scattering.profiles import u_parabolic

    U, w = u_parabolic(0.2, 1, 1), 1.0
    s = np.linspace(0, 1, 401)
    real = s
    Cp = 0.5 + 0.5 * np.exp(1j * np.pi * (1 - s))    # 0 -> 1 through the upper half-plane
    Cm = np.conj(Cp)
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    for z, lab, st in ((real + 0j, "$z = x_2$", "-"), (Cp, "$z = C^+(s)$", "--"), (Cm, "$z = C^-(s)$", ":")):
        a.plot(z.real, z.imag, "k", ls=st, lw=1.6, label=lab)
        k = w / U.f(z)
        b.plot(k.real, k.imag, "k", ls=st, lw=1.6, label=r"$k_1 = \omega/U(" + lab[5:-1] + ")$")
    a.set(xlabel=r"Re $z$", ylabel=r"Im $z$", title=r"(a) Integration paths in the complex $z$-plane")
    a.set_aspect("equal")
    b.set(xlabel=r"Re $k_1$", ylabel=r"Im $k_1$",
          title=r"(b) Critical layer $k_1 = \omega/U(z)$, $U_0 = 0.2U_\infty$, $\omega = 1$")
    for ax in (a, b):
        ax.legend(fontsize=8, loc="upper right")
    save(fig, "fig10_critical_layer_deformation", ["integration_paths", "critical_layer_images"])


# --------------------------------------------------------------------------
# Fig. 11 -- phase of the hard-wall dispersion function
# --------------------------------------------------------------------------
def fig11():
    def compute():
        # D_h via C+ is analytic except on its cut k1 = omega/U(C+(s)).  Evaluate it
        # with whichever of two other paths keeps clear of the critical point
        # z_c (U(z_c) = omega/k1): both agree with C+ wherever z_c is not between
        # them and C+, so the result is D_h(C+) right up to the cut, and the cut
        # shows as a sharp jump rather than a band of numerical noise.
        #   z_c inside the half-disc bounded by C+  -> tall semi-ellipse (height 0.85)
        #   otherwise                               -> flat semi-ellipse (height 0.2)
        from shear_flow_scattering.contours import semiellipse_path
        from shear_flow_scattering.pridmore_brown import comp_dec_aug18a_skip
        from shear_flow_scattering.profiles import c_const, u_parabolic
        U, c = u_parabolic(0.2, 1, 1), c_const(5.0)
        kr, ki = np.linspace(-3, 6, 451), np.linspace(-1.5, 1.5, 151)
        K = (kr[None, :] + 1j * ki[:, None]).ravel()
        with np.errstate(divide="ignore", invalid="ignore"):
            zc = 1 - np.sqrt((1 - 1.0 / K) / 0.8 + 0j)      # root with Re z_c <= 1 (U = 1 - 0.8(1-z)^2)
        inside = (zc.imag > 0) & (np.abs(zc - 0.5) < 0.5)
        # near the wall (|z_c| small) every path passes close to z_c: refine the RK steps
        wall = np.abs(zc) < 0.15
        D = np.empty(K.size, complex)
        for sel, h in ((inside, 0.85), (~inside, 0.2)):
            for sub, N in ((sel & ~wall, 1000), (sel & wall, 12000)):
                if sub.any():
                    D[sub] = comp_dec_aug18a_skip(K[sub], 0.0, 1.0, U, c, 1.0, N,
                                                  semiellipse_path(1.0, 0.0, +1, h), prepare=False).v
        return dict(kr=kr, ki=ki, D=D.reshape(ki.size, kr.size))

    d = cached("fig11", compute)
    import colorcet as cc
    import matplotlib.pyplot as plt
    # only the cut's end points (k1 = omega/U_inf, omega/U0), where every path
    # meets the critical point, remain unreliable: blank a small disc round each
    K = d["kr"][None, :] + 1j * d["ki"][:, None]
    near = (np.abs(K - 1.0) < 0.04) | (np.abs(K - 5.0) < 0.04)
    ph = np.ma.masked_where(near | ~np.isfinite(d["D"]), np.angle(d["D"]))
    fig, ax = plt.subplots(figsize=(8, 3.4), constrained_layout=True)
    m = ax.pcolormesh(d["kr"], d["ki"], ph, shading="auto", cmap=cc.cm["CET_C6"],
                      vmin=-np.pi, vmax=np.pi)
    cb = fig.colorbar(m, ax=ax, pad=0.02, ticks=[-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    cb.ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    cb.set_label(r"arg $D_h$")
    with np.errstate(divide="ignore"):
        logD = np.log10(np.abs(d["D"]))
    lv = np.arange(np.floor(np.nanpercentile(logD, 2)), np.ceil(np.nanpercentile(logD, 98)) + 1)
    ax.contour(d["kr"], d["ki"], np.ma.masked_where(np.ma.getmaskarray(ph), logD), levels=lv, colors="0.25", linewidths=0.5)
    ax.set(xlabel=r"Re $k_1$", ylabel=r"Im $k_1$", aspect="equal",
           title=r"$D_h = v_d(0)$ via $C^+$, $U_0 = 0.2U_\infty$, $\omega = 1$ (contours: $\log_{10}|D_h|$)")
    save(fig, "fig11_dispersion_function_phase", ["dispersion_function_phase"])


# --------------------------------------------------------------------------
# Proposed section 4.2 figure
# --------------------------------------------------------------------------
def figNew42():
    import sec42_studies as s42

    def compute():
        data = s42.compute()
        out, labels = {}, []
        for i, ((study, label), (th, A)) in enumerate(data.items()):
            out[f"theta_{i}"], out[f"A_{i}"] = th, A
            labels.append(f"{study}|{label}")
        out["labels"] = np.array(labels)
        return out

    d = cached("figNew42", compute)
    import matplotlib.pyplot as plt
    labels = [str(x).split("|") for x in d["labels"]]
    studies = list(dict.fromkeys(s for s, _ in labels))
    titles = {"A: asymmetry": "(a) Asymmetry, $c_0 = 5$", "B: Mach number": "(b) Mach number",
              "C: profile shape": "(c) Profile shape, $c_0 = 5$"}
    fig, axes = plt.subplots(1, 3, subplot_kw={"projection": "polar"}, figsize=(15, 5.2), constrained_layout=True)
    for ax, study in zip(axes, studies):
        for i, (s, lab) in enumerate(labels):
            if s == study:
                ax.plot(d[f"theta_{i}"], d[f"A_{i}"], lw=1.8, label=lab)
        ax.set_title(titles.get(study, study), pad=14)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), fontsize=8, ncol=2, frameon=False)
    save(fig, "figNew42_parameter_studies", ["asymmetry", "mach_number", "profile_shape"])


FIGURES = dict(fig03=fig03, fig04=fig04, fig05=fig05, fig06=fig06, fig07=fig07, fig08_09=fig08_09,
               fig10=fig10, fig11=fig11, figNew42=figNew42)

if __name__ == "__main__":
    sys.path.insert(0, HERE)
    os.chdir(HERE)
    for p in (OUT, TIKZ, CACHE):
        os.makedirs(p, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    import matplotlib
    matplotlib.use("Agg")
    names = [a for a in sys.argv[1:] if a in FIGURES] or list(FIGURES)
    for n in names:
        try:
            FIGURES[n]()
        except Exception:
            logging.exception("%s failed", n)
