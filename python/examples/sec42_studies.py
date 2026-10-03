"""Parameter studies for paper section 4.2: asymmetry, Mach number, profile shape.

Setup as paper section 4.3.1 (Fig. 6): hard trailing edge, omega = 5, forced by a
unit wall-pressure gust (P = 1) convecting at Uc = 0.8 U_inf^+, so
kappa1 = omega / (0.8 U_inf^+).  This isolates the transfer from wall
pressure to far-field sound.

Plotted quantity: the physical far-field amplitude

    A(theta) = lim r^(1/2) |p'|
             = |FF| sqrt(K0 / 2 pi) (cos^2 theta + beta^2 sin^2 theta)^(-1/4),

i.e. the saddle-point result *including* the factors the original far-field
scripts dropped: 1/(2 pi) from the Fourier inversion, sqrt(pi/(R K0)), and
the Mach-stretched distance R.  Without them, magnitudes cannot be compared
across Mach numbers.  (Checked against direct near-field inversion to ~1%.)

Studies (all U0 = 0.3, delta = 1 both sides unless stated):

A. asymmetry -- U_inf^- in {1, 1.25, 1.5, 2} below the wake, c0 = 5;
B. Mach      -- c0 in {5, 10, 20, 50} (Ma = 0.2 ... 0.02), symmetric;
C. profiles  -- same slip U0 = 0.3, increasing wall shear: a cubic with
                U'(0) = 1.1, then power laws U = 1 - 0.7 (1 - x2)^n with
                n = 2 (parabolic), 3, 5, 8, i.e. U'(0) = 1.4 ... 5.6; all
                monotonic, inflection-free, symmetric, c0 = 5.  (The paper's
                two-piece profiles are not analytic, so cannot be used with
                the complex-path method.)

Usage::

    uv run --extra plot python examples/sec42_studies.py [--fast] [--replot] [--tikz]
"""
import logging
import sys

import numpy as np

fast = "--fast" in sys.argv
stem = "sec42_studies" + ("_fast" if fast else "")
OMEGA, UC_RATIO = 5.0, 0.8
N_TH = 61 if fast else 181
KW = dict(N_bl=600 if fast else 1500, prec_N=2000 if fast else 4000)


def cases():
    from shear_flow_scattering.profiles import u_cubic, u_parabolic, u_power
    par = u_parabolic(0.3, 1, 1)
    A = {f"$U_\\infty^- = {u}$": dict(Up=par, Um=u_parabolic(0.3, u, 1), c0=5.0) for u in (1.0, 1.25, 1.5, 2.0)}
    B = {f"$c_0 = {c:g}$ (Ma = {1/c:.2g})": dict(Up=par, Um=par, c0=c) for c in (5.0, 10.0, 20.0, 50.0)}
    C = {"cubic, $U'(0)=1.1$": dict(Up=u_cubic(0.3, 1.1), Um=u_cubic(0.3, 1.1), c0=5.0)}
    for n in (2, 3, 5, 8):
        Un = u_power(0.3, 1, 1, n)
        C[f"$n={n}$, $U'(0)={0.7*n:.2g}$"] = dict(Up=Un, Um=Un, c0=5.0)
    # (the paper's two-piece U22 profiles are not analytic, so the complex-path
    # method cannot be used with them -- see README)
    return {"A: asymmetry": A, "B: Mach number": B, "C: profile shape": C}


def amplitude(theta1, theta2, FF1, FF2, Up, Um, c0):
    from shear_flow_scattering.steepest_descent import far_field_saddle
    out = []
    for th, FF, U, lower in ((theta1, FF1, Up, False), (theta2, FF2, Um, True)):
        sd = far_field_saddle(th, OMEGA, U.f(1.0), c0, lower=lower)
        out.append(np.abs(np.ravel(FF)) * np.sqrt(sd.K0 / (2 * np.pi))
                   * (np.cos(th) ** 2 + sd.beta ** 2 * np.sin(th) ** 2) ** -0.25)
    return np.concatenate([theta1, theta2]), np.concatenate(out)


def theta_grids(n=N_TH):
    """Observer angles above (theta1) and below (theta2) the plate: uniform,
    plus log-clustered points next to theta = pi, where the directivity has a
    very narrow grazing null (half-width ~0.1 deg at Ma = 0.02).  The end
    points 0, pi, 2 pi are avoided (saddle on the acoustic branch point)."""
    d = np.logspace(-3.7, -0.6, 40)
    th1 = np.union1d(np.linspace(0, np.pi, n + 2)[1:-1], np.pi - d)
    th2 = np.union1d(np.linspace(np.pi, 2 * np.pi, n + 2)[1:-1], np.pi + d)
    return th1, th2


def compute():
    from shear_flow_scattering.profiles import c_const
    from shear_flow_scattering.scripts.trailing_edge import far_field_trailing_edge_lined

    th1, th2 = theta_grids()
    data = {}
    for study, cs in cases().items():
        for label, p in cs.items():
            logging.info("%s -- %s", study, label)
            c = c_const(p["c0"])
            ff = far_field_trailing_edge_lined(Up=p["Up"], Um=p["Um"], cp=c, cm=c, omega=OMEGA, Z=np.inf,
                                               theta1=th1, theta2=th2, pressure_upstream=1.0,
                                               kappa1=OMEGA / (UC_RATIO * p["Up"].f(1.0)),
                                               factor_contour="tanh", **KW)
            th, A = amplitude(th1, th2, ff.FF1p.p, ff.FF2p.p, p["Up"], p["Um"], p["c0"])
            # use the downstream (-) evaluation where appropriate, as the far-field script does
            th_m, A_m = amplitude(th1, th2, ff.FF1m.p, ff.FF2m.p, p["Up"], p["Um"], p["c0"])
            use_m = (th < np.pi / 2) | (th > 3 * np.pi / 2)
            A[use_m] = A_m[use_m]
            data[study, label] = (th, A, ff.K.W)
    keys = list(data)
    np.savez(stem + ".npz", study=np.array([k[0] for k in keys]), label=np.array([k[1] for k in keys]),
             W=np.array([data[k][2] for k in keys]),
             **{f"theta_{i}": data[k][0] for i, k in enumerate(keys)},
             **{f"A_{i}": data[k][1] for i, k in enumerate(keys)})
    return {k: v[:2] for k, v in data.items()}


def load():
    d = np.load(stem + ".npz")
    return {(str(s), str(l)): (d[f"theta_{i}"], d[f"A_{i}"])
            for i, (s, l) in enumerate(zip(d["study"], d["label"]))}


def plot(data):
    import matplotlib.pyplot as plt

    studies = list(dict.fromkeys(s for s, _ in data))
    fig, axes = plt.subplots(1, len(studies), subplot_kw={"projection": "polar"}, figsize=(15, 5.2),
                             constrained_layout=True)
    for ax, study in zip(axes, studies):
        for (s, label), (th, A) in data.items():
            if s == study:
                ax.plot(th, A, lw=1.8, label=label)
        ax.set_title(study, pad=14)
        ax.set_theta_zero_location("E")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), fontsize=8, ncol=2, frameon=False)
    fig.suptitle(r"Far-field amplitude $\lim r^{1/2}|p'|$ for a unit wall-pressure gust, "
                 r"$\omega = 5$, $\kappa_1 = \omega/(0.8U_\infty^+)$, hard trailing edge", fontsize=11)
    fig.savefig(stem + ".png", dpi=130)
    if "--tikz" in sys.argv:  # one matlab2tikz-style .tex per panel
        from shear_flow_scattering.tikz import save_figure_tikz
        print("wrote", *save_figure_tikz(fig, stem))
    print("saved", stem + ".png")

    print("\nAmplitude at selected angles, and dB relative to the first case of each study:")
    for study in studies:
        rows = [(l, th, A) for (s, l), (th, A) in data.items() if s == study]
        print(f"  {study}")
        ref = None
        for l, th, A in rows:
            vals = [A[np.argmin(np.abs(th - np.radians(d)))] for d in (45, 90, 135, 225, 270, 315)]
            ref = vals if ref is None else ref
            db = 20 * np.log10(np.array(vals) / np.array(ref))
            print(f"    {l:32s} A(45,90,135,225,270,315deg) = " + " ".join(f"{v:.4f}" for v in vals)
                  + "   dB: " + " ".join(f"{x:+5.1f}" for x in db))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    plot(load() if "--replot" in sys.argv else compute())
