"""Corrected paper Fig. 6 (thesis Fig. 6.4.1): transfer of wall pressure to the far field.

Parabolic slipping profile U0 in {0.2, 0.6, 1.0} (1.0 = uniform flow), omega = 5,
unit wall-pressure gust with kappa1 = omega / (0.8 U_inf), hard trailing edge,
(a) Ma = 0.02 (c0 = 50), (b) Ma = 0.2 (c0 = 5).

Solid: corrected (winding cancelled, tanh contours above all wake modes).
Dashed: the original computation (Veitch contour, W = 0) as used for the paper.
Both use the physical far-field amplitude lim r^(1/2)|p'| (see sec42_studies.py).

Usage::

    uv run --extra plot python examples/fig_paper_6.py [--fast] [--replot]
"""
import logging
import sys

import numpy as np

fast = "--fast" in sys.argv
stem = "fig_paper_6" + ("_fast" if fast else "")
OMEGA = 5.0
U0S, C0S = (0.2, 0.6, 1.0), (50.0, 5.0)
N_TH = 61 if fast else 181
KW = dict(N_bl=600 if fast else 1500, prec_N=2000 if fast else 4000)


def compute():
    from sec42_studies import amplitude
    from shear_flow_scattering.profiles import c_const, u_parabolic
    from shear_flow_scattering.scripts.trailing_edge import far_field_trailing_edge_lined

    th1 = np.linspace(0, np.pi, N_TH + 2)[1:-1]
    th2 = np.linspace(np.pi, 2 * np.pi, N_TH + 2)[1:-1]
    out = {}
    for c0 in C0S:
        for U0 in U0S:
            U, c = u_parabolic(U0, 1, 1), c_const(c0)
            for mode, kw in (("corrected", dict(factor_contour="tanh", winding="auto")),
                             ("original", dict(factor_contour="veitch", winding=0))):
                logging.info("c0=%g U0=%g %s", c0, U0, mode)
                ff = far_field_trailing_edge_lined(Up=U, Um=U, cp=c, cm=c, omega=OMEGA, Z=np.inf, theta1=th1,
                                                   theta2=th2, pressure_upstream=1.0, kappa1=OMEGA / 0.8, **kw, **KW)
                th, A = amplitude(th1, th2, ff.FF1p.p, ff.FF2p.p, U, U, c0)
                _, Am = amplitude(th1, th2, ff.FF1m.p, ff.FF2m.p, U, U, c0)
                use_m = (th < np.pi / 2) | (th > 3 * np.pi / 2)
                A[use_m] = Am[use_m]
                out[(c0, U0, mode)] = A
    np.savez(stem + ".npz", theta=th, keys=np.array([f"{c}|{u}|{m}" for (c, u, m) in out]),
             **{f"A_{i}": v for i, v in enumerate(out.values())})
    return th, out


def load():
    d = np.load(stem + ".npz")
    out = {}
    for i, k in enumerate(d["keys"]):
        c, u, m = str(k).split("|")
        out[(float(c), float(u), m)] = d[f"A_{i}"]
    return d["theta"], out


def plot(th, out):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, subplot_kw={"projection": "polar"}, figsize=(11, 5.2), constrained_layout=True)
    for ax, c0, title in zip(axes, C0S, ("(a) Ma = 0.02 (acoustically compact)", "(b) Ma = 0.2")):
        for j, U0 in enumerate(U0S):
            ax.plot(th, out[(c0, U0, "corrected")], color=f"C{j}", lw=2, label=f"$U_0 = {U0}U_\infty$")
            orig = out[(c0, U0, "original")]
            ax.plot(th, orig * np.max(out[(c0, U0, "corrected")]) / np.max(orig), color=f"C{j}", lw=1, ls="--")
        ax.set_title(title, pad=12)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3, fontsize=8, frameon=False)
    fig.suptitle("Wall-pressure-to-far-field transfer, $\omega = 5$, $\kappa_1 = \omega/(0.8U_\infty)$: "
                 "corrected (solid); original, rescaled to the same peak (dashed)", fontsize=10)
    fig.savefig(stem + ".png", dpi=130)
    print("saved", stem + ".png")
    for c0 in C0S:
        for U0 in U0S:
            a, o = out[(c0, U0, "corrected")], out[(c0, U0, "original")]
            i90 = np.argmin(np.abs(th - np.pi / 2))
            print(f"c0={c0:4.0f} U0={U0}: peak corrected {a.max():.4f}, original {o.max():.4f} "
                  f"(x{a.max()/o.max():.1f}); at 90deg corrected {a[i90]:.4f}, original {o[i90]:.4f}")


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    plot(*(load() if "--replot" in sys.argv else compute()))
