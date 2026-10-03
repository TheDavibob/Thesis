"""Compare the thesis (6.3.10) source weighting with the literal Jacobian.

Same setup as ``fig_distributed.py`` (``Distributed_TrailingEdge.m``
defaults).  The expensive, weighting-independent setup runs once; each panel
is then a cheap ``rerun_with_vorticity``.  Shows Re p' (total) for

* rows: ``omega_c = 1`` (the thesis figure) and ``omega_c = y2 (delta - y2)``
  (the example in ``Distributed_TrailingEdge_VarOmegaC.m``);
* columns: ``weighting="thesis"`` (as the MATLAB) and ``"jacobian"``.

All four cases are converged in ``N_distribution``: the script finishes by
comparing N = 40 against N = 80 (relative change < 0.1% at full resolution).
The ``"jacobian"`` weight grows like 1/U'^2 towards the boundary-layer edge,
but the unit-source wall pressure I0 vanishes there like (delta - y2)^4.

Usage::

    uv run --extra plot python examples/fig_weighting_comparison.py [--fast] [--tikz]
    uv run --extra plot python examples/fig_weighting_comparison.py --replot [--fast] [--tikz]

``--replot`` redraws from the saved ``.npz`` without recomputing.
"""
import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

fast = "--fast" in sys.argv
stem = "fig_weighting_comparison" + ("_fast" if fast else "")
DIST_LABELS = (r"$\omega_c = 1$", r"$\omega_c = y_2(\delta - y_2)$")
WEIGHTINGS = ("thesis", "jacobian")


def plot(x1, x2, results):
    fig, axes = plt.subplots(2, 2, figsize=(12, 6), constrained_layout=True)
    for i, dname in enumerate(DIST_LABELS):
        for j, w in enumerate(WEIGHTINGS):
            ax, data = axes[i, j], results[i, w]
            lo, hi = np.percentile(data, [1, 99])
            mesh = ax.pcolormesh(x1, x2, data, shading="auto", cmap="viridis", vmin=lo, vmax=hi)
            fig.colorbar(mesh, ax=ax, pad=0.02)
            ax.plot([-5, 0], [0, 0], "k", lw=4)
            ax.set(xlabel="$x_1$", ylabel="$x_2$", title=f"Re $p'$ total, {dname}, weighting='{w}'")
    fig.savefig(stem + ".png", dpi=120)
    if "--tikz" in sys.argv:  # one matlab2tikz-style .tex per panel
        from shear_flow_scattering.tikz import save_figure_tikz
        print("wrote", *save_figure_tikz(fig, stem))
    print("saved", stem + ".png")


if "--replot" in sys.argv:
    d = np.load(stem + ".npz")
    plot(d["x1"], d["x2"], {(i, w): d[f"p_total_{i}_{w}"] for i in range(2) for w in WEIGHTINGS})
    sys.exit()

from shear_flow_scattering.scripts.trailing_edge import (distributed_trailing_edge,  # noqa: E402
                                                         rerun_with_vorticity)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
grid = dict(x1=np.arange(-5, 5 + 1e-9, 0.1 if fast else 0.05),
            x2=np.arange(-2, 2 + 1e-9, 0.04 if fast else 0.02),
            N_bl=500 if fast else 2000, N_C=200 if fast else 500, N_CL=400 if fast else 1000,
            prec_N=1000 if fast else 2000)
dist_fns = (np.ones_like, lambda y: y * (1 - y))

base = distributed_trailing_edge(N_distribution=30 if fast else 80, **grid)
results = {(i, w): rerun_with_vorticity(base.state, f(base.y2), weighting=w).I_total.p.real.T
           for i, f in enumerate(dist_fns) for w in WEIGHTINGS}
np.savez(stem + ".npz", x1=base.x1, x2=base.x2, y2=base.y2,
         **{f"p_total_{i}_{w}": v for (i, w), v in results.items()})
plot(base.x1, base.x2, results)

if not fast:  # convergence in N_distribution
    coarse = distributed_trailing_edge(N_distribution=40, **grid)
    for w in WEIGHTINGS:
        for i, f in enumerate(dist_fns):
            a = rerun_with_vorticity(coarse.state, f(coarse.y2), weighting=w).I_total.p
            b = rerun_with_vorticity(base.state, f(base.y2), weighting=w).I_total.p
            print(f"{w:9s} case {i}: max|p(N=40) - p(N=80)| / max|p(N=80)| = "
                  f"{np.abs(a - b).max() / np.abs(b).max():.2e}")
