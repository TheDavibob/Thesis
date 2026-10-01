"""Compare the thesis (6.3.10) source weighting with the literal Jacobian.

Same setup as ``fig_distributed.py`` (``Distributed_TrailingEdge.m``
defaults).  The expensive, weighting-independent setup runs once; each panel
is then a cheap ``rerun_with_vorticity``.  Shows Re p' (total) for

* rows: ``omega_c = 1`` (the thesis figure) and ``omega_c = y2 (delta - y2)``
  (the example in ``Distributed_TrailingEdge_VarOmegaC.m``);
* columns: ``weighting="thesis"`` (as the MATLAB) and ``"jacobian"``.

``omega_c = 1`` with ``"jacobian"`` is log-divergent at the boundary-layer
edge (``U'(delta) = 0``); it is plotted to show that, and will change with
``N_distribution``.  The script also prints the change in that panel between
N_distribution = 40 and 80 as a check.

Usage::

    uv run --extra plot python examples/fig_weighting_comparison.py [--fast]
"""
import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

from shear_flow_scattering.scripts.trailing_edge import (distributed_trailing_edge,
                                                         rerun_with_vorticity)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
fast = "--fast" in sys.argv
grid = dict(x1=np.arange(-5, 5 + 1e-9, 0.1 if fast else 0.05),
            x2=np.arange(-2, 2 + 1e-9, 0.04 if fast else 0.02),
            N_bl=500 if fast else 2000, N_C=200 if fast else 500, N_CL=400 if fast else 1000,
            prec_N=1000 if fast else 2000)
N_dist = 30 if fast else 80

base = distributed_trailing_edge(N_distribution=N_dist, **grid)
st, y2 = base.state, base.y2
dists = {r"$\omega_c = 1$": np.ones_like(y2), r"$\omega_c = y_2(\delta - y_2)$": y2 * (1 - y2)}

results = {}
for dname, oc in dists.items():
    for w in ("thesis", "jacobian"):
        results[dname, w] = rerun_with_vorticity(st, oc, weighting=w).I_total.p.real.T

stem = "fig_weighting_comparison" + ("_fast" if fast else "")
np.savez(stem + ".npz", x1=base.x1, x2=base.x2, y2=y2,
         **{f"p_total_{i}_{w}": results[d, w] for i, d in enumerate(dists) for w in ("thesis", "jacobian")})

fig, axes = plt.subplots(2, 2, figsize=(12, 6), constrained_layout=True)
for i, dname in enumerate(dists):
    for j, w in enumerate(("thesis", "jacobian")):
        ax, data = axes[i, j], results[dname, w]
        lo, hi = np.percentile(data, [1, 99])
        mesh = ax.pcolormesh(base.x1, base.x2, data, shading="auto", cmap="viridis", vmin=lo, vmax=hi)
        fig.colorbar(mesh, ax=ax, pad=0.02)
        ax.plot([-5, 0], [0, 0], "k", lw=4)
        note = "  (log-divergent!)" if (i == 0 and w == "jacobian") else ""
        ax.set(xlabel="$x_1$", ylabel="$x_2$", title=f"Re $p'$ total, {dname}, weighting='{w}'{note}")
fig.savefig(stem + ".png", dpi=120)
print("saved", stem + ".png")

# Convergence check in N_distribution for the y2(delta - y2) case, both weightings
if not fast:
    coarse = distributed_trailing_edge(N_distribution=40, **grid)
    for w in ("thesis", "jacobian"):
        for dname, oc_fn in ((r"omega_c=1", np.ones_like), ("omega_c=y2(1-y2)", lambda y: y * (1 - y))):
            a = rerun_with_vorticity(coarse.state, oc_fn(coarse.y2), weighting=w).I_total.p
            b = rerun_with_vorticity(st, oc_fn(y2), weighting=w).I_total.p
            print(f"{w:9s} {dname:18s} max|p(N=40) - p(N=80)| / max|p(N=80)| = "
                  f"{np.abs(a - b).max() / np.abs(b).max():.3f}")
