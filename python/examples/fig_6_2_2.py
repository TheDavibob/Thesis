"""Regenerate thesis Fig. 6.2.2: vortex sheet scattered by a hard trailing edge.

Source: ``TrailingEdge_Lined.m`` (single vortex sheet at ``y2``), with the flow
set to the one described in thesis section 6.2.5 (p. 226):

* hard plate (``Z = inf``), ``omega = 1``, ``delta = 1`` both sides, ``y2 = 0.5``;
* ``Struct_U_Parabolic(0.3, 1, 1)`` above and ``Struct_U_Parabolic(0.3, 1.5, 1)``
  below -- parabolas with ``U'(delta) = 0`` ("shear smooth at the edge of the
  boundary layer");
* ``c0 = 5`` (not stated for this figure in the thesis).

The profile and ``c0`` are confirmed by the wake mode: this setup gives
``kappa_w = 0.8720 +- 0.2618i``, matching the value quoted in the thesis to
all four digits (linear profiles, or c0 = 2 or 10, do not).

Usage::

    uv run --extra plot python examples/fig_6_2_2.py [--fast] [--tikz]
"""
import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

from shear_flow_scattering.profiles import c_const, u_parabolic
from shear_flow_scattering.scripts.trailing_edge import trailing_edge_lined

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
fast = "--fast" in sys.argv

res = trailing_edge_lined(
    Up=u_parabolic(0.3, 1.0, 1.0), Um=u_parabolic(0.3, 1.5, 1.0),
    cp=c_const(5.0), cm=c_const(5.0), deltap=1.0, deltam=1.0,
    omega=1.0, Z=np.inf, y2=0.5,
    incident="original",  # as the thesis; "regularised" fixes the overall phase/amplitude
    x1=np.arange(-5, 10 + 1e-9, 0.1 if fast else 0.05),
    x2=np.arange(-2, 2 + 1e-9, 0.04 if fast else 0.02),
    N_bl=500 if fast else 2000, N_C=200 if fast else 500, N_CL=400 if fast else 1000,
    prec_N=1000 if fast else 2000,
)
print("wake modes:", res.zeros_wake_m)
np.savez("fig_6_2_2" + ("_fast" if fast else "") + ".npz", x1=res.x1, x2=res.x2,
         **{f"{w}_{f}": res[w][f] for w in ("I_total", "I_scattered") for f in ("p", "v")})

panels = [("I_total", "p", r"(a) Re $p'_0$, total"), ("I_scattered", "p", r"(b) Re $p'_0$, scattered"),
          ("I_total", "v", r"(c) Re $v'_0$, total"), ("I_scattered", "v", r"(d) Re $v'_0$, scattered")]
fig, axes = plt.subplots(4, 1, figsize=(6, 11), constrained_layout=True)
for ax, (which, f, title) in zip(axes, panels):
    data = res[which][f].real.T
    lo, hi = np.percentile(data, [1, 99])  # clip the few extreme pixels near the critical layer
    mesh = ax.pcolormesh(res.x1, res.x2, data, shading="auto", cmap="viridis", vmin=lo, vmax=hi)
    fig.colorbar(mesh, ax=ax, pad=0.02)
    ax.plot([-5, 0], [0, 0], "k", lw=4)
    ax.set(xlabel="$x_1$", ylabel="$x_2$", title=title)
out = "fig_6_2_2_fast.png" if fast else "fig_6_2_2.png"
fig.savefig(out, dpi=130)
if "--tikz" in sys.argv:  # one matlab2tikz-style .tex per panel
    from shear_flow_scattering.tikz import save_figure_tikz
    print("wrote", *save_figure_tikz(fig, out[:-4]))
print("saved", out)
