"""Regenerate the distributed-vorticity trailing-edge figure
("a general convected harmonic disturbance is scattered off a hard
trailing-edge, with omega_c = 1 for x2 in (0, 1)").

Source: ``Distributed_TrailingEdge.m``, whose defaults are exactly this case:
the same flow as Fig. 6.2.2 (hard plate, omega = 1, delta = 1,
``Struct_U_Parabolic(0.3, 1, 1)`` above / ``(0.3, 1.5, 1)`` below, c0 = 5),
forced by a uniform vorticity distribution ``omega_c = 1`` across the upper
boundary layer (80 source heights, ``shift_y2 = 0.01``).

Plots Re p and Re v, total and scattered, with colourbars.

Usage::

    uv run --extra plot python examples/fig_distributed.py [--fast]
"""
import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

from shear_flow_scattering.scripts.trailing_edge import distributed_trailing_edge

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
fast = "--fast" in sys.argv

res = distributed_trailing_edge(
    x1=np.arange(-5, 5 + 1e-9, 0.1 if fast else 0.05),
    x2=np.arange(-2, 2 + 1e-9, 0.04 if fast else 0.02),
    N_bl=500 if fast else 2000, N_C=200 if fast else 500, N_CL=400 if fast else 1000,
    N_distribution=30 if fast else 80, prec_N=1000 if fast else 2000,
)
print("wake modes:", res.zeros_wake_m)
np.savez("fig_distributed" + ("_fast" if fast else "") + ".npz", x1=res.x1, x2=res.x2,
         **{f"{w}_{f}": res[w][f] for w in ("I_total", "I_scattered") for f in ("p", "v")})

panels = [("I_total", "p", r"Re $p'$, total"), ("I_scattered", "p", r"Re $p'$, scattered"),
          ("I_total", "v", r"Re $v'$, total"), ("I_scattered", "v", r"Re $v'$, scattered")]
fig, axes = plt.subplots(4, 1, figsize=(6, 11), constrained_layout=True)
for ax, (which, f, title) in zip(axes, panels):
    data = res[which][f].real.T
    lo, hi = np.percentile(data, [1, 99])  # clip the few extreme pixels near the critical layer
    mesh = ax.pcolormesh(res.x1, res.x2, data, shading="auto", cmap="viridis", vmin=lo, vmax=hi)
    fig.colorbar(mesh, ax=ax, pad=0.02)
    ax.plot([-5, 0], [0, 0], "k", lw=4)
    ax.set(xlabel="$x_1$", ylabel="$x_2$", title=title)
out = "fig_distributed_fast.png" if fast else "fig_distributed.png"
fig.savefig(out, dpi=130)
print("saved", out)
