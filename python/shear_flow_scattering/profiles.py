"""Mean-flow profiles ``U(x2)`` and sound-speed profiles ``c(x2)``.

Ports of ``FixedStep/Flow_Structures/Struct_*.m``.  Each returns a
:class:`Profile` with callables ``f`` (value), ``df``, ``d2f`` (and ``d3f``
for the cubic) -- the same fields the MATLAB structs had, so ``U.f(x2)``,
``U.df(x2)`` read identically.

All profiles accept complex ``x2``: the Pridmore-Brown solver integrates along
complex paths to avoid the critical layer.  Piecewise definitions compare the
*real part* of ``x2``, which is what MATLAB's ``>`` does on complex input.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from ._utils import Struct, log


class Profile(Struct):
    """Struct with ``f``, ``df``, ``d2f`` (optionally ``d3f``) callables."""

    def __init__(self, f: Callable, df: Callable, d2f: Optional[Callable] = None,
                 d3f: Optional[Callable] = None):
        super().__init__(f=f, df=df)
        if d2f is not None:
            self.d2f = d2f
        if d3f is not None:
            self.d3f = d3f


def _ones(x):
    return np.ones(np.shape(x))


def u_const(U0) -> Profile:
    """``Struct_U_Const.m``: uniform flow ``U = U0``."""
    return Profile(lambda x2: U0 * _ones(x2), lambda x2: 0 * _ones(x2),
                   lambda x2: 0 * _ones(x2))


def u_linear(U0, Uinf, deltainf) -> Profile:
    """``Struct_U_Linear.m``: ``U(0) = U0`` rising linearly to ``U(deltainf) = Uinf``.

    Note the shear does *not* vanish at ``deltainf``; the outer-region
    formulae elsewhere assume ``U'(delta) = 0`` (see README, "Assumptions").
    """
    s = (Uinf - U0) / deltainf
    return Profile(lambda x2: U0 + s * np.asarray(x2), lambda x2: s * _ones(x2),
                   lambda x2: 0 * _ones(x2))


def u_parabolic(U0, Uinf, deltainf) -> Profile:
    """``Struct_U_Parabolic.m``: parabola with ``U(0)=U0``, ``U(delta)=Uinf``,
    ``U'(delta) = 0`` (continuous shear at the boundary-layer edge)."""
    A = (Uinf - U0) / deltainf ** 2
    return Profile(lambda x2: Uinf - A * (np.asarray(x2) - deltainf) ** 2,
                   lambda x2: -2 * A * (np.asarray(x2) - deltainf),
                   lambda x2: -2 * A * _ones(x2))


def u_parabolic_b(u, delta) -> Profile:
    """``Struct_U_Parabolic_B.m``: parabola through three points ``U(delta_j) = u_j``."""
    C = ((u[2] - u[0]) / (delta[2] - delta[0])
         - (u[1] - u[0]) / (delta[1] - delta[0])) / (delta[2] - delta[1])
    B = (u[1] - u[0]) / (delta[1] - delta[0]) - C * (delta[1] - delta[0])
    if np.sign(B) != np.sign(B + 2 * C * (delta[2] - delta[0])):
        log.warning("Parabola is not monotonic across the boundary-layer")
    return Profile(lambda x: u[0] + B * (np.asarray(x) - delta[0]) + C * (np.asarray(x) - delta[0]) ** 2,
                   lambda x: B + 2 * C * (np.asarray(x) - delta[0]),
                   lambda x: 2 * C * _ones(x))


def u_pw_parabolic(U_points, delta_points) -> Profile:
    """``Struct_U_PWParabolic.m``: two parabolas joined at ``delta_points[1]``
    with continuous value and shear; ``U'(deltainf) = 0``."""
    U0, U1, Uinf = U_points
    d0, d1, dinf = delta_points
    A = (U1 - Uinf) / (d1 - dinf) ** 2
    sigma = -2 * A * (dinf - d1)
    B = ((U0 - U1) - sigma * (d0 - d1)) / (d0 - d1) ** 2

    def pick(x2, above, below):
        x2 = np.asarray(x2)
        up = np.real(x2) > d1
        return np.where(up, above(x2), below(x2))

    return Profile(
        lambda x2: pick(x2, lambda x: U1 + sigma * (x - d1) + A * (x - d1) ** 2,
                        lambda x: U1 + sigma * (x - d1) + B * (x - d1) ** 2),
        lambda x2: pick(x2, lambda x: sigma + 2 * A * (x - d1),
                        lambda x: sigma + 2 * B * (x - d1)),
        lambda x2: pick(x2, lambda x: 2 * A * _ones(x), lambda x: 2 * B * _ones(x)))


def u_cubic(U0, sigma0) -> Profile:
    """``Struct_U_Cubic.m``: cubic with ``U(0)=U0``, ``U'(0)=sigma0`` and
    ``U(1) = 1``, ``U'(1) = 0`` (``delta = Uinf = 1``).

    Monotonic if ``3 U0 + sigma0 < 3`` (and ``sigma0 > 0``); inflection-free if
    additionally ``3 U0 + 2 sigma0 > 3``.
    """
    A, B = U0, sigma0
    C = 3 - 3 * U0 - 2 * sigma0
    D = -2 + 2 * U0 + sigma0
    return Profile(lambda x: A + B * np.asarray(x) + C * np.asarray(x) ** 2 + D * np.asarray(x) ** 3,
                   lambda x: B + 2 * C * np.asarray(x) + 3 * D * np.asarray(x) ** 2,
                   lambda x: 2 * C + 6 * D * np.asarray(x),
                   lambda x: 6 * D * _ones(x))


def u_cubic_b(U_points, delta_points) -> Profile:
    """``Struct_U_Cubic_B.m``: cubic through three points, ``U'(deltainf) = 0``."""
    U0, U1, Uinf = U_points
    d0, d1, dinf = delta_points
    B = ((U1 - Uinf) / ((d1 - dinf) ** 2 * (d1 - d0))
         - (U0 - Uinf) / ((d0 - dinf) ** 2 * (d1 - d0)))
    A = (U1 - Uinf) / (d1 - dinf) ** 2 - B * (d1 - dinf)
    return Profile(lambda x: Uinf + A * (np.asarray(x) - dinf) ** 2 + B * (np.asarray(x) - dinf) ** 3,
                   lambda x: 2 * A * (np.asarray(x) - dinf) + 3 * B * (np.asarray(x) - dinf) ** 2,
                   lambda x: 2 * A + 6 * B * (np.asarray(x) - dinf))


def u_quartic(U_points, delta_points) -> Profile:
    """``Struct_U_Quartic.m``: even quartic about ``deltainf`` through three points."""
    U0, U1, Uinf = U_points
    d0, d1, dinf = delta_points
    den = (d1 - dinf) ** 2 - (d0 - dinf) ** 2
    B = ((U1 - Uinf) / ((d1 - dinf) ** 2 * den)
         - (U0 - Uinf) / ((d0 - dinf) ** 2 * den))
    A = (U1 - Uinf) / (d1 - dinf) ** 2 - B * (d1 - dinf) ** 2
    return Profile(lambda x: Uinf + A * (np.asarray(x) - dinf) ** 2 + B * (np.asarray(x) - dinf) ** 4,
                   lambda x: 2 * A * (np.asarray(x) - dinf) + 4 * B * (np.asarray(x) - dinf) ** 3,
                   lambda x: 2 * A + 12 * B * (np.asarray(x) - dinf) ** 2)


def c_const(c0) -> Profile:
    """``Struct_c_Const.m``: uniform sound speed (``f`` and ``df`` only)."""
    return Profile(lambda x2: c0 * _ones(x2), lambda x2: 0 * _ones(x2))


def smooth_parabola_plot(ax=None):
    """``SmoothParabolaPlot.m``: the family of parabolic profiles, U0 = 0.5 ... 0.

    Requires matplotlib (``uv sync --extra plot``).
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    x2a = np.arange(0, 1.0001, 0.01)
    x2b = np.arange(1.01, 2.0001, 0.01)
    for j, U0 in enumerate(np.arange(0.5, -0.05, -0.1)):
        U = u_parabolic(U0, 1, 1)
        prof = np.concatenate([U.f(x2a), U.f(1) * np.ones_like(x2b)])
        ax.plot(prof, np.concatenate([x2a, x2b]), lw=2,
                ls="-" if j == 0 else "--", color="C0")
    ax.axis([0, 1.2, 0, 2])
    ax.set_xlabel("$U(x_2)$")
    ax.set_ylabel("$x_2$")
    return ax
