"""Complex logarithms / powers with movable branch cuts.

Ports of ``CommonTools/LogA.m``, ``SqrtA.m``, ``PowerA.m``, ``PowerB.m``,
``Log_CurvedBranch.m``, ``LogCts.m``, ``LogCtsF.m`` and ``NaN2Zero.m``.

Branch-cut convention: ``cut`` is the angle (radians) of the ray from the
origin along which the cut lies.  ``cut = pi`` recovers the principal branch.
Note that ``-pi/2`` and ``3*pi/2`` describe the same ray but give values that
differ by ``2*pi*i`` in the log (this matters, and is used deliberately in
the Wiener-Hopf factorisations: ``3*pi/2`` for "+" functions analytic in the
upper half-plane, ``pi/2`` for "-" functions).

.. note::
   NumPy honours signed zeros on the principal cut (``sqrt(-4-0j) = -2j``),
   whereas MATLAB returns ``+2i``.  Points lying *exactly* on a cut may
   therefore differ from the MATLAB results; points off the cut agree.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from .quadrature import adaptive_integral


def log_a(x2, cut):
    """``LogA.m`` -- ``log(x2)`` with the branch cut along ``arg z = cut``."""
    alpha = 1j * (cut - np.pi)
    return np.log(np.exp(-alpha) * np.asarray(x2, dtype=complex)) + alpha


def sqrt_a(x2, cut):
    """``SqrtA.m`` -- ``sqrt(x2)`` with the branch cut along ``arg z = cut``."""
    alpha = 1j * (cut - np.pi)
    return np.exp(alpha / 2) * np.sqrt(np.exp(-alpha) * np.asarray(x2, dtype=complex))


def power_a(x2, power, cut):
    """``PowerA.m`` -- ``x2**power`` using :func:`log_a` (constant power)."""
    return np.exp(power * log_a(x2, cut))


def power_b(x2, power: Callable, cut):
    """``PowerB.m`` -- ``x2**power(x2)``: exponent is itself a function."""
    return np.exp(power(x2) * log_a(x2, cut))


def log_curved_branch(z, f: Callable):
    """``Log_CurvedBranch.m`` -- log with a *curved* cut ``theta = f(|z|)``.

    The returned argument lies in ``(f(r), f(r) + 2*pi]``.
    """
    z = np.asarray(z, dtype=complex)
    r = np.abs(z)
    theta = np.pi + f(r) + np.angle(z * np.exp(-1j * (f(r) + np.pi)))
    return np.log(r) + 1j * theta


def nan_to_zero(f: Callable, x):
    """``NaN2Zero.m`` -- evaluate ``f(x)`` and replace NaNs by 0."""
    g = np.array(f(x), dtype=complex)
    g[np.isnan(g)] = 0
    return g


def log_cts(t, C: Callable, dC: Callable):
    """``LogCts.m`` -- ``log(C(t))`` continuous along the contour ``z = C(t)``.

    Requires ``C(0) = 1``; computed as ``int_0^t C'(eta)/C(eta) d eta`` so it
    is independent of any branch choice (fails if the contour hits 0).
    """
    return adaptive_integral(lambda eta: dC(eta) / C(eta), 0.0, t)


def log_cts_f(t, F: Callable, dF: Callable, Cinv: Callable, dCinv: Callable):
    """``LogCtsF.m`` -- ``log(F)`` continuous along a contour, ``F -> 1`` at ``-inf``.

    Integrates ``F'(Cinv(eta)) Cinv'(eta) / F(Cinv(eta))`` from ``-inf`` to
    ``t``.  (The MATLAB comment and variable names are ambiguous about whether
    ``Cinv`` is the contour or its inverse; this is a literal port.)
    """
    g = lambda eta: dF(Cinv(eta)) * dCinv(eta) / F(Cinv(eta))
    return adaptive_integral(g, -np.inf, t)
