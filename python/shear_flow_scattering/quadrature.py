"""Quadrature rules and (adaptive) complex line integrals.

Ports of ``lgwt/lgwt.m`` and ``CommonTools/GaussQuad/Gauss{Legendre,
Hermite,Laguerre}.m``, plus NumPy replacements for MATLAB's adaptive
``integral`` (used by ``LogCts``, ``BoxIntegral``, ``ZerosPoles`` ...).
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np


# --------------------------------------------------------------------------
# Fixed rules
# --------------------------------------------------------------------------
def lgwt(N: int, a: float, b: float):
    """``lgwt.m`` -- N-point Gauss-Legendre nodes/weights on ``[a, b]``.

    **Ordering matches MATLAB:** nodes run from ``b`` down to ``a``.  Much of
    the original code immediately reverses them (``t(end:-1:1)``); those
    reversals are kept in the port, so this ordering must be preserved.
    If ``a > b`` the weights are negative (as in MATLAB), which
    ``kappa1_quad`` relies upon.
    """
    y, w = np.polynomial.legendre.leggauss(N)  # ascending on [-1, 1]
    y, w = y[::-1], w[::-1]                    # MATLAB lgwt: descending
    x = (a * (1 - y) + b * (1 + y)) / 2
    return x, w * (b - a) / 2


def gauss_legendre(n: int):
    """``GaussLegendre.m`` (Van Damme) -- Golub-Welsch on ``[-1, 1]``, ascending."""
    i = np.arange(1, n)
    a = i / np.sqrt(4 * i ** 2 - 1)
    x, V = np.linalg.eigh(np.diag(a, 1) + np.diag(a, -1))
    return x, 2 * V[0, :] ** 2


def gauss_hermite(n: int):
    """``GaussHermite.m`` -- nodes/weights for ``int exp(-x^2) f(x)`` on the real line."""
    i = np.arange(1, n)
    a = np.sqrt(i / 2)
    x, V = np.linalg.eigh(np.diag(a, 1) + np.diag(a, -1))
    return x, np.sqrt(np.pi) * V[0, :] ** 2


def gauss_laguerre(n: int, alpha: float = 0.0):
    """``GaussLaguerre.m`` -- nodes/weights for ``int_0^inf x^alpha e^-x f(x)``."""
    from math import gamma

    i = np.arange(1, n + 1)
    a = (2 * i - 1) + alpha
    b = np.sqrt(i[:-1] * (np.arange(1, n) + alpha))
    x, V = np.linalg.eigh(np.diag(a) + np.diag(b, 1) + np.diag(b, -1))
    return x, gamma(alpha + 1) * V[0, :] ** 2


# --------------------------------------------------------------------------
# Adaptive integration (replacement for MATLAB ``integral``)
# --------------------------------------------------------------------------
_GL_X, _GL_W = np.polynomial.legendre.leggauss(15)


def _gl(f, a, b):
    m, r = (a + b) / 2, (b - a) / 2
    return r * np.sum(_GL_W * f(m + r * _GL_X))


def _adaptive(f, a, b, whole, abstol, reltol, depth):
    m = (a + b) / 2
    left, right = _gl(f, a, m), _gl(f, m, b)
    both = left + right
    if depth <= 0 or abs(both - whole) <= max(abstol, reltol * abs(both)):
        return both
    return (_adaptive(f, a, m, left, abstol / 2, reltol, depth - 1)
            + _adaptive(f, m, b, right, abstol / 2, reltol, depth - 1))


def adaptive_integral(f: Callable, a, b, abstol: float = 1e-10,
                      reltol: float = 1e-6, max_depth: int = 30) -> complex:
    """Adaptive Gauss-Legendre integral of a (vectorised, complex) ``f``.

    Mimics MATLAB ``integral(f, a, b)`` defaults (AbsTol 1e-10, RelTol 1e-6).
    ``a``/``b`` may be ``+-inf`` (handled by a rational substitution).  The
    exact adaptive scheme differs from MATLAB's, so results agree to the
    tolerance, not bit-for-bit.
    """
    if np.isinf(a) and np.isinf(b):
        g = lambda u: f(u / (1 - u ** 2)) * (1 + u ** 2) / (1 - u ** 2) ** 2
        lo, hi, sgn = -1.0, 1.0, (1 if b > a else -1)
    elif np.isinf(a):
        g = lambda u: f(b - u / (1 - u)) / (1 - u) ** 2
        lo, hi, sgn = 0.0, 1.0, (1 if a < 0 else -1)
    elif np.isinf(b):
        g = lambda u: f(a + u / (1 - u)) / (1 - u) ** 2
        lo, hi, sgn = 0.0, 1.0, (1 if b > 0 else -1)
    else:
        return _adaptive(f, a, b, _gl(f, a, b), abstol, reltol, max_depth)
    return sgn * _adaptive(g, lo, hi, _gl(g, lo, hi), abstol, reltol, max_depth)


def line_integral(f: Callable, z1: complex, z2: complex, **kw) -> complex:
    """``int f(z) dz`` along the straight segment ``z1 -> z2`` (adaptive)."""
    dz = z2 - z1
    return adaptive_integral(lambda t: f(z1 + t * dz) * dz, 0.0, 1.0, **kw)


def polyline_integral(f: Callable, points: Sequence[complex], **kw) -> complex:
    """Integral along a polyline -- MATLAB ``integral(f, a, b, 'Waypoints', ...)``."""
    return sum(line_integral(f, points[i], points[i + 1], **kw)
               for i in range(len(points) - 1))


def box_contour(x, y):
    """Anticlockwise corner list of the box ``[x1,x2] x [y1,y2]`` (closed).

    This is the path the MATLAB code integrates with
    ``integral(f, x1+iy1, x1+iy1, 'Waypoints', [x2+iy1, x2+iy2, x1+iy2])``.
    """
    return [x[0] + 1j * y[0], x[1] + 1j * y[0], x[1] + 1j * y[1],
            x[0] + 1j * y[1], x[0] + 1j * y[0]]
