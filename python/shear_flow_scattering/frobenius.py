"""Frobenius series solutions of ``z^2 y'' + alpha(z) z y' + beta(z) y = 0``.

Port of ``CommonTools/Two_Solutions_Complex.m`` and ``param_setup.m``.

.. warning:: **Partly reconstructed.**  ``Two_Solutions_Complex.m`` calls
   three helpers that are *not in the repository*:

   * ``TaylorCoeff_Complex(f, N, r, m)`` -- Taylor coefficients of ``f``;
   * ``Two_Series(N, a, b)`` -- the Frobenius recurrences;
   * ``Taylor_Series(c)`` -- evaluate ``sum c_n z^n``.

   They are re-implemented below from their call signatures and from how
   their outputs are used (e.g. ``dy2`` fixes the form
   ``y2 = C log(z) y1 + z^s2 sum d_n z^n``).  The maths is standard, but the
   exact conventions of the lost originals -- in particular the meaning of the
   integer flag ``I`` and the normalisation ``d_m`` for integer root
   separation -- are guesses.  Nothing else in the repository calls this
   function.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from ._utils import Struct
from .complex_funcs import log_a, power_a


def param_setup() -> Struct:
    """``param_setup.m`` -- default series parameters."""
    return Struct(N=20, m=100, r=1.0, cut=np.pi)


def taylor_coeff_complex(f: Callable, N: int, r: float, m: int) -> np.ndarray:
    """*(Reconstructed)* ``TaylorCoeff_Complex`` -- coefficients ``c_0..c_N`` of
    ``f(z) = sum c_n z^n`` from ``m`` samples on ``|z| = r`` (FFT / trapezium
    rule on the Cauchy integral).  Needs ``m > N``; ``r`` inside the radius
    of convergence."""
    theta = 2 * np.pi * np.arange(m) / m
    samples = np.asarray(f(r * np.exp(1j * theta)), dtype=complex) * np.ones(m)
    return (np.fft.fft(samples) / m)[: N + 1] / r ** np.arange(N + 1)


def taylor_series(c) -> Callable:
    """*(Reconstructed)* ``Taylor_Series`` -- ``z -> sum_n c_n z^n``."""
    c = np.asarray(c, dtype=complex)
    return lambda z: np.polynomial.polynomial.polyval(z, c)


def two_series(N: int, a, b):
    """*(Reconstructed)* ``Two_Series`` -- Frobenius coefficients.

    With ``F(s) = s(s-1) + a0 s + b0`` (indicial polynomial), roots ordered
    so ``Re s1 >= Re s2`` and ``m = s1 - s2``:

    * ``y1 = z^s1 sum c_n z^n``, ``c_0 = 1``;
    * ``y2 = C log(z) y1 + z^s2 sum d_n z^n``.

    Returns ``(c, d, C, s1, s2, I)`` where the (guessed) type flag ``I`` is
    1 = roots not differing by an integer (``C = 0``), 2 = repeated root
    (``C = 1``, ``d_0 = 0``), 3 = roots differing by a positive integer
    (``d_0 = 1``, ``C`` possibly 0, ``d_m = 0``).
    """
    a = np.asarray(a, dtype=complex)
    b = np.asarray(b, dtype=complex)
    F = lambda s: s * (s - 1) + a[0] * s + b[0]
    roots = np.roots([1, a[0] - 1, b[0]])
    s1, s2 = sorted(roots, key=lambda s: (s.real, s.imag), reverse=True)
    m = s1 - s2

    def rhs(coef, s, n):
        return sum(((s + n - k) * a[k] + b[k]) * coef[n - k] for k in range(1, n + 1))

    c = np.zeros(N + 1, dtype=complex)
    c[0] = 1
    for n in range(1, N + 1):
        c[n] = -rhs(c, s1, n) / F(s1 + n)

    # g_j: coefficient of z^(j+s1) in L[y1 log z] - log z L[y1]
    g = np.array([(2 * (j + s1) - 1) * c[j] + sum(a[k] * c[j - k] for k in range(j + 1))
                  for j in range(N + 1)])
    d = np.zeros(N + 1, dtype=complex)
    tol = 1e-6  # a, b come from sampled Taylor coefficients, so not exact
    m_int = int(round(m.real))
    if abs(m) < tol:  # repeated root
        I, C = 2, 1.0
        for n in range(1, N + 1):
            d[n] = -(rhs(d, s2, n) + C * g[n]) / F(s2 + n)
    elif abs(m - m_int) < tol and m_int > 0:  # integer separation
        I, C = 3, 0.0
        d[0] = 1
        for n in range(1, N + 1):
            if n < m_int:
                d[n] = -rhs(d, s2, n) / F(s2 + n)
            elif n == m_int:
                C = -rhs(d, s2, n) / g[0]  # g[0] = s1 - s2 != 0
                d[n] = 0  # free constant (adds a multiple of y1)
            else:
                d[n] = -(rhs(d, s2, n) + C * g[n - m_int]) / F(s2 + n)
    else:
        I, C = 1, 0.0
        d[0] = 1
        for n in range(1, N + 1):
            d[n] = -rhs(d, s2, n) / F(s2 + n)
    return c, d, C, s1, s2, I


def two_solutions_complex(alpha: Callable, beta: Callable, param: Struct = None) -> Struct:
    """``Two_Solutions_Complex.m`` -- two independent solutions as callables.

    Returns ``y`` with ``y1``, ``dy1``, ``y2``, ``dy2`` (functions of ``z``)
    and ``s1``, ``s2``, ``C``, ``c``, ``d``, ``I``, ``a``, ``b``.
    ``param`` has ``N`` (terms), ``r`` (sampling radius), ``m`` (samples,
    ``m >> N``) and ``cut`` (branch-cut angle for logs/powers).
    """
    param = param or param_setup()
    a = taylor_coeff_complex(alpha, param.N, param.r, param.m)
    b = taylor_coeff_complex(beta, param.N, param.r, param.m)
    c, d, C, s1, s2, I = two_series(param.N, a, b)
    # "Clean up the functions a little bit" (original thresholds)
    C = 0 if abs(C) < 1e-5 else C
    s1 = 0 if abs(s1) < 1e-5 else s1
    s2 = 0 if abs(s2) < 1e-5 else s2
    n = np.arange(param.N + 1)
    f1, f2 = taylor_series(c), taylor_series(d)
    nf1, nf2 = taylor_series(c * n), taylor_series(d * n)
    cut = param.cut
    y = Struct(s1=s1, s2=s2, C=C, c=c, d=d, I=I, a=a, b=b)
    y.y1 = lambda z: power_a(z, s1, cut) * f1(z)
    y.dy1 = lambda z: power_a(z, s1, cut) / z * (s1 * f1(z) + nf1(z))
    y.y2 = lambda z: C * log_a(z, cut) * y.y1(z) + power_a(z, s2, cut) * f2(z)
    y.dy2 = lambda z: (C * y.y1(z) / z + C * y.dy1(z) * log_a(z, cut)
                       + power_a(z, s2, cut) / z * (s2 * f2(z) + nf2(z)))
    return y
