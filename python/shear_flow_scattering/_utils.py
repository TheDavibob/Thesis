"""Small shared helpers that replace MATLAB idioms.

MATLAB structs are used everywhere in the original code (``U.f``, ``dec.phi``,
``C.k``/``C.wk`` ...).  They are replaced here by :class:`Struct`, a thin
``SimpleNamespace`` that also supports ``s["field"]`` access (the MATLAB
``s.(name)`` dynamic-field idiom used in the later scattering scripts).
"""
from __future__ import annotations

import logging
from types import SimpleNamespace
from typing import Callable, Iterable

import numpy as np

log = logging.getLogger("shear_flow_scattering")

#: The four physical variables carried through every inversion script:
#: velocity potential-like ``phi``, its wall-normal derivative ``dphi``,
#: pressure ``p`` and wall-normal velocity ``v``.
FIELDS = ("phi", "dphi", "p", "v")

#: Sign flip applied to each of :data:`FIELDS` when a solution computed in
#: ``-x2 > 0`` is mapped back to ``x2 < 0`` (odd quantities change sign).
FIELD_SIGN_BELOW = {"phi": 1, "dphi": -1, "p": 1, "v": -1}


class Struct(SimpleNamespace):
    """MATLAB-struct replacement: attribute *and* item access."""

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def fields(self):
        return list(vars(self))

    def update(self, **kw):
        vars(self).update(kw)
        return self

    def copy(self) -> "Struct":
        return Struct(**vars(self))


class Path(Struct):
    """A parametrised integration path ``x2 = f(s)``, ``s`` in ``[0, 1]``.

    Equivalent to the MATLAB ``z.f`` / ``z.df`` structs passed to the
    ``Comp_*`` Pridmore-Brown routines.
    """

    def __init__(self, f: Callable, df: Callable):
        super().__init__(f=f, df=df)


def heaviside(x):
    """``heaviside.m``: ``x >= 0`` (on the real part, as MATLAB compares), as 0./1.

    MATLAB logicals act as 0/1 in arithmetic, scalar or array.  (The first
    port used ``(np.real(x) >= 0).astype(float)``, which fails for a Python
    scalar: the comparison then returns a plain ``bool``.)
    """
    return np.where(np.real(x) >= 0, 1.0, 0.0)


def trapezoid(y, x):
    """``trapz(x, y)`` for complex data (``np.trapz`` renamed in numpy 2)."""
    y = np.asarray(y)
    x = np.asarray(x)
    return np.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1]))


def row(a) -> np.ndarray:
    """MATLAB ``a(:).'`` -- flatten to a 1-D array (complex-safe)."""
    return np.ravel(np.asarray(a))


def progress(iterable: Iterable, desc: str = "", every: int = 0):
    """Stand-in for MATLAB ``waitbar``.

    Logs at DEBUG level; set ``every`` to log only every n-th item.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    items = list(iterable)
    n = len(items)
    for i, item in enumerate(items, 1):
        if every and (i % every == 0 or i == n):
            log.debug("%s: %d/%d", desc, i, n)
        yield item


def pchip(x, y, xi):
    """Shape-preserving piecewise cubic Hermite interpolation.

    Replacement for ``interp1(x, y, xi, 'pchip')`` without SciPy.  ``y`` may
    be complex (real and imaginary parts are interpolated separately, as
    MATLAB does).  ``x`` may be increasing or decreasing.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y)
    xi = np.asarray(xi, dtype=float)
    if np.iscomplexobj(y):
        return pchip(x, y.real, xi) + 1j * pchip(x, y.imag, xi)
    if x[0] > x[-1]:
        x, y = x[::-1], y[::-1]
    h = np.diff(x)
    delta = np.diff(y) / h
    n = len(x)
    d = np.zeros(n)
    if n == 2:
        d[:] = delta[0]
    else:
        # Fritsch-Carlson interior slopes (same scheme as MATLAB/SciPy)
        for k in range(1, n - 1):
            if delta[k - 1] * delta[k] > 0:
                w1 = 2 * h[k] + h[k - 1]
                w2 = h[k] + 2 * h[k - 1]
                d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])
        # one-sided, shape-preserving end slopes
        def _end(h0, h1, d0, d1):
            s = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
            if np.sign(s) != np.sign(d0):
                s = 0.0
            elif np.sign(d0) != np.sign(d1) and abs(s) > abs(3 * d0):
                s = 3 * d0
            return s
        d[0] = _end(h[0], h[1], delta[0], delta[1])
        d[-1] = _end(h[-1], h[-2], delta[-1], delta[-2])
    k = np.clip(np.searchsorted(x, xi) - 1, 0, n - 2)
    t = (xi - x[k]) / h[k]
    h00 = (1 + 2 * t) * (1 - t) ** 2
    h10 = t * (1 - t) ** 2
    h01 = t ** 2 * (3 - 2 * t)
    h11 = t ** 2 * (t - 1)
    return h00 * y[k] + h10 * h[k] * d[k] + h01 * y[k + 1] + h11 * h[k] * d[k + 1]
