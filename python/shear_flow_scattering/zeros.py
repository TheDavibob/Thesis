"""Locating zeros and poles of analytic functions in the complex plane.

Ports of the ``ZP finder`` directory and ``CommonTools/PrincipalValue``.

Two approaches coexist in the original:

* **argument principle** -- ``(1/2 pi i) oint f'/f dz`` counts zeros minus
  poles in a box; ``oint z^j f'/f`` gives power sums of the zeros, from which
  up to two zeros can be recovered in closed form;
* **Newton-Raphson** from initial guesses, with deflation
  ``f -> f / (z - z0)`` so the same zero is not found twice.

Most routines take the function *and* its derivative (``fun``, ``dfun``);
:func:`numerical_derivatives` / :func:`pointwise_derivative` build the latter.

.. note::
   MATLAB anonymous functions capture variables *by value* when created;
   Python closures capture *by reference*.  Deflation loops that rebind
   ``F`` therefore pass the previous function in explicitly here.
"""
from __future__ import annotations

from typing import Callable, List, Sequence

import numpy as np

from ._utils import Struct, log, trapezoid
from .quadrature import box_contour, gauss_legendre, lgwt, polyline_integral


# ==========================================================================
# Basic tools
# ==========================================================================
def newton_raphson(f: Callable, df: Callable, z0, tol: float = 1e-10, max_iter: int = 1000):
    """``NewtonRaphson.m`` -- scalar Newton iteration to ``|f| <= 1e-10``.

    Returns NaN if ``|f|`` ever exceeds its starting value (the original's
    "BREAK CLAUSE NEEDS CHECKING").  ``max_iter`` is new: the MATLAB loop
    could spin forever on a stagnating iteration.
    """
    val = lambda z: complex(np.ravel(f(z))[0])
    Z = complex(z0)
    f0 = abs(val(Z))
    for _ in range(max_iter):
        fz = val(Z)
        if abs(fz) <= tol:
            return Z
        if abs(fz) > f0:
            return complex(np.nan, np.nan)
        Z = Z - fz / complex(np.ravel(df(Z))[0])
    log.warning("newton_raphson: no convergence after %d iterations", max_iter)
    return Z


def pointwise_derivative(f: Callable, z, h: float = 1e-10):
    """``PointwiseDerivative.m`` -- central difference with real step ``h``."""
    return (f(z + h) - f(z - h)) / (2 * h)


def numerical_derivatives(f: Callable, acc: float) -> Struct:
    """``Numerical_Derivatives.m`` -- struct with ``f``, ``df``, ``d2f``."""
    return Struct(f=f, df=lambda z: (f(z + acc) - f(z - acc)) / (2 * acc),
                  d2f=lambda z: (f(z + acc) - 2 * f(z) + f(z - acc)) / acc ** 2)


def _box(fun, x, y):
    return polyline_integral(fun, box_contour(x, y))


def box_integral(fun: Callable, x, y):
    """``BoxIntegral.m`` -- ``(1/2 pi i) oint fun`` round the box (sum of residues)."""
    return _box(lambda z: fun(z) / (2j * np.pi), x, y)


def zeros_poles(fun: Callable, dfun: Callable, x, y) -> float:
    """``ZerosPoles.m`` -- (#zeros - #poles) in ``x1<Re z<x2, y1<Im z<y2``."""
    return float(np.real(_box(lambda z: dfun(z) / (2j * np.pi * fun(z)), x, y)))


def sum_zeros_j(fun: Callable, dfun: Callable, x, y, j: int):
    """``SumZerosj.m`` -- ``sum_n z_n^j`` over the zeros inside the box."""
    return _box(lambda z: z ** j * dfun(z) / (2j * np.pi * fun(z)), x, y)


def _two_zeros(fun, dfun, x, y):
    d1 = sum_zeros_j(fun, dfun, x, y, 1)
    d2 = sum_zeros_j(fun, dfun, x, y, 2)
    r = np.sqrt(complex(2 * d2 - d1 ** 2))
    return [(d1 + r) / 2, (d1 - r) / 2]


def simple_zero_finder(fun, dfun, x, y):
    """``SimpleZeroFinder.m`` -- up to two zeros of a pole-free ``fun`` in a box."""
    if abs(box_integral(fun, x, y)) > 1e-3:
        raise ValueError("the function contains singularities within the box")
    Z = zeros_poles(fun, dfun, x, y)
    if Z > 2.001:
        raise ValueError("more than two zeros found within the box - try moving or shrinking the box")
    if Z > 1.001:
        return np.array(_two_zeros(fun, dfun, x, y))
    return sum_zeros_j(fun, dfun, x, y, 1)


def complete_zero_finder(fun, dfun, x, y, plot=False, min_size=1e-5):
    """``CompleteZeroFinder.m`` -- all zeros of a pole-free ``fun`` in a box,
    by quartering any sub-box holding more than two zeros.

    Boxes shrunk below ``min_size`` are treated as a cluster of multiple
    zeros at their corner (original "break clause").
    """
    if abs(box_integral(fun, x, y)) > 0.1:
        raise ValueError("the function contains singularities within the box")
    boxes = [(tuple(x), tuple(y))]
    Z = [zeros_poles(fun, dfun, x, y)]
    n_zeros = int(round(Z[0]))
    log.info("Number of zeros in region %d", n_zeros)
    while max(Z) > 2.001:
        I = int(np.argmax(Z))
        (x1, x2), (y1, y2) = boxes[I]
        xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
        new = [((x1, xm), (y1, ym)), ((x1, xm), (ym, y2)),
               ((xm, x2), (y1, ym)), ((xm, x2), (ym, y2))]
        newZ = []
        for bx, by in new:
            zz = zeros_poles(fun, dfun, bx, by)
            if min(abs(bx[1] - bx[0]), abs(by[1] - by[0])) < min_size:
                zz = -1
            newZ.append(zz)
        boxes[I:I + 1] = new
        Z[I:I + 1] = newZ
    z: List[complex] = []
    for (bx, by), zz in zip(boxes, Z):
        if zz == -1:
            mult = int(round(zeros_poles(fun, dfun, bx, by)))
            z += [bx[0] + 1j * by[0]] * mult
        elif zz < 0.001:
            continue
        elif zz < 1.001:
            z.append(sum_zeros_j(fun, dfun, bx, by, 1))
        else:
            z += _two_zeros(fun, dfun, bx, by)
    z = np.array(z)
    if plot:
        import matplotlib.pyplot as plt

        plt.scatter(z.real, z.imag, marker="o")
    return z


# ==========================================================================
# Gauss-Legendre box integrals (ZP finder/Zeros)
# ==========================================================================
def complex_integral_line(fun, z1, z2, N):
    """``ComplexIntegral_Line.m`` -- N-point Gauss-Legendre ``int fun dz``, z1 -> z2."""
    s, w = gauss_legendre(N)
    z = (z1 + z2) / 2 + s * (z2 - z1) / 2
    return np.sum(fun(z) * w * (z2 - z1) / 2)


def complex_integral_box(fun, x, y, Nh, Nv):
    """``ComplexIntegral_Box.m`` -- anticlockwise round a box; ``Nh``/``Nv``
    points on horizontal/vertical sides."""
    if np.size(x) != 2 or np.size(y) != 2:
        raise ValueError("Box must be given by two x and y coordinates")
    X, Y = np.sort(x), np.sort(y)
    Z = X[:, None] + 1j * Y[None, :]
    return (complex_integral_line(fun, Z[0, 0], Z[1, 0], Nh)
            + complex_integral_line(fun, Z[1, 0], Z[1, 1], Nv)
            + complex_integral_line(fun, Z[1, 1], Z[0, 1], Nh)
            + complex_integral_line(fun, Z[0, 1], Z[0, 0], Nv))


def number_zeros_box(fun: Struct, x, y, Nh, Nv, precision):
    """``NumberZeros_Box.m`` -- zero count of analytic ``fun.f`` (with
    ``fun.df``) in a box, after a Cauchy-theorem analyticity check."""
    if abs(complex_integral_box(fun.f, x, y, Nh, Nv)) > precision:
        raise ValueError("function is not analytic in the box")
    N2 = complex_integral_box(lambda z: fun.df(z) / (2j * np.pi * fun.f(z)), x, y, Nh, Nv)
    if abs(N2 - round(N2.real)) > precision:
        raise ValueError("function seems to have a fractional number of zeros: either accuracy "
                         "is insufficient or branch cuts have been crossed")
    return int(round(N2.real))


def number_zeros_bisection(fun: Struct, x, y, Nh, Nv, recurse, precision):
    """``NumberZeros_Bisection.m`` -- zero counts on a ``2^recurse`` grid of
    sub-boxes, refined only where zeros exist.

    Returns ``(N, Z, sz, X, Y)``: counts per box, box centres, and ``sz``
    (for single-zero boxes, the zero itself via ``oint z f'/f``; if any box
    holds several zeros ``sz`` is just the centres).
    """
    fine = 2 ** recurse
    x, y = np.sort(x), np.sort(y)
    X = np.linspace(x[0], x[1], fine + 1)
    Y = np.linspace(y[0], y[1], fine + 1)
    Z2 = X[:, None] + 1j * Y[None, :]
    Z = 0.5 * (Z2[:-1, :-1] + Z2[1:, 1:])
    N = np.full((fine, fine), number_zeros_box(fun, x, y, Nh, Nv, precision))
    for j in range(1, recurse + 1):
        step = fine // 2 ** j
        idx = np.arange(0, fine + 1, step)
        for k in range(len(idx) - 1):
            for l in range(len(idx) - 1):
                if N[idx[k], idx[l]] != 0:
                    N[idx[k]:idx[k + 1], idx[l]:idx[l + 1]] = number_zeros_box(
                        fun, [X[idx[k]], X[idx[k + 1]]], [Y[idx[l]], Y[idx[l + 1]]], Nh, Nv, precision)
        log.info("Completed bisection %d", j)
    if N.max() > 1:
        log.info("Some boxes have multiple zeros")
        return N, Z, Z, X, Y
    zdff = lambda z: z * fun.df(z) / (2j * np.pi * fun.f(z))
    sz = np.zeros(N.shape, dtype=complex)
    for k in range(fine):
        for l in range(fine):
            if N[k, l] != 0:
                sz[k, l] = complex_integral_box(zdff, [X[k], X[k + 1]], [Y[l], Y[l + 1]], Nh, Nv)
    return N, Z, sz, X, Y


def zeros_box(f: Callable, x, y, N=500, prec=1e-2) -> Struct:
    """``Zeros_Box.m`` -- counts zeros of ``f`` in a box: total (``z``), on
    the real line (``r``), upper (``U``) and lower (``L``) half-planes.

    .. warning:: **Bug fixed.**  ``N = varargin(1)`` in the 4-argument call
       produced a cell, not a number.  Also: counts for half-planes not in the
       box are ``0`` here (MATLAB used ``[]``).
    """
    prec2 = prec * 1e-5
    x1, x2 = min(x), max(x)
    y1, y2 = min(y), max(y)
    df = lambda z: (f(z + prec2) - f(z - prec2)) / (2 * prec2)
    g = lambda z: df(z) / (2j * np.pi * f(z))
    t, wt = lgwt(N, 0, 1)
    t, wt = t[::-1], wt[::-1]

    def count(ll, lu, ul, uu):
        side = lambda a, b: np.sum(g(a + (b - a) * t) * (b - a) * wt)
        return np.round(side(ll, ul) + side(ul, uu) - side(lu, uu) - side(ll, lu), 1)

    Z = Struct(z=count(x1 + 1j * y1, x1 + 1j * y2, x2 + 1j * y1, x2 + 1j * y2), r=0, U=0, L=0)
    log.info("%s zeros in box", Z.z)
    if y1 < 0 < y2:
        Z.r = count(x1 - 1j * prec, x1 + 1j * prec, x2 - 1j * prec, x2 + 1j * prec)
    if y2 > 0:
        Z.U = count(x1 + 1j * prec, x1 + 1j * y2, x2 + 1j * prec, x2 + 1j * y2)
    if y1 < 0:
        Z.L = count(x1 + 1j * y1, x1 - 1j * prec, x2 + 1j * y1, x2 - 1j * prec)
    if Z.r + Z.U + Z.L != Z.z:
        log.warning("Zeros don't add up? Possible precision failure")
    return Z


# ==========================================================================
# Newton with deflation
# ==========================================================================
def _multiplicity(f, df, z0, N_int=100, rad=1e-3):
    t, wt = lgwt(N_int, 0, 2 * np.pi)
    t, wt = t[::-1], wt[::-1]
    s = rad * np.exp(1j * t)
    return np.round(np.real(np.sum(wt * s * df(z0 + s) / (2 * np.pi * f(z0 + s)))), 1)


def _deflate(f, df, z0, M):
    """Divide out ``(k - z0)^M`` and return the *correct* derivative.

    .. warning:: **Bug fixed.**  The MATLAB rebuilt ``F.df`` *after*
       reassigning ``F.f``, so ``F.df = F.df/(k-z) - F.f/(k-z)^2`` used the
       already-deflated ``f`` (wrong by a factor ``1/(k-z)^M`` in the second
       term).  Newton still often converged, but slowly.
    """
    g = lambda k: f(k) / (k - z0) ** M
    dg = lambda k: df(k) / (k - z0) ** M - M * f(k) / (k - z0) ** (M + 1)
    return g, dg


def zero_finder_recursive(f: Callable, z0: Sequence[complex], Z: Struct = None):
    """``ZeroFinder_Recursive.m`` -- Newton from each guess, deflating as it goes.

    Guesses are split into real / upper / lower half-plane groups; ``Z`` (e.g.
    from :func:`zeros_box`) optionally gives the expected counts.  Returns a
    Struct with ``r``, ``U``, ``L`` and ``z`` (all of them).

    .. warning:: **Bug fixed.**  The lower-half-plane check was
       ``imag(z) < prec_real`` -- true for *every* LHP zero, so all LHP zeros
       were rejected and the routine errored after 10 nudges.  It now rejects
       ``imag(z) > -prec_real``.
    """
    z0 = np.asarray(z0, dtype=complex)
    init = Struct(r=z0[z0.imag == 0], U=z0[z0.imag > 0], L=z0[z0.imag < 0])
    if Z is None:
        Z = Struct(z=z0.size, r=init.r.size, U=init.U.size, L=init.L.size)
    for name, what in (("r", "real"), ("U", "positive imaginary"), ("L", "negative imaginary")):
        if init[name].size != Z[name]:
            log.warning("guess has different number of %s zeros to the function", what)
            break
    if z0.size < Z.z:
        log.warning("initial guess lower than the total number of zeros, some will be missed")
    elif z0.size > Z.z:
        log.info("More elements in initial guess than zeros of the function: ignoring latter guesses")

    prec_diff, prec_real = 1e-10, 1e-5
    F = f
    dF = lambda k: (f(k + prec_diff) - f(k - prec_diff)) / (2 * prec_diff)
    n_r = min(Z.r, z0.size)
    n_U = min(Z.U, z0.size - n_r)
    n_L = min(Z.L, z0.size - n_r - n_U)
    out = Struct()
    checks = {"r": lambda z: abs(z.imag) > prec_real,
              "U": lambda z: z.imag < prec_real,
              "L": lambda z: z.imag > -prec_real}
    for name, n in (("r", n_r), ("U", n_U), ("L", n_L)):
        guesses = np.sort_complex(init[name])
        if n > guesses.size and guesses.size:
            guesses = np.sort_complex(np.concatenate([guesses, np.repeat(guesses[-1], n - guesses.size)]))
        guesses = list(guesses[:n])
        found = np.zeros(n, dtype=complex)
        j, count = 0, 0
        while j < n:
            zj = newton_raphson(F, dF, guesses[j])
            bad = np.isnan(zj) or checks[name](zj)
            if bad:
                guesses[j] *= 1.1  # "give it a nudge"
                count += 1
                if count > 10:
                    raise RuntimeError(f"zero finding failed for group '{name}'")
                continue
            M = int(_multiplicity(F, dF, zj))
            if M < 1:
                raise RuntimeError("Something fishy about the multiplicity")
            found[j:j + M] = zj
            F, dF = _deflate(F, dF, zj, M)
            j += M
            count = 0
            log.info("Found %d zeros in group '%s'", min(j, n), name)
        out[name] = found
    out.z = np.concatenate([out.r, out.U, out.L])
    return out


def poles(f: Callable, p0: Sequence[complex]):
    """``Poles.m`` -- poles of ``f`` as zeros of ``1/f``, Newton from guesses
    ``p0`` with deflation.  Guesses must be good."""
    p0 = np.asarray(p0, dtype=complex)
    p = np.zeros(p0.shape, dtype=complex)
    F = lambda z: 1 / f(z)
    for j in range(p0.size):
        if j > 0:
            F = (lambda Fold, pj: (lambda z: Fold(z) / (z - pj)))(F, p[j - 1])
        p[j] = newton_raphson(F, (lambda Fj: (lambda z: pointwise_derivative(Fj, z)))(F), p0[j])
    return p


def simple_zero_tracker(f: Callable, p: Sequence, z0, plot=False):
    """``SimpleZeroTracker.m`` -- follow a zero of ``f(z, p)`` as ``p`` varies,
    using each result as the next guess.  Stops (NaN-filling) on failure."""
    p = np.asarray(p)
    z = np.full(p.shape, np.nan, dtype=complex)
    for j, pj in enumerate(p):
        F = (lambda pj: (lambda zz: f(zz, pj)))(pj)
        z[j] = newton_raphson(F, lambda zz, F=F: pointwise_derivative(F, zz),
                              z0 if j == 0 else z[j - 1])
        if np.isnan(z[j]):
            break
    if plot:
        import matplotlib.pyplot as plt

        plt.plot(z.real, z.imag)
    return z


# ==========================================================================
# Pointwise (non-vectorised) integrals
# ==========================================================================
def pointwise_trapezium(f: Callable, x):
    """``PointwiseTrapezium.m`` -- trapezium rule for a scalar-only ``f``."""
    x = np.asarray(x)
    return trapezoid(np.array([f(xi) for xi in x.ravel()]), x.ravel())


def pointwise_trapezium_box(f: Callable, x, y):
    """``PointwiseTrapeziumBox.m`` -- anticlockwise round ``z = x + iy``."""
    x, y = np.asarray(x), np.asarray(y)
    return (pointwise_trapezium(f, x + 1j * y[0]) + pointwise_trapezium(f, x[-1] + 1j * y)
            - pointwise_trapezium(f, x + 1j * y[-1]) - pointwise_trapezium(f, x[0] + 1j * y))


def pointwise_zp(f: Callable, x, y):
    """``PointwiseZP.m`` -- zeros minus poles using the pointwise tools."""
    return pointwise_trapezium_box(
        lambda z: pointwise_derivative(f, z) / (2j * np.pi * f(z)), x, y)


def pv_trapz(f: Callable, a, b, p, epsilon, N=1000):
    """``PVTrapz.m`` -- principal-value integral of ``f`` over ``[a, b]``,
    excising ``+-epsilon`` round each pole in ``p`` (brute-force trapezium)."""
    p = np.sort(np.atleast_1d(p))
    starts = np.concatenate([[a], p + epsilon])
    ends = np.concatenate([p - epsilon, [b]])
    total = 0
    for s, e in zip(starts, ends):
        r = np.linspace(s, e, N)
        total += trapezoid(f(r), r)
    return total


# ==========================================================================
# Branch cuts -- incomplete in the original
# ==========================================================================
def branch_cut_integral(f: Callable, B, delta):
    """``BranchCutIntegral.m`` -- trapezium integrals along both sides of a
    polyline ``B`` offset by ``+-delta`` normals.  Returns ``(I1, I2)``;
    ``I1 + I2 ~ 0`` if there is no cut."""
    B = np.asarray(B, dtype=complex)
    nrm = -1j * (B[1:] - B[:-1])
    nrm = np.append(nrm, nrm[-1])
    nrm = nrm / np.abs(nrm)
    C2 = B + delta * nrm
    C1 = B - delta * nrm
    if C1[0].imag < C1[-1].imag:
        C1 = C1[::-1]
    elif C1[0].imag > C1[-1].imag:
        C2 = C2[::-1]
    elif C1[0].real > C1[-1].real:
        C1 = C1[::-1]
    else:
        C2 = C2[::-1]
    return trapezoid(f(C1), C1), trapezoid(f(C2), C2)


def zeros_poles_branch_cuts(fun, dfun, x, y):
    """``ZerosPolesBranchCuts.m`` -- **cannot run: dependencies missing.**

    Needs ``ContourSweep`` (to detect branch cuts crossing the box) and
    ``BranchCuts`` (to assemble them into polylines), neither of which is in
    the repository.  The original also contains ``numel(B,1)``, which is not
    valid MATLAB for this purpose (``size(B,1)`` was presumably meant).
    """
    raise NotImplementedError(
        "ZerosPolesBranchCuts depends on ContourSweep.m and BranchCuts.m, which "
        "are missing from the repository.")
