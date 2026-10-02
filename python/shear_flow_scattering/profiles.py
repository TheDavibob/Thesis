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
    """Struct with ``f``, ``df``, ``d2f`` (optionally ``d3f``) callables.

    ``analytic`` records whether ``f`` is an analytic function of complex
    ``x2``.  The deformed (complex) integration paths used to step round the
    critical layer are only valid for analytic profiles; for piecewise ones
    (``breaks`` lists the joins, each piece analytic) the solvers switch to
    real-axis paths with a small local indentation round the critical point.
    """

    def __init__(self, f: Callable, df: Callable, d2f: Optional[Callable] = None,
                 d3f: Optional[Callable] = None, analytic: bool = True, breaks=()):
        super().__init__(f=f, df=df, analytic=analytic, breaks=tuple(breaks))
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


def u_power(U0, Uinf, deltainf, n) -> Profile:
    """``U = Uinf - (Uinf - U0) (1 - x2/delta)^n`` (new; not in the MATLAB).

    Analytic, monotonic and inflection-free, with ``U'(delta) = 0`` for
    ``n >= 2``; ``n = 2`` is :func:`u_parabolic`.  Wall shear
    ``U'(0) = n (Uinf - U0) / delta``, so ``n`` varies the near-wall shear at
    fixed slip velocity.
    """
    D = Uinf - U0
    xi = lambda x2: 1 - np.asarray(x2) / deltainf
    return Profile(lambda x2: Uinf - D * xi(x2) ** n,
                   lambda x2: n * D / deltainf * xi(x2) ** (n - 1),
                   lambda x2: -n * (n - 1) * D / deltainf ** 2 * xi(x2) ** (n - 2))


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
    with continuous value and shear; ``U'(deltainf) = 0``.

    .. warning:: Not analytic (``U''`` jumps at the join, and the branch is
       chosen by ``Re x2``), so results from complex/semicircular integration
       paths are wrong -- by O(1) for the paper's Fig. 5 "base" profile.
       The solvers therefore replace complex paths by real-axis paths with a
       local indentation (:func:`.contours.indented_path`); an analytic fit
       (:func:`fit_analytic`) is the alternative.
    """
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
        analytic=False, breaks=(d1,),  # piecewise: see warning above
        f=lambda x2: pick(x2, lambda x: U1 + sigma * (x - d1) + A * (x - d1) ** 2,
                        lambda x: U1 + sigma * (x - d1) + B * (x - d1) ** 2),
        df=lambda x2: pick(x2, lambda x: sigma + 2 * A * (x - d1),
                           lambda x: sigma + 2 * B * (x - d1)),
        d2f=lambda x2: pick(x2, lambda x: 2 * A * _ones(x), lambda x: 2 * B * _ones(x)))


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


def fit_analytic(target: Profile, U0=None, Uinf=1.0, deltainf=1.0, degree=8, deriv_weight=0.3,
                 n_fit=400) -> Profile:
    """Analytic polynomial fit to a (possibly piecewise) profile.

    ``U = Uinf - (Uinf - U0) sum_{j=2}^{degree} b_j xi^j`` with
    ``xi = 1 - x2/deltainf`` and ``sum b_j = 1``, so ``U(0) = U0``,
    ``U(delta) = Uinf`` and ``U'(delta) = 0`` hold exactly.  The ``b_j`` are a
    least-squares fit to ``target.f`` and ``target.df`` (weight
    ``deriv_weight``) on ``[0, deltainf]``.

    Intended for the two-piece profiles of the paper's Fig. 5, which are not
    analytic and so cannot be used with complex integration paths.  The
    result records ``fit_rms``, ``fit_max`` (errors in U), ``fit_max_dU`` and
    whether it is ``monotonic`` / ``inflection_free`` on ``[0, delta]``.
    """
    U0 = float(np.real(target.f(0.0))) if U0 is None else U0
    D = Uinf - U0
    x = np.linspace(0, deltainf, n_fit)
    xi = 1 - x / deltainf
    J = np.arange(3, degree + 1)
    # b_2 = 1 - sum_{j>=3} b_j  =>  U = Uinf - D [xi^2 + sum_j b_j (xi^j - xi^2)]
    rhs_U = (Uinf - np.real(target.f(x))) / D - xi ** 2
    A_U = xi[:, None] ** J[None, :] - xi[:, None] ** 2
    # dU/dx = (D/delta) [2 xi + sum_j b_j (j xi^(j-1) - 2 xi)]
    rhs_dU = np.real(target.df(x)) * deltainf / D - 2 * xi
    A_dU = J[None, :] * xi[:, None] ** (J[None, :] - 1) - 2 * xi[:, None]
    A = np.vstack([A_U, deriv_weight * A_dU])
    rhs = np.concatenate([rhs_U, deriv_weight * rhs_dU])
    bj = np.linalg.lstsq(A, rhs, rcond=None)[0]
    b = np.concatenate([[1 - bj.sum()], bj])            # b_2 ... b_degree
    powers = np.arange(2, degree + 1)

    # Horner evaluation (much faster than explicit powers inside the RK loop)
    from numpy.polynomial import polynomial as npoly
    c_f = np.concatenate([[0.0, 0.0], b])                  # coefficients of xi^0 .. xi^degree
    c_d = npoly.polyder(c_f)                               # d/dxi
    c_dd = npoly.polyder(c_f, 2)

    def f(x2):
        return Uinf - D * npoly.polyval(1 - np.asarray(x2) / deltainf, c_f)

    def df(x2):   # d/dx2 = -(1/delta) d/dxi
        return D / deltainf * npoly.polyval(1 - np.asarray(x2) / deltainf, c_d)

    def d2f(x2):
        return -D / deltainf ** 2 * npoly.polyval(1 - np.asarray(x2) / deltainf, c_dd)

    P = Profile(f, df, d2f)
    xx = np.linspace(0, deltainf, 2001)
    err = np.real(f(xx)) - np.real(target.f(xx))
    P.update(coeffs=b, fit_rms=float(np.sqrt(np.mean(err ** 2))), fit_max=float(np.max(np.abs(err))),
             fit_max_dU=float(np.max(np.abs(np.real(df(xx)) - np.real(target.df(xx))))),
             monotonic=bool(np.all(np.real(df(xx)) >= -1e-12)),
             inflection_free=bool(np.all(np.real(d2f(xx)) <= 1e-9)))
    return P


def afshari_profiles():
    """The paper's Fig. 5 two-piece profiles (``Struct_U_PWParabolic``):
    ``(base, sheared)``, inspired by Afshari et al. (2016)."""
    return (u_pw_parabolic([0.6, 0.7, 1], [0, 0.1, 1]),
            u_pw_parabolic([0.2, 0.75, 1], [0, 0.3, 1]))
