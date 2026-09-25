"""Pieces repeated (copy-pasted) across the MATLAB driver scripts."""
from __future__ import annotations

import numpy as np

from .._utils import FIELDS, Struct, row
from ..contours import append_contour, finite_cl_ellipse, finite_symmetric_hyperbola, pole_circle
from ..inversion import fourier_nov2018, integral_nov2018
from ..pridmore_brown import branch_points


def wall_operators(omega, U, sigma=None):
    """Boundary operators ``L phi = l1 phi' + l2 phi`` at ``x2 = 0``.

    Returns a Struct of ``(l1, l2)`` pairs: ``hard`` (``v = 0``),
    ``pressure`` (``p = 0``, "pressure release") and ``impedance(Z)``
    (``i omega Z v - C p = 0``).  ``sigma`` overrides ``U'(0)``.
    """
    U0 = U.f(0)
    s0 = U.df(0) if sigma is None else sigma
    C0 = lambda k1: 1j * (omega - U0 * k1)
    dC0 = lambda k1: -1j * s0 * k1
    V = Struct(l1=lambda k1: -C0(k1) ** 2, l2=lambda k1: -3 * C0(k1) * dC0(k1))
    P = Struct(l1=lambda k1: np.zeros(np.shape(k1)), l2=lambda k1: -C0(k1) ** 3)

    def impedance(Z):
        if np.isinf(Z):
            return V
        Zed = 1j * omega * Z
        return Struct(l1=lambda k1: Zed * V.l1(k1) - C0(k1) * P.l1(k1),
                      l2=lambda k1: Zed * V.l2(k1) - C0(k1) * P.l2(k1))

    return Struct(hard=V, pressure=P, impedance=impedance, C0=C0, dC0=dC0)


def sd_contours(N_C, omega, Uinf, cinf, k3, apex, curvature, angle):
    """The hyperbolic "steepest-descent-like" pair: ``Cp`` round ``kb2``
    (``x1 < 0``) and ``Cm`` round ``kb1`` (``x1 >= 0``)."""
    kb1, kb2 = branch_points(omega, Uinf, cinf, k3)
    Cp = finite_symmetric_hyperbola(N_C, kb2, -apex, curvature, angle)
    Cm = finite_symmetric_hyperbola(N_C, kb1, apex, curvature, angle)
    return Cp, Cm


def cl_contours(N_CL, omega, Ubl_edge, U0, shift1, shift2, fudge=100.0):
    """Ellipses round the critical layer ``[omega/U(delta), omega/U(0)]``.

    For a no-slip wall (``U0 = 0``) the cut extends to infinity; the original
    truncates at ``100 a`` ("a fudge, but useful" / "this crosses the
    critical layer").
    """
    a = omega / Ubl_edge
    b = omega / U0 if U0 != 0 else fudge * a
    return finite_cl_ellipse(N_CL, a, b, shift1, shift2)


def add_pole_circles(C_CLp, C_CLm, zeros, radius=0.01, step=0.01):
    """Append small clockwise circles round known dispersion zeros:
    upper-half-plane zeros to ``C_CLp``, the rest to ``C_CLm``."""
    for z in np.atleast_1d(zeros):
        circ = pole_circle(z, radius, step)
        if np.imag(z) > 0:
            C_CLp = append_contour(C_CLp, circ)
        else:
            C_CLm = append_contour(C_CLm, circ)
    return C_CLp, C_CLm


def upstream_field(x1, Y: Struct, kappa1, Q0) -> Struct:
    """Convected incident field ``Q0 exp(-i kappa1 x1) Y(x2)`` on the full grid."""
    X = np.exp(-1j * row(x1) * kappa1)
    return Struct(**{f: Q0 * np.outer(X, Y[f][:, 0]) for f in FIELDS})


def acoustic_and_cl(x1, Cp, Cm, C_CLp, C_CLm, Kp, Km, K_CLp, K_CLm, Yp, Ym, Y_CLp, Y_CLm,
                    x1_shift=0.0) -> Struct:
    """The standard assembly used by every single-junction script::

        I     = [Ip; Im]                 (x1 < 0 rows from Cp, x1 >= 0 from Cm)
        I_CL  = I_CLp + I_CLm            (x1 >= 0 only)
        I_tot = [Ip; Im + I_CL]

    ``x1_shift`` translates the junction to ``x1 = x1_shift``.  Rows are
    placed by mask, so ``x1`` need not be sorted.
    """
    xs = row(x1) - x1_shift
    up, dn = xs < 0, xs >= 0
    Xp, Xm = fourier_nov2018(xs[up], Cp.k), fourier_nov2018(xs[dn], Cm.k)
    X_CLp, X_CLm = fourier_nov2018(xs[dn], C_CLp.k), fourier_nov2018(xs[dn], C_CLm.k)
    n2 = Yp.phi.shape[0]
    res = Struct(I=Struct(), I_CL=Struct(), I_tot=Struct())
    for f in FIELDS:
        I = np.zeros((xs.size, n2), dtype=complex)
        I[up] = integral_nov2018(Cp.wk, Kp, Xp, Yp[f])
        I[dn] = integral_nov2018(Cm.wk, Km, Xm, Ym[f])
        I_CL = np.zeros_like(I)
        I_CL[dn] = (integral_nov2018(C_CLp.wk, K_CLp, X_CLp, Y_CLp[f])
                    + integral_nov2018(C_CLm.wk, K_CLm, X_CLm, Y_CLm[f]))
        res.I[f], res.I_CL[f], res.I_tot[f] = I, I_CL, I + I_CL
    return res


def add_fields(*structs) -> Struct:
    """Field-wise sum of several ``phi/dphi/p/v`` structs."""
    return Struct(**{f: sum(s[f] for s in structs) for f in FIELDS})


def far_field(K, Theta, Y: Struct) -> Struct:
    """Saddle-point directivity ``K sin(Theta) Y`` for each field."""
    return Struct(**{f: K * np.sin(Theta) * Y[f] for f in FIELDS})


def blend_far_field(FFp: Struct, FFm: Struct, use_minus) -> Struct:
    """Use the ``+`` evaluation upstream and the ``-`` one where ``use_minus``
    (downstream, ``theta < pi/2``) -- the two agree between the contours."""
    out = Struct()
    for f in FIELDS:
        v = np.array(FFp[f], dtype=complex)
        v[use_minus] = FFm[f][use_minus]
        out[f] = v
    return out
