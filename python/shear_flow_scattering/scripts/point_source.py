"""Point (vortex-sheet) source in a sheared boundary layer -- no scattering.

Ports of the scripts in ``GenericInversion/PointSource_Cts``:
``PointSource_Generic_Nov18.m``, ``PointSource_Hard_Nov18.m``,
``PointSource_Wake_Feb19.m``, ``PressureZeroTracker.m`` and
``FarField/FarField_{Generic,Wake}PointSource*.m``.

The MATLAB scripts had parameters hard-coded at the top; they are keyword
arguments here with the same defaults.  Each returns a :class:`Struct` with
the grids and the ``I_tot`` field (``phi, dphi, p, v`` over ``x1 x x2``) or
far-field directivity ``FF``.
"""
from __future__ import annotations

import numpy as np

from .._utils import Struct, log
from ..contours import straight_path
from ..inversion import (dispersion_generic_nov2018, k_generic, k_hard, wake_dec_dphi,
                         wake_dec_phi, y_boundary_generic, y_greens_generic, y_greens_hard,
                         y_greens_wake)
from ..profiles import c_const, u_linear, u_parabolic
from ..steepest_descent import far_field_saddle
from ..zeros import newton_raphson
from .common import acoustic_and_cl, add_pole_circles, cl_contours, far_field, sd_contours, wall_operators


def _driver(x1, x2, Yfun, Kfun, Cp, Cm, C_CLp, C_CLm):
    Yp, Ym = Yfun(x2, Cp.k, 0), Yfun(x2, Cm.k, 0)
    Y_CLp, Y_CLm = Yfun(x2, C_CLp.k, 1), Yfun(x2, C_CLm.k, -1)
    Kp, Km = Kfun(Cp.k, 0), Kfun(Cm.k, 0)
    K_CLp, K_CLm = Kfun(C_CLp.k, 1), Kfun(C_CLm.k, -1)
    res = acoustic_and_cl(x1, Cp, Cm, C_CLp, C_CLm, Kp, Km, K_CLp, K_CLm, Yp, Ym, Y_CLp, Y_CLm)
    res.x1, res.x2 = x1, x2
    res.contours = Struct(Cp=Cp, Cm=Cm, C_CLp=C_CLp, C_CLm=C_CLm)
    return res


def point_source_generic(bc="hard", Z=0.1 + 0.1j, U=None, c=None, delta=1.0, omega=1.0,
                         k3=0.0, N_bl=1000, y2=0.5, x1=None, x2=None, N_C=100,
                         angle=3 * np.pi / 8, curvature=0.01, apex=0.1, N_CL=1000,
                         shift1=None, shift2=None, known_zeros=(), radius=0.01, step=0.01):
    """``PointSource_Generic_Nov18.m`` -- source at ``y2`` above a wall with a
    linear BC: ``bc`` is ``"hard"`` (default, as in the file), ``"pressure"``
    (release) or ``"impedance"`` (uses ``Z``).

    ``known_zeros`` are dispersion-relation zeros to encircle (the "Extra, if
    zeros known" block).  The scaling constant ``A`` is "to be fixed" in the
    original and is not applied.
    """
    U = U or u_parabolic(0.3, 1, delta)
    c = c or c_const(5)
    x1 = np.arange(-2, 4 + 1e-9, 0.02) if x1 is None else np.asarray(x1)
    x2 = np.arange(0, 3 + 1e-9, 0.02) if x2 is None else np.asarray(x2)
    shift1 = 0.1 * omega if shift1 is None else shift1
    shift2 = 0.1 * omega if shift2 is None else shift2
    ops = wall_operators(omega, U)
    L = {"hard": ops.hard, "pressure": ops.pressure}.get(bc) or ops.impedance(Z)
    Cp, Cm = sd_contours(N_C, omega, U.f(delta), c.f(delta), k3, apex, curvature, angle)
    C_CLp, C_CLm = cl_contours(N_CL, omega, U.f(delta), U.f(0), shift1, shift2)
    C_CLp, C_CLm = add_pole_circles(C_CLp, C_CLm, known_zeros, radius, step)
    Yfun = lambda xx, k, HP: y_greens_generic(xx, y2, k, k3, omega, U, c, L.l1, L.l2, delta, N_bl, HP)
    Kfun = lambda k, HP: k_generic(y2, k, k3, omega, U, c, L.l1, L.l2, delta, N_bl, HP)
    res = _driver(x1, x2, Yfun, Kfun, Cp, Cm, C_CLp, C_CLm)
    res.kappa1 = omega / U.f(y2)
    return res


def point_source_hard(U=None, c=None, delta=1.0, omega=0.1, k3=0.0, N_bl=1000, y2=0.2,
                      x1=None, x2=None, N_C=100, angle=3 * np.pi / 8, curvature=0.01,
                      apex=None, N_CL=1000, shift1=None, shift2=None):
    """``PointSource_Hard_Nov18.m`` -- source above a hard wall (linear profile)."""
    U = U or u_linear(0.2, 1, delta)
    c = c or c_const(5)
    x1 = np.arange(-2, 5 + 1e-9, 0.01) if x1 is None else np.asarray(x1)
    x2 = np.arange(0, 2 + 1e-9, 0.01) if x2 is None else np.asarray(x2)
    apex = 0.1 * omega if apex is None else apex
    shift1 = 0.4 * omega if shift1 is None else shift1
    shift2 = 0.4 * omega if shift2 is None else shift2
    Cp, Cm = sd_contours(N_C, omega, U.f(delta), c.f(delta), k3, apex, curvature, angle)
    C_CLp, C_CLm = cl_contours(N_CL, omega, U.f(delta), U.f(0), shift1, shift2)
    Yfun = lambda xx, k, HP: y_greens_hard(xx, y2, k, k3, omega, U, c, delta, N_bl, HP)
    Kfun = lambda k, HP: k_hard(y2, k, k3, omega, U, c, delta, N_bl, HP)
    res = _driver(x1, x2, Yfun, Kfun, Cp, Cm, C_CLp, C_CLm)
    res.kappa1 = omega / U.f(y2)
    return res


def _wake_bc(omega, k3, Up, Um, cm, deltam, N_bl):
    """Matching to the lower-side decaying solution, expressed as a wall BC
    ``l1 phi' + l2 phi = 0`` for the upper side (``PointSource_Wake_Feb19``)."""
    C0 = lambda k1: 1j * (omega - Up.f(0) * k1)
    dC0p = lambda k1: -1j * Up.df(0) * k1
    dC0m = lambda k1: -1j * Um.df(0) * k1
    z_below = straight_path(deltam, 0.0)
    phi_b = lambda k1: wake_dec_phi(k1, k3, omega, Um, cm, deltam, N_bl, z_below)
    dphi_b = lambda k1: wake_dec_dphi(k1, k3, omega, Um, cm, deltam, N_bl, z_below)
    l1 = lambda k1: phi_b(k1)
    l2 = lambda k1: dphi_b(k1) + 3 * (dC0m(k1) + dC0p(k1)) / C0(k1) * phi_b(k1)
    return l1, l2


def point_source_wake(U0=0.3, Uinfp=1.0, deltap=1.0, Uinfm=1.0, deltam=1.0, c0=5.0,
                      omega=1.0, k3=0.0, N_bl=1000, y2=0.5, x1=None, x2=None, N_C=200,
                      angle=np.pi / 4, curvature=0.01, apex=0.1, N_CL=1000,
                      shift1=None, shift2=None, zero_guesses=(1.5 + 0.4j, 1.5 - 0.4j),
                      radius=0.01, step=0.01):
    """``PointSource_Wake_Feb19.m`` -- source above a (possibly asymmetric)
    wake with parabolic profiles either side.  The two dispersion zeros are
    found by Newton from ``zero_guesses`` and encircled."""
    Up, Um = u_parabolic(U0, Uinfp, deltap), u_parabolic(U0, Uinfm, deltam)
    cp, cm = c_const(c0), c_const(c0)
    x1 = np.arange(-2, 10 + 1e-9, 0.05) if x1 is None else np.asarray(x1)
    x2 = np.arange(-2, 2 + 1e-9, 0.05) if x2 is None else np.asarray(x2)
    shift1 = 0.1 * omega if shift1 is None else shift1
    shift2 = 0.1 * omega if shift2 is None else shift2
    l1, l2 = _wake_bc(omega, k3, Up, Um, cm, deltam, N_bl)

    Cp, Cm = sd_contours(N_C, omega, Up.f(deltap), cp.f(deltap), k3, apex, curvature, angle)
    C_CLp, C_CLm = cl_contours(N_CL, omega, Up.f(deltap), Up.f(0), shift1, shift2)
    D_l = lambda k1: dispersion_generic_nov2018(k1, k3, omega, Up, cp, l1, l2, deltap, N_bl, 0)
    smol = 1e-12
    dD_l = lambda k1: (D_l(k1 + smol) - D_l(k1 - smol)) / (2 * smol)
    zeds = np.array([newton_raphson(D_l, dD_l, g) for g in zero_guesses])
    log.info("Wake dispersion zeros: %s", zeds)
    C_CLp, C_CLm = add_pole_circles(C_CLp, C_CLm, zeds, radius, step)

    Yfun = lambda xx, k, HP: y_greens_wake(xx, y2, k, k3, omega, Up, cp, Um, cm, l1, l2,
                                           deltap, deltam, N_bl, HP)
    Kfun = lambda k, HP: k_generic(y2, k, k3, omega, Up, cp, l1, l2, deltap, N_bl, HP)
    res = _driver(x1, x2, Yfun, Kfun, Cp, Cm, C_CLp, C_CLm)
    res.dispersion_zeros = zeds
    return res


def pressure_zero_tracker(c0=None, U=None, delta=1.0, omega=1.0, k3=0.0, N_bl=1000,
                          guess=1.4 + 0.3j):
    """``PressureZeroTracker.m`` -- follow the pressure-release dispersion zero
    as the sound speed ``c0`` varies (default 6, 7, ..., 100)."""
    c0 = np.arange(6, 101) if c0 is None else np.asarray(c0)
    U = U or u_linear(0.2, 1, 1)
    P = wall_operators(omega, U).pressure  # BCs don't depend on c
    smol = 1e-14
    z = np.zeros(c0.size, dtype=complex)
    for j, cj in enumerate(c0):
        c = c_const(cj)
        D = lambda k1, c=c: dispersion_generic_nov2018(k1, k3, omega, U, c, P.l1, P.l2, delta, N_bl, 1)
        dD = lambda k1, D=D: (D(k1 + smol) - D(k1 - smol)) / (2 * smol)
        z[j] = newton_raphson(D, dD, guess if j == 0 else z[j - 1])
    return c0, z


def far_field_generic_point_source(bc="pressure", Z=np.exp(-1j * np.pi / 4), U=None, c=None,
                                   delta=1.0, omega=20.0, k3=0.0, N_bl=1000, y2=0.01,
                                   N_theta=1000, plot=False):
    """``FarField_GenericPointSource.m`` -- far-field directivity of a source
    inside the boundary layer (default: pressure-release wall).

    The commented-out parameter sweeps at the end of the MATLAB file (over
    ``c0``, ``omega``, ``y2``) are simply loops over this function.
    """
    U = U or u_parabolic(0.3, 1, delta)
    c = c or c_const(5)
    ops = wall_operators(omega, U)
    L = {"hard": ops.hard, "pressure": ops.pressure}.get(bc) or ops.impedance(Z)
    theta = np.linspace(0, np.pi, N_theta)
    sd = far_field_saddle(theta, omega, U.f(delta), c.f(delta), k3)
    Y = y_boundary_generic(y2, sd.saddle, k3, omega, U, c, L.l1, L.l2, delta, N_bl, 0)
    K = k_generic(y2, sd.saddle, k3, omega, U, c, L.l1, L.l2, delta, N_bl, 0)
    FF = far_field(K, sd.Theta, Y)
    if plot:
        _polar(theta, FF.p)
    return Struct(theta=theta, Theta=sd.Theta, saddle=sd.saddle, FF=FF)


def far_field_wake_point_source(U0=1.0, Uinfp=1.0, deltap=1.0, Uinfm=1.0, deltam=1.0,
                                c0=5.0, omega=10.0, k3=0.0, N_bl=1000, y2=0.5,
                                N_theta_p=200, N_theta_m=200, plot=False):
    """``FarField_WakePointSource_Feb19.m`` -- far field of a source above a wake.

    .. warning:: **Suspicious in the original (kept).**  The lower half
       (``theta in [pi, 2 pi]``) evaluates ``Y_BoundaryGeneric`` and
       ``K_Generic`` with the *upper* flow ``Up, cp, deltap``; only the saddle
       points use the lower flow.  For the symmetric default this makes no
       difference.
    """
    Up, Um = u_parabolic(U0, Uinfp, deltap), u_parabolic(U0, Uinfm, deltam)
    cp, cm = c_const(c0), c_const(c0)
    l1, l2 = _wake_bc(omega, k3, Up, Um, cm, deltam, N_bl)
    theta_p = np.linspace(0, np.pi, N_theta_p)
    theta_m = np.linspace(np.pi, 2 * np.pi, N_theta_m)
    sp = far_field_saddle(theta_p, omega, Up.f(deltap), cp.f(deltap), k3)
    sm = far_field_saddle(theta_m, omega, Um.f(deltam), cm.f(deltam), k3, lower=True)
    Y_p = y_boundary_generic(y2, sp.saddle, k3, omega, Up, cp, l1, l2, deltap, N_bl, 0)
    K_p = k_generic(y2, sp.saddle, k3, omega, Up, cp, l1, l2, deltap, N_bl, 0)
    Y_m = y_boundary_generic(y2, sm.saddle, k3, omega, Up, cp, l1, l2, deltap, N_bl, 0)
    K_m = k_generic(y2, sm.saddle, k3, omega, Up, cp, l1, l2, deltap, N_bl, 0)
    FFp, FFm = far_field(K_p, sp.Theta, Y_p), far_field(K_m, sm.Theta, Y_m)
    FF = Struct(**{f: np.concatenate([FFp[f], FFm[f]]) for f in FFp.fields()})
    theta = np.concatenate([theta_p, theta_m])
    if plot:
        _polar(theta, FF.p)
    return Struct(theta=theta, FF=FF)


def _polar(theta, values):
    import matplotlib.pyplot as plt

    ax = plt.subplot(projection="polar")
    ax.plot(theta, np.abs(values))
    return ax
