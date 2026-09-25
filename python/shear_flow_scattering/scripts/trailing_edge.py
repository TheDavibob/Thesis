"""Trailing-edge scattering: a (lined) plate for ``x1 < 0`` and a wake beyond.

Ports of:

* ``Scattering/TrailingEdge_Lined.m``            -> :func:`trailing_edge_lined`
* ``Scattering/FarField_TrailingEdge_Lined.m``   -> :func:`far_field_trailing_edge_lined`
* ``Scattering/WallPressure_FarField_TrailingEdge_Lined.m``
  -> :func:`far_field_trailing_edge_lined` with ``kappa1=omega/Uc`` and
  ``pressure_upstream=1`` (see :func:`wall_pressure_far_field_trailing_edge_lined`)
* ``Distributed_Source/Distributed_TrailingEdge.m`` -> :func:`distributed_trailing_edge`
* ``Distributed_Source/Distributed_TrailingEdge_VarOmegaC.m``
  -> :func:`rerun_with_vorticity` (re-uses everything that does not depend on
  the vorticity distribution -- the original relied on the MATLAB workspace)
* ``Distributed_Source/WallPressure_PointSource_{Scaled,FarField,FarField_routine}.m``

Geometry: flow ``Up`` above the plate/wake (``x2 > 0``) and ``Um`` below, the
latter written as a function of ``-x2``.  Fields for ``x2 < 0`` are computed in
the reflected coordinate and mapped back with the sign flips of
:data:`.._utils.FIELD_SIGN_BELOW`.  The plate has impedance ``Z``
(``inf`` = hard).
"""
from __future__ import annotations

import numpy as np

from .._utils import FIELD_SIGN_BELOW, FIELDS, Struct, log, row
from ..contours import contour_semicircle, tanh_contour, veitch_contour
from ..distributed import (distributed_wall_pressure, fbar, kappa1_quad, upstream_distributed_a,
                           upstream_distributed_b, upstream_distributed_integrated)
from ..factorisation import fact_mult_cont_scaled_finite
from ..inversion import (dispersion_generic_nov2018, fourier_nov2018, integral_nov2018,
                         y_continuous_scattering, y_greens_generic)
from ..kernels import ff_wake_kernel, wake_kernel
from ..profiles import c_const, u_linear, u_parabolic, u_pw_parabolic
from ..steepest_descent import far_field_saddle
from ..zeros import newton_raphson
from .common import cl_contours, sd_contours


# ---------------------------------------------------------------------------
# Boundary operators for the plate
# ---------------------------------------------------------------------------
def plate_operators(omega, Z, Up, Um):
    """``V1``/``V2`` (velocity), ``P`` (pressure) and the lined-plate
    operators ``L1`` (upper face) / ``L2`` (lower face).

    Sign convention follows ``TrailingEdge_Lined.m`` / ``Wake_Kernel.m``
    (``L2 = -Zed V2 - C0 P``).  ``Distributed_TrailingEdge.m`` instead used
    ``L2 = Zed V2 - C0 P`` with ``DZ2 = -D(L2)``: identical for a hard plate
    (the only case it supported) but inconsistent with ``Wake_Kernel`` for
    finite ``Z``.
    """
    U0 = Up.f(0)
    C0 = lambda k1: 1j * (omega - U0 * k1)
    dC1 = lambda k1: -1j * Up.df(0) * k1
    dC2 = lambda k1: -1j * Um.df(0) * k1
    V1 = Struct(l1=lambda k1: -C0(k1) ** 2, l2=lambda k1: -3 * C0(k1) * dC1(k1))
    V2 = Struct(l1=lambda k1: -C0(k1) ** 2, l2=lambda k1: -3 * C0(k1) * dC2(k1))
    P = Struct(l1=lambda k1: np.zeros(np.shape(k1)), l2=lambda k1: -C0(k1) ** 3)
    if np.isinf(Z):
        L1 = V1
        L2 = Struct(l1=lambda k1: -V2.l1(k1), l2=lambda k1: -V2.l2(k1))
        Zed = np.inf
    else:
        Zed = 1j * omega * Z
        L1 = Struct(l1=lambda k1: Zed * V1.l1(k1) - C0(k1) * P.l1(k1),
                    l2=lambda k1: Zed * V1.l2(k1) - C0(k1) * P.l2(k1))
        L2 = Struct(l1=lambda k1: -Zed * V2.l1(k1) - C0(k1) * P.l1(k1),
                    l2=lambda k1: -Zed * V2.l2(k1) - C0(k1) * P.l2(k1))
    return Struct(V1=V1, V2=V2, P=P, L1=L1, L2=L2, Zed=Zed, C0=C0)


def _dispersions(ops, k3, omega, Up, cp, deltap, Um, cm, deltam, N_bl):
    D = lambda L, U, c, d: (lambda k1, HP: dispersion_generic_nov2018(
        k1, k3, omega, U, c, L.l1, L.l2, d, N_bl, HP))
    DZ1, DZ2 = D(ops.L1, Up, cp, deltap), D(ops.L2, Um, cm, deltam)
    Dh1, Dp1 = D(ops.V1, Up, cp, deltap), D(ops.P, Up, cp, deltap)
    Dh2, Dp2 = D(ops.V2, Um, cm, deltam), D(ops.P, Um, cm, deltam)
    Dw = lambda k1, HP: -Dp1(k1, HP) * Dh2(k1, HP) - Dp2(k1, HP) * Dh1(k1, HP)
    return Struct(DZ1=DZ1, DZ2=DZ2, Dh1=Dh1, Dp1=Dp1, Dh2=Dh2, Dp2=Dp2, Dw=Dw)


def _wake_factorisation(k3, omega, Up, cp, deltap, Um, cm, deltam, Z, N_bl, contour, prec_N,
                        branch_cut):
    z1, z2 = contour_semicircle(deltap, 1), contour_semicircle(deltam, 1)
    Kernel = lambda k1: wake_kernel(k1, k3, omega, Up, cp, deltap, Um, cm, deltam, Z, N_bl, z1, z2)
    Kff = ff_wake_kernel(k3, omega, Up.f(0), cp.f(0), Up.df(0), -Um.df(0), Z, 0)
    return fact_mult_cont_scaled_finite(Kernel, Kff, contour, Struct(N=prec_N), branch_cut)


# ---------------------------------------------------------------------------
# Full-field solver: setup (expensive, forcing-independent) + assembly
# ---------------------------------------------------------------------------
def _te_setup(Up, cp, deltap, Um, cm, deltam, omega, Z, k3, x1, x2, N_bl, N_C, angle,
              curvature, apex, N_CL, shift1, shift2, cl_a, prec_N, branch_cut, mode_guess):
    st = Struct(Up=Up, cp=cp, deltap=deltap, Um=Um, cm=cm, deltam=deltam, omega=omega, Z=Z,
                k3=k3, x1=row(x1), x2=row(x2), N_bl=N_bl)
    st.ops = plate_operators(omega, Z, Up, Um)
    st.up, st.dn = st.x1 < 0, st.x1 >= 0
    st.above, st.below = st.x2 >= 0, st.x2 < 0
    xa, xb = st.x2[st.above], -st.x2[st.below]  # reflected coordinate below

    st.Cp1, st.Cm1 = sd_contours(N_C, omega, Up.f(deltap), cp.f(deltap), k3, apex, curvature, angle)
    st.Cp2, st.Cm2 = sd_contours(N_C, omega, Um.f(deltam), cm.f(deltam), k3, apex, curvature, angle)
    st.C_CLp, st.C_CLm = cl_contours(N_CL, omega, omega / cl_a, Up.f(0), shift1, shift2)
    log.info("Contours set up")

    xu, xd = st.x1[st.up], st.x1[st.dn]
    st.X = Struct(p1=fourier_nov2018(xu, st.Cp1.k), m1=fourier_nov2018(xd, st.Cm1.k),
                  p2=fourier_nov2018(xu, st.Cp2.k), m2=fourier_nov2018(xd, st.Cm2.k),
                  CLp=fourier_nov2018(xd, st.C_CLp.k), CLm=fourier_nov2018(xd, st.C_CLm.k))
    Ysc = lambda xx, k, U, c, d, HP=0: y_continuous_scattering(xx, k, k3, omega, U, c, d, N_bl, HP, HP != 0)
    st.Y = Struct(p1=Ysc(xa, st.Cp1.k, Up, cp, deltap), m1=Ysc(xa, st.Cm1.k, Up, cp, deltap),
                  p2=Ysc(xb, st.Cp2.k, Um, cm, deltam), m2=Ysc(xb, st.Cm2.k, Um, cm, deltam),
                  CLp1=Ysc(xa, st.C_CLp.k, Up, cp, deltap, 1), CLm1=Ysc(xa, st.C_CLm.k, Up, cp, deltap, -1),
                  CLp2=Ysc(xb, st.C_CLp.k, Um, cm, deltam, 1), CLm2=Ysc(xb, st.C_CLm.k, Um, cm, deltam, -1))
    log.info("Computed x2 dependence")

    st.D = D = _dispersions(st.ops, k3, omega, Up, cp, deltap, Um, cm, deltam, N_bl)
    smol = 1e-12
    dDw = lambda k1: (D.Dw(k1 + smol, 0) - D.Dw(k1 - smol, 0)) / (2 * smol)
    # "Assuming, for now, the only modes are in the gLHP and thus only the wake
    # ones are important."
    zg = newton_raphson(lambda k1: D.Dw(k1, 0), dDw, mode_guess)
    st.zeros_wake_m = np.array([zg, np.conj(zg)])
    st.dDw_wake_m = dDw(st.zeros_wake_m)
    log.info("Located wake modes %s", st.zeros_wake_m)

    st.Dc = Struct(DZ1p=D.DZ1(st.Cp1.k, 0), DZ2p=D.DZ2(st.Cp2.k, 0),
                   DZ1m=D.DZ1(st.Cm1.k, 0), DZ2m=D.DZ2(st.Cm2.k, 0),
                   Dw1m=D.Dw(st.Cm1.k, 0), Dw2m=D.Dw(st.Cm2.k, 0),
                   DZ1_CLp=D.DZ1(st.C_CLp.k, 1), DZ1_CLm=D.DZ1(st.C_CLm.k, -1),
                   DZ2_CLp=D.DZ2(st.C_CLp.k, 1), DZ2_CLm=D.DZ2(st.C_CLm.k, -1),
                   Dw_CLp=D.Dw(st.C_CLp.k, 1), Dw_CLm=D.Dw(st.C_CLm.k, -1))
    log.info("Computed dispersion relationships along integration contours")

    Cp, dCp = tanh_contour(0, 1, -1, 1)
    Cm, dCm = tanh_contour(0, 1, -1, 1)
    st.K = K = _wake_factorisation(k3, omega, Up, cp, deltap, Um, cm, deltam, Z, N_bl,
                                   Struct(Cp=Cp, dCp=dCp, Cm=Cm, dCm=dCm), prec_N, branch_cut)
    st.Kf = Struct(p1=K.p(st.Cp1.k), m1=K.m(st.Cm1.k), p2=K.p(st.Cp2.k), m2=K.m(st.Cm2.k),
                   CLp=K.m(st.C_CLp.k), CLm=K.m(st.C_CLm.k), modes=K.m(st.zeros_wake_m))
    log.info("Kernel found and factorised")

    zw = st.zeros_wake_m
    st.mode = Struct(Y1=Ysc(xa, zw, Up, cp, deltap), Y2=Ysc(xb, zw, Um, cm, deltam),
                     X=fourier_nov2018(xd, zw), DZ1=D.DZ1(zw, 0), DZ2=D.DZ2(zw, 0))
    return st


def _place(st, up_rows, blocks_above, blocks_below):
    """Build full ``(x1, x2)`` arrays; ``blocks_*`` map row masks to arrays."""
    out = Struct()
    for f in FIELDS:
        full = np.zeros((st.x1.size, st.x2.size), dtype=complex)
        for rows, blk in blocks_above:
            full[np.ix_(rows, st.above)] += blk[f]
        for rows, blk in blocks_below:
            full[np.ix_(rows, st.below)] += FIELD_SIGN_BELOW[f] * blk[f]
        out[f] = full
    return out


def _te_assemble(st, cauchy, A, A_upstream, I_upstream_above) -> Struct:
    """Combine the precomputed pieces for a given forcing ``cauchy(k)``."""
    Dc, Kf, Y, X = st.Dc, st.Kf, st.Y, st.X
    c_ = lambda C: cauchy(C.k)
    K = Struct(p1=A_upstream * c_(st.Cp1) / (Dc.DZ1p * Kf.p1),
               p2=A_upstream * c_(st.Cp2) / (Dc.DZ2p * Kf.p2),
               m1=A * c_(st.Cm1) * Dc.DZ2m * Kf.m1 / Dc.Dw1m,
               m2=A * c_(st.Cm2) * Dc.DZ1m * Kf.m2 / Dc.Dw2m,
               CLp1=A * c_(st.C_CLp) * Dc.DZ2_CLp * Kf.CLp / Dc.Dw_CLp,
               CLm1=A * c_(st.C_CLm) * Dc.DZ2_CLm * Kf.CLm / Dc.Dw_CLm,
               CLp2=A * c_(st.C_CLp) * Dc.DZ1_CLp * Kf.CLp / Dc.Dw_CLp,
               CLm2=A * c_(st.C_CLm) * Dc.DZ1_CLm * Kf.CLm / Dc.Dw_CLm)
    I = lambda C, k, x, y: Struct(**{f: integral_nov2018(C.wk, k, x, y[f]) for f in FIELDS})

    I_acoustic = _place(st, None,
                        [(st.up, I(st.Cp1, K.p1, X.p1, Y.p1)), (st.dn, I(st.Cm1, K.m1, X.m1, Y.m1))],
                        [(st.up, I(st.Cp2, K.p2, X.p2, Y.p2)), (st.dn, I(st.Cm2, K.m2, X.m2, Y.m2))])
    I_crit = _place(st, None,
                    [(st.dn, I(st.C_CLp, K.CLp1, X.CLp, Y.CLp1)), (st.dn, I(st.C_CLm, K.CLm1, X.CLm, Y.CLm1))],
                    [(st.dn, I(st.C_CLp, K.CLp2, X.CLp, Y.CLp2)), (st.dn, I(st.C_CLm, K.CLm2, X.CLm, Y.CLm2))])

    zw, m = st.zeros_wake_m, st.mode
    common = -2j * np.pi * A * Kf.modes * cauchy(zw) / st.dDw_wake_m
    K1w, K2w = common * m.DZ2, common * m.DZ1
    modal_above = Struct(**{f: sum(K1w[k] * np.outer(m.X[:, k], m.Y1[f][:, k]) for k in range(zw.size))
                            for f in FIELDS})
    modal_below = Struct(**{f: sum(K2w[k] * np.outer(m.X[:, k], m.Y2[f][:, k]) for k in range(zw.size))
                            for f in FIELDS})
    I_modal = _place(st, None, [(st.dn, modal_above)], [(st.dn, modal_below)])

    I_upstream = _place(st, None, [(np.ones(st.x1.size, bool), I_upstream_above)], [])
    I_scattered = Struct(**{f: I_modal[f] + I_crit[f] + I_acoustic[f] for f in FIELDS})
    I_total = Struct(**{f: I_scattered[f] + I_upstream[f] for f in FIELDS})
    return Struct(x1=st.x1, x2=st.x2, I_acoustic=I_acoustic, I_crit=I_crit, I_modal=I_modal,
                  I_upstream=I_upstream, I_scattered=I_scattered, I_total=I_total,
                  zeros_wake_m=zw, K=st.K)


def _leading(I0, Km_kappa1, Z, Zed):
    """Leading coefficients ``(A, A_upstream)`` for a plate of impedance ``Z``."""
    if np.isinf(Z):
        A = I0 / (1j * Km_kappa1)
        return A, A
    if Z == 0:
        raise NotImplementedError(
            "Z = 0 is not usable in the original code: A = I0/(1i*Zed*Km) divides by Zed = 0, "
            "and the velocity forcing (Velocity_Upstream) is undefined in the far-field scripts.")
    A = I0 / (1j * Zed * Km_kappa1)
    return A, Zed * A


def trailing_edge_lined(Up=None, Um=None, cp=None, cm=None, deltap=1.0, deltam=1.0,
                        omega=5.0, Z=2 - 2j, k3=0.0, y2=0.5, x1=None, x2=None, N_bl=2000,
                        N_C=500, angle=np.pi / 4, curvature=0.01, apex=None, N_CL=1000,
                        shift1=None, shift2=None, Q0=1e-10, prec_N=2000,
                        branch_cut=3 * np.pi / 2, mode_guess=0.9 + 0.35j):
    """``TrailingEdge_Lined.m`` -- vortex sheet at ``y2`` above a lined plate
    ending at ``x1 = 0`` (asymmetric linear profiles by default).

    Returns ``I_total`` (full field on the ``x1 x x2`` grid) together with
    its acoustic, critical-layer, modal and incident parts.  (Very expensive
    at the default resolution.)
    """
    Up = Up or u_linear(0.3, 1, deltap)
    Um = Um or u_linear(0.3, 1.5, deltam)
    cp, cm = cp or c_const(5), cm or c_const(5)
    x1 = np.arange(-10, 10 + 1e-9, 0.05) if x1 is None else np.asarray(x1)
    x2 = np.arange(-10, 10 + 1e-9, 0.05) if x2 is None else np.asarray(x2)
    apex = 0.1 * omega if apex is None else apex
    shift1 = 0.1 * omega if shift1 is None else shift1
    shift2 = 0.1 * omega if shift2 is None else shift2
    kappa1 = omega / Up.f(y2)
    cl_a = min(omega / Up.f(deltap), omega / Um.f(deltam))  # longest critical layer
    st = _te_setup(Up, cp, deltap, Um, cm, deltam, omega, Z, k3, x1, x2, N_bl, N_C, angle,
                   curvature, apex, N_CL, shift1, shift2, cl_a, prec_N, branch_cut, mode_guess)
    ops = st.ops

    log.info("Computing upstream behaviour")
    xa = st.x2[st.above]
    L1 = ops.L1
    Y_up = y_greens_generic(xa, y2, kappa1, k3, omega, Up, cp, L1.l1, L1.l2, deltap, 10 * N_bl, -1)
    wall = y_greens_generic(0, y2, kappa1, k3, omega, Up, cp, L1.l1, L1.l2, deltap, 10 * N_bl, -1)
    Xu = np.exp(-1j * st.x1 * kappa1)
    I_up = Struct(**{f: Q0 * np.outer(Xu, Y_up[f][:, 0]) for f in FIELDS})
    at_wall = lambda L: Q0 * (L.l1(kappa1) * wall.dphi[0, 0] + L.l2(kappa1) * wall.phi[0, 0])
    I0 = at_wall(ops.P) if Z != 0 else at_wall(ops.V1)

    A, A_up = _leading(I0, st.K.m(kappa1), Z, ops.Zed)
    res = _te_assemble(st, lambda k: 1 / (k - kappa1), A, A_up, I_up)
    res.kappa1 = kappa1
    return res


def far_field_trailing_edge_lined(U0=0.3, Uinfp=1.0, Uinfm=1.0, deltap=1.0, deltam=1.0, c0=5.0,
                                  omega=5.0, Z=np.inf, k3=0.0, q0=1.0, y2=0.1, N_bl=1000,
                                  N_theta_p=300, N_theta_m=300, prec_N=10000, scale=1.0,
                                  branch_cut=3 * np.pi / 2, shift_y2=0.001, Up=None, Um=None,
                                  cp=None, cm=None, theta1=None, theta2=None, kappa1=None,
                                  pressure_upstream=None):
    """``FarField_TrailingEdge_Lined.m`` -- far-field directivity of the
    trailing-edge problem ("Only really optimised for hard-wall").

    .. warning:: **Missing input in the original.**  ``U0`` was commented out
       (``%U0 = 0.3;``) yet used to build the profiles, so the script only ran
       with ``U0`` already in the workspace.  It defaults to 0.3 here (the
       value in the commented line).

    ``kappa1`` may be an array (one forcing wavenumber per source height, as in
    ``WallPressure_PointSource_FarField.m``); results then gain a trailing
    axis.  ``pressure_upstream`` overrides the computed wall-pressure forcing.
    """
    Up = Up or u_parabolic(U0, Uinfp, deltap)
    Um = Um or u_parabolic(U0, Uinfm, deltam)
    cp, cm = cp or c_const(c0), cm or c_const(c0)
    theta1 = np.linspace(0, np.pi, N_theta_p) if theta1 is None else np.atleast_1d(theta1)
    theta2 = np.linspace(np.pi, 2 * np.pi, N_theta_m) if theta2 is None else np.atleast_1d(theta2)
    y2 = np.atleast_1d(y2)
    kappa1 = omega / Up.f(y2) if kappa1 is None else np.atleast_1d(kappa1)
    s1 = far_field_saddle(theta1, omega, Up.f(deltap), cp.f(deltap), k3)
    s2 = far_field_saddle(theta2, omega, Um.f(deltam), cm.f(deltam), k3, lower=True)
    S1, S2 = s1.saddle, s2.saddle
    ops = plate_operators(omega, Z, Up, Um)

    shift = 0.5 * (0.5 * (S1[0] - S1[-1]))
    Cp, dCp = veitch_contour(S1[0] - shift, -scale * omega, 0.5 * omega / Up.f(0))
    Cm, dCm = veitch_contour(S1[-1] + shift, -scale * omega, 0.5 * omega / Up.f(0))
    K = _wake_factorisation(k3, omega, Up, cp, deltap, Um, cm, deltam, Z, N_bl,
                            Struct(Cp=Cp, dCp=dCp, Cm=Cm, dCm=dCm), prec_N, branch_cut)
    Km_kappa1 = K.m(kappa1)

    if pressure_upstream is None:
        wp = distributed_wall_pressure(y2, omega / Up.f(y2), shift_y2, omega, k3, Up, cp, deltap,
                                       N_bl, ops.L1.l1, ops.L1.l2, -1)
        pressure_upstream = -q0 * wp * cp.f(y2) ** 2 / Up.df(y2)
    if Z == 0:
        raise NotImplementedError("Velocity_Upstream is undefined in the original far-field script")
    A, A_up = _leading(np.asarray(pressure_upstream), Km_kappa1, Z, ops.Zed)

    D = _dispersions(ops, k3, omega, Up, cp, deltap, Um, cm, deltam, N_bl)
    col = lambda v: np.asarray(v)[:, None]
    C1 = 1 / (col(S1) - kappa1[None, :])
    C2 = 1 / (col(S2) - kappa1[None, :])
    K1_p = A_up * C1 / col(D.DZ1(S1, 0) * K.p(S1))
    K2_p = A_up * C2 / col(D.DZ2(S2, 0) * K.p(S2))
    K1_m = A * C1 * col(D.DZ2(S1, 0) * K.m(S1) / D.Dw(S1, 0))
    K2_m = A * C2 * col(D.DZ1(S2, 0) * K.m(S2) / D.Dw(S2, 0))

    Y1 = y_continuous_scattering(deltap, S1, k3, omega, Up, cp, deltap, 10)
    Y2 = y_continuous_scattering(deltam, S2, k3, omega, Um, cm, deltam, 10)
    ff = lambda Kx, Th, Y: Struct(**{f: np.squeeze(Kx * col(np.sin(Th)) * col(Y[f][0]))
                                     for f in ("phi", "p", "v")})
    FF1p, FF1m = ff(K1_p, s1.Theta, Y1), ff(K1_m, s1.Theta, Y1)
    FF2p, FF2m = ff(K2_p, s2.Theta, Y2), ff(K2_m, s2.Theta, Y2)
    theta = np.concatenate([theta1, theta2])
    use_m = (theta < np.pi / 2) | (theta > 3 * np.pi / 2)
    FF = Struct()
    for f in ("phi", "p", "v"):
        fp = np.concatenate([np.atleast_1d(FF1p[f]).reshape(theta1.size, -1),
                             np.atleast_1d(FF2p[f]).reshape(theta2.size, -1)])
        fm = np.concatenate([np.atleast_1d(FF1m[f]).reshape(theta1.size, -1),
                             np.atleast_1d(FF2m[f]).reshape(theta2.size, -1)])
        fp[use_m] = fm[use_m]
        FF[f] = np.squeeze(fp)
    return Struct(theta=theta, FF=FF, FF1p=FF1p, FF1m=FF1m, FF2p=FF2p, FF2m=FF2m, K=K)


def wall_pressure_far_field_trailing_edge_lined(Uc_ratio=0.8, Uinfp=1.0, omega=5.0, **kw):
    """``WallPressure_FarField_TrailingEdge_Lined.m`` -- as
    :func:`far_field_trailing_edge_lined` but forced by a unit wall pressure
    convecting at ``Uc = 0.8 Uinf``."""
    return far_field_trailing_edge_lined(Uinfp=Uinfp, omega=omega,
                                         kappa1=omega / (Uc_ratio * Uinfp),
                                         pressure_upstream=1.0, **kw)


# ---------------------------------------------------------------------------
# Distributed vorticity
# ---------------------------------------------------------------------------
def distributed_trailing_edge(Up=None, Um=None, cp=None, cm=None, deltap=1.0, deltam=1.0,
                              omega=1.0, Z=np.inf, k3=0.0, x1=None, x2=None, N_bl=2000,
                              N_distribution=80, shift_y2=0.01, omega_c=None, N_C=500,
                              angle=3 * np.pi / 8, curvature=0.01, apex=None, N_CL=1000,
                              shift1=None, shift2=None, prec_N=2000, branch_cut=np.pi,
                              mode_guess=0.9 + 0.35j):
    """``Distributed_TrailingEdge.m`` -- hard plate with a vorticity
    distribution ``omega_c(y2)`` across the upper boundary layer.

    Returns the result Struct; ``res.state`` holds everything independent of
    ``omega_c`` so :func:`rerun_with_vorticity` can cheaply try other
    distributions (the ``_VarOmegaC`` script).

    .. warning:: **Bug fixed.**  The critical-layer start was
       ``min(omega/Up.f(deltap), omega/Up.f(deltap))`` -- ``Up`` twice; the
       lined version (and the comment) use ``Um`` for the second term.
    """
    Up = Up or u_parabolic(0.3, 1, deltap)
    Um = Um or u_parabolic(0.3, 1.5, deltam)
    cp, cm = cp or c_const(5), cm or c_const(5)
    x1 = np.arange(-10, 10 + 1e-9, 0.05) if x1 is None else np.asarray(x1)
    x2 = np.arange(-2, 2 + 1e-9, 0.02) if x2 is None else np.asarray(x2)
    apex = 0.1 * omega if apex is None else apex
    shift1 = 0.1 * omega if shift1 is None else shift1
    shift2 = 0.1 * omega if shift2 is None else shift2
    if not np.isinf(Z):
        log.warning("Distributed_TrailingEdge was only written for a hard plate (Z = inf)")

    kappa1, dkappa1, y2 = kappa1_quad(omega, Up, deltap, N_distribution)
    cl_a = min(omega / Up.f(deltap), omega / Um.f(deltam))
    st = _te_setup(Up, cp, deltap, Um, cm, deltam, omega, Z, k3, x1, x2, N_bl, N_C, angle,
                   curvature, apex, N_CL, shift1, shift2, cl_a, prec_N, branch_cut, mode_guess)
    L1 = st.ops.L1
    log.info("Computing distributed upstream behaviour")
    G_A = upstream_distributed_a(y2, kappa1, shift_y2, omega, k3, Up, cp, deltap, 10 * N_bl,
                                 L1.l1, L1.l2, 1)
    st.G_01 = upstream_distributed_b(st.x2[st.above], G_A, omega, k3, Up, cp, deltap)
    st.I0 = distributed_wall_pressure(y2, kappa1, shift_y2, omega, k3, Up, cp, deltap, N_bl,
                                      L1.l1, L1.l2, 1)
    st.Km_distributed = st.K.m(kappa1)
    st.update(kappa1=kappa1, dkappa1=dkappa1, y2=y2)
    return rerun_with_vorticity(st, np.ones(y2.size) if omega_c is None else omega_c)


def rerun_with_vorticity(state: Struct, omega_c) -> Struct:
    """``Distributed_TrailingEdge_VarOmegaC.m`` -- recompute only what depends
    on the vorticity distribution ``omega_c`` (length ``N_distribution``).

    Example from the original: ``omega_c = y2 * (deltap - y2)``.
    """
    st = state
    omega_c = row(omega_c)
    I_up = upstream_distributed_integrated(st.x1, omega_c, st.G_01, st.dkappa1)
    Int = (-st.omega * st.cp.f(st.y2) ** 2 * st.I0 * omega_c
           / (st.Km_distributed * st.Up.f(st.y2) ** 2))
    F_bar = lambda k: fbar(k, st.kappa1, st.dkappa1, Int)
    A = 1 / 1j  # "Much simpler than before"
    res = _te_assemble(st, F_bar, A, A, I_up)
    res.update(state=st, omega_c=omega_c, y2=st.y2, kappa1=st.kappa1)
    return res


# ---------------------------------------------------------------------------
# Wall-pressure studies (base vs sheared boundary layer)
# ---------------------------------------------------------------------------
def _study_profiles():
    base = u_pw_parabolic([0.6, 0.7, 1], [0, 0.1, 1])
    sheared = u_pw_parabolic([0.2, 0.75, 1], [0, 0.3, 1])
    return base, sheared


def wall_pressure_point_source_scaled(y2=None, omega=None, shift_y2=0.005, cinf=35.0,
                                      N_bl=3000, k3=0.0, q0=1.0, delta=1.0, plot=False):
    """``WallPressure_PointSource_Scaled.m`` -- wall pressure beneath a point
    vortex at each ``y2``, for a "base" and a "sheared" piecewise-parabolic
    boundary layer, over a range of frequencies.

    ``f = 10 omega / (2 pi 0.02)`` is the dimensional frequency used for
    plotting.  The original called ``PlotRI`` (plot real/imag parts), which
    is **missing from the repository**; ``plot=True`` draws an equivalent.
    """
    y2 = np.arange(0.1, 0.5 + 1e-9, 0.1) if y2 is None else np.asarray(y2)
    omega = np.logspace(0, 3, 30) if omega is None else np.asarray(omega)
    c0 = c_const(cinf)
    base, sheared = _study_profiles()
    out = {"base": np.zeros((omega.size, y2.size), complex),
           "sheared": np.zeros((omega.size, y2.size), complex)}
    for j, Om in enumerate(omega):
        for name, U in (("base", base), ("sheared", sheared)):
            C0 = lambda k1: 1j * (Om - U.f(0) * k1)
            dC1 = lambda k1: -1j * U.df(0) * k1
            l1 = lambda k1: -C0(k1) ** 2
            l2 = lambda k1: -3 * C0(k1) * dC1(k1)
            wp = distributed_wall_pressure(y2, Om / U.f(y2), shift_y2, Om, k3, U, c0, delta,
                                           N_bl, l1, l2, -1)
            out[name][j] = -q0 * wp * c0.f(y2) ** 2 / U.df(y2)
        log.info("%d of %d", j + 1, omega.size)
        if plot:
            _plot_ri(y2, np.vstack([out["base"][j], out["sheared"][j]]))
    return Struct(y2=y2, omega=omega, f=10 * omega / (2 * np.pi * 0.02),
                  Full_Pressure_base=out["base"], Full_Pressure_sheared=out["sheared"])


def wall_pressure_point_source_far_field(Omega, U, y2, cinf=35.0, delta=1.0, k3=0.0,
                                         N_bl=3000, prec_N=1000, Z=np.inf):
    """``WallPressure_PointSource_FarField.m`` -- far field directly above
    (``theta = pi/2``) and below (``-pi/2``) the trailing edge for a unit
    wall-pressure source convecting at ``U(y2)``.  Returns ``FF1p.p`` (one
    value per ``y2``), which is what the ``_routine`` collects."""
    c0 = c_const(cinf)
    y2 = np.atleast_1d(y2)
    r = far_field_trailing_edge_lined(Up=U, Um=U, cp=c0, cm=c0, deltap=delta, deltam=delta,
                                      omega=Omega, Z=Z, k3=k3, N_bl=N_bl, prec_N=prec_N,
                                      theta1=[np.pi / 2], theta2=[-np.pi / 2], y2=y2,
                                      kappa1=Omega / U.f(y2), pressure_upstream=1.0)
    return np.ravel(r.FF1p.p)


def wall_pressure_point_source_far_field_routine(y2=None, omega=None, cinf=35.0, N_bl=3000):
    """``WallPressure_PointSource_FarField_routine.m`` -- loop the above over
    frequency for the base and sheared profiles."""
    y2 = np.linspace(0, 1, 40) if y2 is None else np.asarray(y2)
    omega = np.logspace(-1, 3, 40) if omega is None else np.asarray(omega)
    base, sheared = _study_profiles()
    FFb = np.zeros((omega.size, y2.size), complex)
    FFs = np.zeros((omega.size, y2.size), complex)
    for j, Om in enumerate(omega):
        FFb[j] = wall_pressure_point_source_far_field(Om, base, y2, cinf, N_bl=N_bl)
        FFs[j] = wall_pressure_point_source_far_field(Om, sheared, y2, cinf, N_bl=N_bl)
        log.info("%d of %d", j + 1, omega.size)
    return Struct(y2=y2, omega=omega, f=10 * omega / (2 * np.pi * 0.02),
                  FF_Pressure_base=FFb, FF_Pressure_sheared=FFs)


def _plot_ri(x, Y):
    """Stand-in for the missing ``PlotRI``: real (solid) and imaginary
    (dashed) parts of each row of ``Y`` against ``x``."""
    import matplotlib.pyplot as plt

    plt.figure()
    for j, y in enumerate(np.atleast_2d(Y)):
        plt.plot(x, y.real, "-", color=f"C{j}")
        plt.plot(x, y.imag, "--", color=f"C{j}")
