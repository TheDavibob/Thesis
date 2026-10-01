"""Scattering of a convected vortex-sheet disturbance by wall-BC junctions.

Ports of the ``GenericInversion/Scattering`` scripts:

* ``HardSoft_Continuous_Nov2018.m``      -> :func:`hard_soft_continuous`
* ``FarField_HardSoft_Continuous_Nov2018.m`` -> :func:`far_field_hard_soft`
* ``HardSoftHard_Continuous_Feb2019.m``  -> :func:`hard_soft_hard`
* ``FarField_HardSoftHard_Continuous_Feb2019.m`` -> :func:`far_field_hard_soft_hard`
* ``ImpImp_Continuous_Jan2019.m``        -> :func:`imp_imp_continuous`
* ``FarField_ImpImp_Continuous_Feb2019.m`` -> :func:`far_field_imp_imp`

Common recipe: an incident disturbance convecting at ``kappa1 = omega/U(y2)``
(from the point-source Green's function upstream) meets a change in wall
condition at ``x1 = 0``.  The Wiener-Hopf kernel (ratio of dispersion
relations) is factorised numerically, and the scattered field is inverted
along the hyperbolic contours (acoustic part), round the critical layer, and
-- where instability modes exist -- via residues.

Incident field: every script takes ``incident="regularised"`` (default; a
unit sheet via :func:`..inversion.y_greens_convected` -- finite, resolution
independent, continuous across the sheet, causal branch) or ``"original"``
(the MATLAB construction, with the arbitrary amplitude ``Q0 = 1e-10``, which
"scales out" only because ``phi`` is evaluated at its own critical level and
grows like ``N_bl^3``).  Use ``"original"`` to reproduce the thesis figures;
see :func:`.common.incident_sheet`.
"""
from __future__ import annotations

import numpy as np

from .._utils import FIELDS, Struct, log, row
from ..contours import contour_semicircle, straight_path, tanh_contour
from ..factorisation import fact_mult_cont_scaled_finite
from ..inversion import dispersion_generic_nov2018, fourier_nov2018, y_continuous_scattering
from ..kernels import (ff_kernel_hardsoft_noslip, ff_kernel_hardsoft_slip, ff_kernel_impedance,
                       impedance_up_down, kernel_generic, kernel_hardsoft)
from ..pridmore_brown import comp_dec_aug18a_skip
from ..profiles import c_const, u_parabolic
from ..steepest_descent import far_field_saddle
from ..zeros import newton_raphson
from .common import (acoustic_and_cl, add_fields, blend_far_field, cl_contours, far_field,
                     incident_sheet, sd_contours, wall_operators)


# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------
def _contour(shift=0.0, height=1.0):
    Cp, dCp = tanh_contour(0, height, -1 + shift, 1 + shift)
    Cm, dCm = tanh_contour(0, height, -1 - shift, 1 - shift)
    return Struct(Cp=Cp, dCp=dCp, Cm=Cm, dCm=dCm)


def _hardsoft_kernel(k3, omega, U, c, delta, N_bl, contour, prec_N, branch_cut):
    """Factorised hard -> soft kernel ``D_p / D_h`` with the slip / no-slip far field."""
    z_kernel = contour_semicircle(delta, 1)  # so the branch cut is in the LHP
    Kernel = lambda k1: kernel_hardsoft(k1, k3, omega, U, c, delta, N_bl, z_kernel)
    if U.f(0) != 0:
        Kff = ff_kernel_hardsoft_slip(omega, U.f(0), c.f(0), k3, 0)
    else:
        Kff = ff_kernel_hardsoft_noslip(omega, U.df(0), c.f(0), k3, 0)
    K = fact_mult_cont_scaled_finite(Kernel, Kff, contour, Struct(N=prec_N), branch_cut)
    return Kernel, K


def _dec(k, k3, omega, U, c, delta, N_bl, HP=0):
    z = straight_path(delta, 0.0) if HP == 0 else contour_semicircle(delta, HP)
    return comp_dec_aug18a_skip(k, k3, omega, U, c, delta, N_bl, z)


def _upstream_hard(x1, x2, y2, kappa1, k3, omega, U, c, delta, N_bl, N_orig, Q0, incident):
    """Incident field and wall pressure over a hard wall (see
    :func:`.common.incident_sheet`; ``N_orig`` is the step count the
    original used, ``10 N_bl`` in the full-field scripts)."""
    L = wall_operators(omega, U).hard
    N = N_bl if incident == "regularised" else N_orig
    I, wall = incident_sheet(x1, x2, y2, kappa1, k3, omega, U, c, delta, N, L, incident, Q0,
                             hard=True)
    return I, wall.p


# ---------------------------------------------------------------------------
# Hard -> soft
# ---------------------------------------------------------------------------
def hard_soft_continuous(U=None, c=None, delta=1.0, omega=20.0, k3=0.0, N_bl=1000, y2=0.1,
                         x1=None, x2=None, N_C=1000, angle=3 * np.pi / 8, curvature=0.01,
                         apex=None, shift1=0.5, shift2=0.5, Q0=1e-10, prec_N=2000,
                         branch_cut=np.pi / 2, incident="regularised"):
    """``HardSoft_Continuous_Nov2018.m`` -- hard wall (``x1 < 0``) to pressure
    release (``x1 > 0``), slipping or not.

    Returns a Struct with ``I_tot`` (scattered), ``I_upstream`` (incident),
    ``I_full`` (sum) and the factorised kernel ``K``.  The commented-out
    residue block of the original (for dispersion zeros) is not ported.
    """
    U = U or u_parabolic(0.3, 1, delta)
    c = c or c_const(5)
    x1 = np.arange(-5, 5 + 1e-9, 0.05) if x1 is None else np.asarray(x1)
    x2 = np.arange(0, 5 + 1e-9, 0.05) if x2 is None else np.asarray(x2)
    apex = 0.1 * omega if apex is None else apex
    kappa1 = omega / U.f(y2)

    Cp, Cm = sd_contours(N_C, omega, U.f(delta), c.f(delta), k3, apex, curvature, angle)
    C_CLp, C_CLm = cl_contours(N_C, omega, U.f(delta), U.f(0), shift1, shift2)
    I_upstream, P_up = _upstream_hard(x1, x2, y2, kappa1, k3, omega, U, c, delta, N_bl, 10 * N_bl, Q0,
                                   incident)

    Yp = y_continuous_scattering(x2, Cp.k, k3, omega, U, c, delta, N_bl)
    Ym = y_continuous_scattering(x2, Cm.k, k3, omega, U, c, delta, N_bl)
    Y_CLp = y_continuous_scattering(x2, C_CLp.k, k3, omega, U, c, delta, N_bl, 1, True)
    Y_CLm = y_continuous_scattering(x2, C_CLm.k, k3, omega, U, c, delta, N_bl, -1, True)

    Dh_p = _dec(Cp.k, k3, omega, U, c, delta, N_bl).v
    Dp_m = _dec(Cm.k, k3, omega, U, c, delta, N_bl).p
    Dp_CLp = _dec(C_CLp.k, k3, omega, U, c, delta, N_bl, 1).p
    Dp_CLm = _dec(C_CLm.k, k3, omega, U, c, delta, N_bl, -1).p

    _, K = _hardsoft_kernel(k3, omega, U, c, delta, N_bl, _contour(), prec_N, branch_cut)
    A = -1j * P_up / K.m(kappa1)
    pole = lambda k: 1 / (k - kappa1)
    Kp = A * pole(Cp.k) / (K.p(Cp.k) * Dh_p)
    Km = A * pole(Cm.k) * K.m(Cm.k) / Dp_m  # "For analyticity reasons, this is the best one."
    K_CLp = A * pole(C_CLp.k) * K.m(C_CLp.k) / Dp_CLp
    K_CLm = A * pole(C_CLm.k) * K.m(C_CLm.k) / Dp_CLm

    res = acoustic_and_cl(x1, Cp, Cm, C_CLp, C_CLm, Kp, Km, K_CLp, K_CLm, Yp, Ym, Y_CLp, Y_CLm)
    res.I_upstream = I_upstream
    res.I_full = add_fields(res.I_tot, I_upstream)
    res.update(x1=x1, x2=x2, kappa1=kappa1, K=K, A=A)
    return res


def far_field_hard_soft(U=None, c=None, delta=1.0, omega=1.0, k3=0.0, N_bl=2000, y2=0.5,
                        N_theta=1000, prec_N=5000, Q0=1e-10, plot=False,
                        incident="regularised"):
    """``FarField_HardSoft_Continuous_Nov2018.m`` -- far-field directivity of the
    hard -> soft problem (default: no-slip parabolic profile)."""
    U = U or u_parabolic(0, 1, delta)
    c = c or c_const(5)
    kappa1 = omega / U.f(y2)
    theta = np.linspace(0, np.pi, N_theta)
    sd = far_field_saddle(theta, omega, U.f(delta), c.f(delta), k3)
    S = sd.saddle
    shift = 0.5 * (0.5 * (S[0] - S[-1]))
    _, K = _hardsoft_kernel(k3, omega, U, c, delta, N_bl, _contour(shift, 3), prec_N, np.pi)

    _, P_up = _upstream_hard(None, None, y2, kappa1, k3, omega, U, c, delta, N_bl, N_bl, Q0,
                             incident)
    A = -1j * P_up / K.m(kappa1)
    dec = _dec(S, k3, omega, U, c, delta, N_bl)
    pole = 1 / (S - kappa1)
    Kp = A * pole / (K.p(S) * dec.v)
    Km = A * pole * K.m(S) / dec.p  # agree between the contours
    Y = y_continuous_scattering(delta, S, k3, omega, U, c, delta, 10)
    Y = Struct(**{f: Y[f][0] for f in FIELDS})
    FF = blend_far_field(far_field(Kp, sd.Theta, Y), far_field(Km, sd.Theta, Y), theta < np.pi / 2)
    if plot:
        from .point_source import _polar
        _polar(theta, FF.p)
    return Struct(theta=theta, Theta=sd.Theta, saddle=S, FF=FF, K=K)


# ---------------------------------------------------------------------------
# Hard -> soft -> hard (two junctions)
# ---------------------------------------------------------------------------
def _hsh_modes(k3, omega, U, c, delta, N_bl, require_unstable=True):
    P = wall_operators(omega, U).pressure
    D_l = lambda k1: dispersion_generic_nov2018(k1, k3, omega, U, c, P.l1, P.l2, delta, N_bl, 1)
    smol = 1e-12
    dD_l = lambda k1: (D_l(k1 + smol) - D_l(k1 - smol)) / (2 * smol)
    z0 = newton_raphson(D_l, dD_l, omega / U.f(delta) + 0.5j)
    zeds = np.array([z0, np.conj(z0)])
    kappap = zeds[0]
    if not np.imag(kappap) > 0:
        if require_unstable:
            raise RuntimeError("Unstable pressure-release mode not found")
        log.info("No pressure-release mode found, ignoring it")
        kappap = 0.0
    return zeds, kappap, dD_l


def hard_soft_hard(U=None, c=None, delta=1.0, omega=1.0, k3=0.0, N_bl=1000, y2=0.4,
                   x1=None, x2=None, x1_u=-5.0, x1_d=0.0, N_C=1000, angle=np.pi / 4,
                   curvature=0.01, apex=None, N_CL=1000, shift1=None, shift2=None,
                   Q0=1e-10, prec_N=2000, incident="regularised"):
    """``HardSoftHard_Continuous_Feb2019.m`` -- hard wall, pressure-release
    segment ``[x1_u, x1_d)``, hard wall again.

    The soft segment supports a conjugate pair of pressure-release modes
    (``kappap`` unstable); the analysis *assumes* this pair exists, so the
    stable case raises.  Residues give the modal field on the segment.

    .. warning:: **Bug fixed.**  The MATLAB summed
       ``I_s1 + I_s2_1 + I_s2_1`` (and the same for the critical-layer part),
       dropping the downstream forcing term ``I_s2_2`` and double-counting
       ``I_s2_1``.  Here it is ``I_s1 + I_s2_1 + I_s2_2``.
    """
    U = U or u_parabolic(0.3, 1, delta)
    c = c or c_const(5)
    x1 = np.arange(-15, 15 + 1e-9, 0.1) if x1 is None else np.asarray(x1)
    x2 = np.arange(0, 4 + 1e-9, 0.02) if x2 is None else np.asarray(x2)
    apex = 0.1 * omega if apex is None else apex
    shift1 = 0.1 * omega if shift1 is None else shift1
    shift2 = 0.1 * omega if shift2 is None else shift2
    kappa1 = omega / U.f(y2)
    l_junction = x1_d - x1_u

    Cp, Cm = sd_contours(N_C, omega, U.f(delta), c.f(delta), k3, apex, curvature, angle)
    C_CLp, C_CLm = cl_contours(N_CL, omega, U.f(delta), U.f(0), shift1, shift2)
    log.info("Contours set up")

    zeds, kappap, dD_l = _hsh_modes(k3, omega, U, c, delta, N_bl)
    dDp_kappap, dDp_zeds = dD_l(kappap), dD_l(zeds)
    log.info("Unstable modes determined, computing upstream behaviour")
    I_upstream, P_up = _upstream_hard(x1, x2, y2, kappa1, k3, omega, U, c, delta, N_bl, 10 * N_bl, Q0,
                                   incident)

    log.info("Computing wall-normal behaviour")
    Yp = y_continuous_scattering(x2, Cp.k, k3, omega, U, c, delta, N_bl)
    Ym = y_continuous_scattering(x2, Cm.k, k3, omega, U, c, delta, N_bl)
    Y_CLp = y_continuous_scattering(x2, C_CLp.k, k3, omega, U, c, delta, N_bl, 1, True)
    Y_CLm = y_continuous_scattering(x2, C_CLm.k, k3, omega, U, c, delta, N_bl, -1, True)

    decp, decm = _dec(Cp.k, k3, omega, U, c, delta, N_bl), _dec(Cm.k, k3, omega, U, c, delta, N_bl)
    dec_CLp = _dec(C_CLp.k, k3, omega, U, c, delta, N_bl, 1)
    dec_CLm = _dec(C_CLm.k, k3, omega, U, c, delta, N_bl, -1)
    Dh_kappap = _dec(kappap, k3, omega, U, c, delta, N_bl, 1).v[0]

    Kernel, K = _hardsoft_kernel(k3, omega, U, c, delta, N_bl, _contour(), prec_N, np.pi)
    Kpp, Kmm = K.p(Cp.k), K.m(Cm.k)
    Km_CLp, Km_CLm = K.m(C_CLp.k), K.m(C_CLm.k)
    Km_kappa1 = K.m(kappa1)
    Kp_kappa1 = Kernel(kappa1) / Km_kappa1
    Km_kappap = K.m(kappap)
    Km_zeds = K.m(zeds)
    log.info("Kernel found and factorised, setting up integrals")

    A_s1 = -P_up / (1j * Km_kappa1)
    A_s2_1 = (P_up * Dh_kappap * Km_kappap ** 2 * np.exp(-1j * kappap * l_junction)
              / (1j * Km_kappa1 * (kappap - kappa1) * dDp_kappap))
    A_s2_2 = P_up * np.exp(-1j * kappa1 * l_junction) / (1j * Kp_kappa1)
    ca = lambda k, a: 1 / (k - a)

    # Upstream junction
    s1 = acoustic_and_cl(
        x1, Cp, Cm, C_CLp, C_CLm,
        A_s1 * ca(Cp.k, kappa1) / (Kpp * decp.v), A_s1 * ca(Cm.k, kappa1) * Kmm / decm.p,
        A_s1 * ca(C_CLp.k, kappa1) * Km_CLp / dec_CLp.p, A_s1 * ca(C_CLm.k, kappa1) * Km_CLm / dec_CLm.p,
        Yp, Ym, Y_CLp, Y_CLm, x1_shift=x1_u)

    # Downstream junction: modal "unscattering" (a = kappap) and forcing (a = kappa1)
    def downstream(A, a):
        return acoustic_and_cl(
            x1, Cp, Cm, C_CLp, C_CLm,
            A * ca(Cp.k, a) * Kpp / decp.p, A * ca(Cm.k, a) / (Kmm * decm.v),
            A * ca(C_CLp.k, a) / (Km_CLp * dec_CLp.v), A * ca(C_CLm.k, a) / (Km_CLm * dec_CLm.v),
            Yp, Ym, Y_CLp, Y_CLm, x1_shift=x1_d)

    s2_1, s2_2 = downstream(A_s2_1, kappap), downstream(A_s2_2, kappa1)

    log.info("Computing modal contributions")
    x1 = row(x1)
    mid = (x1 >= x1_u) & (x1 < x1_d)
    X_mid = fourier_nov2018(x1[mid] - x1_u, zeds)
    Y_mid = y_continuous_scattering(x2, zeds, k3, omega, U, c, delta, N_bl)
    Res_weight = -2j * np.pi * A_s1 * Km_zeds / ((zeds - kappa1) * dDp_zeds)
    I_res = Struct()
    for f in FIELDS:
        full = np.zeros((x1.size, row(x2).size), dtype=complex)
        full[mid] = sum(Res_weight[j] * np.outer(X_mid[:, j], Y_mid[f][:, j]) for j in range(zeds.size))
        I_res[f] = full

    I_acoustic = add_fields(s1.I, s2_1.I, s2_2.I)
    I_CL = add_fields(s1.I_CL, s2_1.I_CL, s2_2.I_CL)
    I_scattered = add_fields(I_res, I_CL, I_acoustic)
    I_scattered1 = add_fields(I_res, s1.I_CL, s1.I)
    log.info("All done - well done")
    return Struct(x1=x1, x2=x2, kappa1=kappa1, kappap=kappap, zeds=zeds, K=K,
                  I_upstream=I_upstream, s1=s1, s2_1=s2_1, s2_2=s2_2, I_res=I_res,
                  I_acoustic=I_acoustic, I_CL=I_CL, I_scattered=I_scattered,
                  I_full=add_fields(I_upstream, I_scattered),
                  I_scattered1=I_scattered1, I_full1=add_fields(I_upstream, I_scattered1))


def far_field_hard_soft_hard(U=None, c=None, delta=1.0, omega=4.0, k3=0.0, N_bl=1000, y2=0.4,
                             x1_u=-7.5, x1_d=0.0, N_theta=1000, prec_N=5000, Q0=1e-10,
                             incident="regularised"):
    """``FarField_HardSoftHard_Continuous_Feb2019.m`` -- far field of the
    double junction; also returns each contribution (``FF_s1``, ``FF_s2_1``,
    ``FF_s2_2``).  Without an unstable mode ``kappap`` is set to 0 and its
    term dropped (as in the original)."""
    U = U or u_parabolic(0.3, 1, delta)
    c = c or c_const(5)
    kappa1 = omega / U.f(y2)
    l_junction = x1_d - x1_u
    theta = np.linspace(0, np.pi, N_theta)
    sd = far_field_saddle(theta, omega, U.f(delta), c.f(delta), k3)
    S = sd.saddle

    zeds, kappap, dD_l = _hsh_modes(k3, omega, U, c, delta, N_bl, require_unstable=False)
    dDp_kappap = dD_l(kappap)
    _, P_up = _upstream_hard(None, None, y2, kappa1, k3, omega, U, c, delta, N_bl, N_bl, Q0,
                             incident)
    Y = y_continuous_scattering(delta, S, k3, omega, U, c, delta, 10)
    Y = Struct(**{f: Y[f][0] for f in FIELDS})

    dec0 = _dec(S, k3, omega, U, c, delta, N_bl)
    Dh_kappap = _dec(kappap, k3, omega, U, c, delta, N_bl).v[0]
    shift = 0.5 * (0.5 * (S[0] - S[-1]))
    Kernel, K = _hardsoft_kernel(k3, omega, U, c, delta, N_bl, _contour(shift, 3), prec_N, np.pi)
    Kp0, Km0 = K.p(S), K.m(S)
    Km_kappa1 = K.m(kappa1)
    Kp_kappa1 = Kernel(kappa1) / Km_kappa1
    Km_kappap = K.m(kappap)

    A_s1 = -P_up / (1j * Km_kappa1)
    A_s2_1 = 0.0
    if kappap != 0:
        A_s2_1 = (P_up * Dh_kappap * Km_kappap ** 2 * np.exp(-1j * kappap * l_junction)
                  / (1j * Km_kappa1 * (kappap - kappa1) * dDp_kappap))
    A_s2_2 = P_up * np.exp(-1j * kappa1 * l_junction) / (1j * Kp_kappa1)
    C1, Cp_ = 1 / (S - kappa1), 1 / (S - kappap)

    parts = {
        "s1": (A_s1 * C1 / (Kp0 * dec0.v), A_s1 * C1 * Km0 / dec0.p),
        "s2_1": (A_s2_1 * Cp_ * Kp0 / dec0.p, A_s2_1 * Cp_ / (Km0 * dec0.v)),
        "s2_2": (A_s2_2 * C1 * Kp0 / dec0.p, A_s2_2 * C1 / (Km0 * dec0.v)),
    }
    down = theta < np.pi / 2
    out = Struct(theta=theta, Theta=sd.Theta, saddle=S, K=K, kappap=kappap)
    for name, (kp, km) in parts.items():
        out["FF_" + name] = blend_far_field(far_field(kp, sd.Theta, Y), far_field(km, sd.Theta, Y), down)
    Kp = sum(p[0] for p in parts.values())
    Km = sum(p[1] for p in parts.values())
    out.FF = blend_far_field(far_field(Kp, sd.Theta, Y), far_field(Km, sd.Theta, Y), down)
    return out


# ---------------------------------------------------------------------------
# Impedance -> impedance
# ---------------------------------------------------------------------------
def _impedance_kernel(Zu, Zd, k3, omega, U, c, delta, N_bl, contour, prec_N):
    Lu, Ld = impedance_up_down(Zu, Zd, omega, U.f(0), U.df(0))
    z_kernel = contour_semicircle(delta, 1)
    Kernel = lambda k1: kernel_generic(k1, k3, omega, U, c, delta, Lu, Ld, N_bl, z_kernel)
    Kff = ff_kernel_impedance(Zu, Zd, omega, U.f(0), c.f(0), k3, 0)
    return Lu, Ld, fact_mult_cont_scaled_finite(Kernel, Kff, contour, Struct(N=prec_N), np.pi)


def imp_imp_continuous(Zu=1 - 1j, Zd=2 * (1 - 1j), U=None, c=None, delta=1.0, omega=1.0,
                       k3=0.0, N_bl=1000, y2=0.5, x1=None, x2=None, N_C=1000,
                       angle=np.pi / 4, curvature=0.01, apex=None, N_CL=1000,
                       shift1=None, shift2=None, Q0=1e-10, prec_N=2000,
                       incident="regularised"):
    """``ImpImp_Continuous_Jan2019.m`` -- junction between two impedance
    linings ``Zu`` (``x1 < 0``) and ``Zd`` (either may be ``inf``)."""
    U = U or u_parabolic(0, 1, delta)
    c = c or c_const(5)
    x1 = np.arange(-5, 5 + 1e-9, 0.02) if x1 is None else np.asarray(x1)
    x2 = np.arange(0, 5 + 1e-9, 0.02) if x2 is None else np.asarray(x2)
    apex = 0.1 * omega if apex is None else apex
    shift1 = 0.1 * omega if shift1 is None else shift1
    shift2 = 0.1 * omega if shift2 is None else shift2
    kappa1 = omega / U.f(y2)

    Cp, Cm = sd_contours(N_C, omega, U.f(delta), c.f(delta), k3, apex, curvature, angle)
    C_CLp, C_CLm = cl_contours(N_CL, omega, U.f(delta), U.f(0), shift1, shift2)
    Lu, Ld, K = _impedance_kernel(Zu, Zd, k3, omega, U, c, delta, N_bl, _contour(), prec_N)

    N_inc = N_bl if incident == "regularised" else 10 * N_bl
    I_upstream, wall = incident_sheet(x1, x2, y2, kappa1, k3, omega, U, c, delta, N_inc, Lu,
                                      incident, Q0)
    forcing = Ld.l1(kappa1) * wall.dphi + Ld.l2(kappa1) * wall.phi

    Yp = y_continuous_scattering(x2, Cp.k, k3, omega, U, c, delta, N_bl)
    Ym = y_continuous_scattering(x2, Cm.k, k3, omega, U, c, delta, N_bl)
    Y_CLp = y_continuous_scattering(x2, C_CLp.k, k3, omega, U, c, delta, N_bl, 1, True)
    Y_CLm = y_continuous_scattering(x2, C_CLm.k, k3, omega, U, c, delta, N_bl, -1, True)

    def D(L, k, HP):
        dec = _dec(k, k3, omega, U, c, delta, N_bl, HP)
        return L.l1(k) * dec.dphi + L.l2(k) * dec.phi

    Du_p, Dd_m = D(Lu, Cp.k, 0), D(Ld, Cm.k, 0)
    Dd_CLp, Dd_CLm = D(Ld, C_CLp.k, 1), D(Ld, C_CLm.k, -1)

    A = -1j * forcing / K.m(kappa1)
    pole = lambda k: 1 / (k - kappa1)
    res = acoustic_and_cl(
        x1, Cp, Cm, C_CLp, C_CLm,
        A * pole(Cp.k) / (K.p(Cp.k) * Du_p), A * pole(Cm.k) * K.m(Cm.k) / Dd_m,
        A * pole(C_CLp.k) * K.m(C_CLp.k) / Dd_CLp, A * pole(C_CLm.k) * K.m(C_CLm.k) / Dd_CLm,
        Yp, Ym, Y_CLp, Y_CLm)
    res.I_upstream = I_upstream
    res.I_full = add_fields(res.I_tot, I_upstream)
    res.update(x1=x1, x2=x2, kappa1=kappa1, K=K, A=A)
    return res


def far_field_imp_imp(R, Z=1 - 1j, U=None, c=None, delta=1.0, omega=1.0, k3=0.0, N_bl=1000,
                      y2=0.5, N_theta=1000, prec_N=5000, Q0=1e-10, incident="regularised"):
    """``FarField_ImpImp_Continuous_Feb2019.m`` -- far field of an impedance
    junction with ``Zu = Z / R`` and ``Zd = R Z``.

    .. warning:: **Missing input in the original.**  ``R`` was never defined
       in the script (it must have been set in the workspace beforehand); it
       is a required argument here.
    """
    U = U or u_parabolic(0, 1, delta)
    c = c or c_const(5)
    Zu, Zd = Z / R, R * Z
    kappa1 = omega / U.f(y2)
    theta = np.linspace(0, np.pi, N_theta)
    sd = far_field_saddle(theta, omega, U.f(delta), c.f(delta), k3)
    S = sd.saddle
    Y = y_continuous_scattering(delta, S, k3, omega, U, c, delta, 10)
    Y = Struct(**{f: Y[f][0] for f in FIELDS})
    shift = 0.5 * (0.5 * (S[0] - S[-1]))
    Lu, Ld, K = _impedance_kernel(Zu, Zd, k3, omega, U, c, delta, N_bl, _contour(shift, 3), prec_N)

    _, wall = incident_sheet(None, None, y2, kappa1, k3, omega, U, c, delta, N_bl, Lu, incident, Q0)
    forcing = Ld.l1(kappa1) * wall.dphi + Ld.l2(kappa1) * wall.phi
    dec = _dec(S, k3, omega, U, c, delta, N_bl)
    Du = Lu.l1(S) * dec.dphi + Lu.l2(S) * dec.phi
    Dd = Ld.l1(S) * dec.dphi + Ld.l2(S) * dec.phi
    A = -forcing / (1j * K.m(kappa1))
    C1 = 1 / (S - kappa1)
    FF = blend_far_field(far_field(A * C1 / (K.p(S) * Du), sd.Theta, Y),
                         far_field(A * C1 * K.m(S) / Dd, sd.Theta, Y), theta < np.pi / 2)
    return Struct(theta=theta, Theta=sd.Theta, saddle=S, FF=FF, K=K)
