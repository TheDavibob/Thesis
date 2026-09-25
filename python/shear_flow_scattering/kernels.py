"""Wiener-Hopf kernels and their analytic far-field factorisations.

Ports of ``Factorisation/Generic_BC``, ``Generic_Wake`` and ``NoSlip_HardSoft``.

A kernel ``K(k1)`` is a ratio of dispersion relations evaluated with the
boundary-layer solver.  The numerical factorisation (:mod:`.factorisation`)
needs ``K -> 1`` at infinity, so each kernel comes with ``FF*`` -- a
far-field approximation ``Kff`` whose ``+``/``-`` factors are known in closed
form (``Kff.f = Kff.p * Kff.m``).  ``W`` is an integer "winding number" that
may be added to the exponents to keep ``log(K/Kff)`` from wrapping round 0.
"""
from __future__ import annotations

import numpy as np

from ._utils import Path, Struct, log
from .complex_funcs import power_a
from .pridmore_brown import comp_dec_aug18a_skip, gamma_ff


def _ff_power(A, alpha, beta, kbp, kbm) -> Struct:
    Kff = Struct(p=lambda k1: A * power_a(k1 - kbp, alpha, 3 * np.pi / 2),
                 m=lambda k1: power_a(k1 - kbm, beta, np.pi / 2))
    Kff.f = lambda k1: Kff.m(k1) * Kff.p(k1)
    return Kff


def _kb(omega, c0, k3):
    kb = np.sqrt(complex((omega / c0) ** 2 - k3 ** 2))
    return kb, -kb


# ---------------------------------------------------------------------------
# Generic linear boundary conditions / impedance
# ---------------------------------------------------------------------------
def impedance_up_down(Zu, Zd, omega, U0, sigma0):
    """``ImpedanceUpDown.m`` -- wall operators ``L phi = l1 phi' + l2 phi``
    for an impedance ``Zu`` upstream and ``Zd`` downstream (``inf`` = hard).

    ``U0 = U(0)`` and ``sigma0 = U'(0)``.  Returns ``(Lu, Ld)``.
    """
    C0 = lambda k1: 1j * (omega - U0 * k1)
    dC0 = lambda k1: -1j * sigma0 * k1

    def op(Z):
        if np.isinf(Z):
            return Struct(l1=lambda k1: -C0(k1) ** 2, l2=lambda k1: -3 * C0(k1) * dC0(k1))
        return Struct(l1=lambda k1: 1j * omega * Z * (-C0(k1) ** 2),
                      l2=lambda k1: 1j * omega * Z * (-3 * C0(k1) * dC0(k1)) + C0(k1) ** 4)

    return op(Zu), op(Zd)


def kernel_generic(k1, k3, omega, U, c, delta, Lu, Ld, N, z: Path):
    """``Kernel_Generic.m`` -- ``K = D_d / D_u`` for a change of wall BC from
    ``Lu phi = 0`` to ``Ld phi = 0``.  ``z`` should run from ``delta`` to 0
    (an upper-half-plane semicircle "works swimmingly")."""
    dec = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N, z)
    Du = Lu.l1(k1) * dec.dphi + Lu.l2(k1) * dec.phi
    Dd = Ld.l1(k1) * dec.dphi + Ld.l2(k1) * dec.phi
    return np.reshape(Dd / Du, np.shape(k1))


def ff_kernel_impedance(Zu, Zd, omega, U0, c0, k3, W) -> Struct:
    """``FFKernel_Impedance.m`` -- far-field factors for :func:`kernel_generic`
    with two impedances.

    With flow (``U0 != 0``) the kernel tends to 1 at both ends and the
    "tricky exponent" ``L`` is zero; without flow the constant is ``Zd/Zu``.
    """
    A = Zd / Zu if U0 == 0 else 1.0
    L = 0.0
    kbp, kbm = _kb(omega, c0, k3)
    return _ff_power(A, L + W, -L - W, kbp, kbm)


def ff_kernel_impedance_historic(Zu, Zd, omega, U0, c0, k3, W) -> Struct:
    """``FFKernel_Impedance_Historic.m`` -- **known wrong, kept for reference.**

    The original is headed "THIS IS WRONG WRONG WRONG DO NOT USE" (the
    boundary conditions were mis-specified).  Use :func:`ff_kernel_impedance`.
    """
    log.warning("ff_kernel_impedance_historic is documented as incorrect in the original code")
    beta = np.sqrt(1 - U0 ** 2 / c0 ** 2)
    redU = U0 / (beta * omega)
    Kp = (Zd - redU) / (Zu - redU)
    Km = (Zd + redU) / (Zu + redU)
    L = 1 / (2j * np.pi) * np.log(complex(Km / Kp))  # complex: MATLAB's log(-x) is complex
    kbp, kbm = _kb(omega, c0, k3)
    return _ff_power(Kp, L + W, -L - W, kbp, kbm)


# ---------------------------------------------------------------------------
# Hard -> soft (pressure release) transition
# ---------------------------------------------------------------------------
def kernel_hardsoft(k1, k3, omega, U, c, delta, N, z: Path):
    """``Kernel_HardSoft.m`` -- ``K = D_p / D_h = p(0) / v(0)`` of the
    decaying solution."""
    dec = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N, z)
    return np.reshape(dec.p / dec.v, np.shape(k1))


def ff_kernel_hardsoft_noslip(omega, sigma, c0, k3, W) -> Struct:
    """``FFKernel_HardSoft_NoSlip.m`` -- far field for a no-slip wall
    (``sigma = U'(0)``), factorised "as per Singh"."""
    A = 1j * (sigma - omega)
    # complex(): for sigma > omega the argument is negative; MATLAB's log then
    # returns log|x| + i pi, whereas np.log of a negative float gives NaN.
    L = 1 / (2j * np.pi) * np.log(complex((sigma + omega) / (omega - sigma)))
    kbp, kbm = _kb(omega, c0, k3)
    return _ff_power(A, -0.5 + L + W, -0.5 - L - W, kbp, kbm)


def ff_kernel_hardsoft_slip(omega, U0, c0, k3, W) -> Struct:
    """``FFKernel_HardSoft_Slip.m`` -- far field with slip velocity ``U0 = U(0)``.

    ``Kff = -C / gamma`` (evaluated with wall values); ``W != 0`` inserts the
    factor ``((k1 - kp)/(k1 - km))^W`` to move the log's winding number.
    """
    gamma, gammap, gammam = gamma_ff(k3, omega, c0, U0)
    C = lambda k1: 1j * (omega - U0 * k1)
    M = U0 / c0
    k0 = omega / c0
    # NOTE: same k3 sign slip as the scripts (see pridmore_brown.branch_points);
    # irrelevant when W == 0, which is how every script calls this.
    kb1 = -M / (1 - M ** 2) * k0 + np.sqrt(k0 ** 2 + k3 ** 2 * (1 - M ** 2)) / (1 - M ** 2)
    kb2 = -M / (1 - M ** 2) * k0 - np.sqrt(k0 ** 2 + k3 ** 2 * (1 - M ** 2)) / (1 - M ** 2)
    kp, km = kb2, kb1
    return Struct(
        f=lambda k1: -(C(k1) / gamma(k1)) * (k1 - kp) ** W / (k1 - km) ** W,
        p=lambda k1: -(C(k1) / gammap(k1)) / (k1 - km) ** W,
        m=lambda k1: (1 / gammam(k1)) * (k1 - kp) ** W)


# ---------------------------------------------------------------------------
# Trailing edge / wake
# ---------------------------------------------------------------------------
def wake_kernel(k1, k3, omega, U1, c1, delta1, U2, c2, delta2, Z, N, z1: Path, z2: Path):
    """``Wake_Kernel.m`` -- lined plate (impedance ``Z``; ``inf`` hard, ``0``
    soft) upstream, a wake (continuity of ``p`` and ``v``) downstream.

    Flow above: ``U1, c1, delta1``; below: ``U2, c2, delta2`` as functions of
    ``-x2``.  ``z1``/``z2`` are the boundary-layer integration paths.
    """
    dec1 = comp_dec_aug18a_skip(k1, k3, omega, U1, c1, delta1, N, z1)
    dec2 = comp_dec_aug18a_skip(k1, k3, omega, U2, c2, delta2, N, z2)
    k = np.ravel(k1)
    C0 = 1j * (omega - U1.f(0) * k)
    Dw = -dec1.p * dec2.v - dec2.p * dec1.v
    if np.isinf(Z):
        K = Dw / (dec1.v * -dec2.v)
    elif Z == 0:  # a fundamentally different case
        K = -Dw / (dec1.p * dec2.p)
    else:
        DZ1 = 1j * omega * Z * dec1.v - C0 * dec1.p
        DZ2 = -1j * omega * Z * dec2.v - C0 * dec2.p
        K = 1j * omega * Z * Dw / (DZ1 * DZ2)
    return np.reshape(K, np.shape(k1))


def ff_wake_kernel(k3, omega, U0, c0, sigma1, sigma2, Z, W) -> Struct:
    """``FFWake_Kernel.m`` -- far-field factors of :func:`wake_kernel` for
    ``Z = inf``, ``Z = 0`` or finite, with or without slip.

    ``sigma1 = U1'(0)``, ``sigma2 = -U2'(0)`` (note the sign in the callers).
    """
    sigma0 = (sigma1 - sigma2) / 2
    M = U0 / c0
    kbp = np.sqrt(complex((omega / c0) ** 2 - k3 ** 2))
    kbm = -kbp
    beta0 = np.sqrt(1 - M ** 2)
    if np.isinf(Z):
        if U0 != 0:
            Kp, Km, dep = 2j * U0 / beta0, -2j * U0 / beta0, 0
        else:
            Kp, Km, dep = 2j * (sigma0 - omega), 2j * (sigma0 + omega), -1
    elif Z == 0:
        if U0 != 0:
            Kp, Km, dep = -2j * beta0 / U0, 2j * beta0 / U0, 0
        else:
            Kp = -1j * (1 / (sigma1 - omega) + 1 / (-sigma2 - omega))
            Km = 1j * (1 / (sigma1 + omega) + 1 / (-sigma2 + omega))
            dep = 1
    else:
        Zed = 1j * omega * Z
        if U0 != 0:
            Kp, Km, dep = -2j * beta0 * Zed / U0 ** 3, 2j * beta0 * Zed / U0 ** 3, -2
        else:
            Kp, Km, dep = 2j * (sigma0 - omega) / Zed, -2j * (sigma0 + omega) / Zed, -1
    L = 1 / (2j * np.pi) * np.log(Km / Kp)
    return _ff_power(Kp, L + dep + W, -L - W, kbp, kbm)
