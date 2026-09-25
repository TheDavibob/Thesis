"""Steepest-descent contours and far-field (saddle-point) evaluation.

Ports of ``FixedStep/PridmoreBrown/Steepest_Descent_Contour*.m`` and the two
exploratory ``Sep_18_Inversion`` files.

Far field: in Mach-stretched polar angle ``Theta = atan2(beta x2, x1)``
(``beta = sqrt(1 - M^2)``) the saddle point is

    k1* = K0 cos(Theta) - M k0 / (1 - M^2),
    K0  = k0 / (1 - M^2) * sqrt(1 - beta^2 k3^2 / k0^2),

and the far-field directivity is ``K(k1*) * sin(Theta) * Y(k1*)`` (up to the
usual ``sqrt(pi / (R K0))`` factor, omitted throughout the original).
"""
from __future__ import annotations

import numpy as np

from ._utils import Path, Struct, log
from .complex_funcs import sqrt_a
from .pridmore_brown import (comp_dec_aug18a_skip, comp_hard_aug18a_skip,
                             dec_aug18a_gridx_skip, gamma_ff, hard_aug18a_gridx_skip)
from .quadrature import lgwt


def steepest_descent_contour(x1, x2, omega, Uinf, cinf, k3, N) -> Struct:
    """``Steepest_Descent_Contour.m`` -- SD contour for every observer point.

    ``x1``, ``x2`` are either vectors (meshgridded here) or meshgrid arrays.
    Returns ``k1`` with ``SD`` (contour points, shape ``(N2, N1, N)``),
    quadrature ``t``, ``w``, ``dt`` and ``wdt`` (weights x Jacobian),
    ``X1``, ``X2``, ``Theta`` and the physical constants ``k0, K0, beta, M``.
    Not yet valid for negative ``x2`` (original comment).
    """
    x1 = np.asarray(x1, dtype=float)
    x2 = np.asarray(x2, dtype=float)
    if x1.ndim <= 1:
        X1, X2 = np.meshgrid(np.atleast_1d(x1), np.atleast_1d(x2))
    else:
        X1, X2 = x1, x2
    M = Uinf / cinf
    k0 = omega / cinf
    beta = np.sqrt(1 - M ** 2)
    Theta = np.arctan2(beta * X2, X1)
    K0 = k0 / (1 - M ** 2) * np.sqrt(complex(1 - beta ** 2 * k3 ** 2 / k0 ** 2))
    K0 = np.real_if_close(K0)
    t, w = lgwt(N, -1, 1)
    t = t[::-1]
    T = t[None, None, :]
    Th = Theta[..., None]
    rt = np.sqrt(1 - T ** 2)
    K1_SD = K0 * (np.cos(Th) + T / rt * np.sin(Th)
                  - 1j * (T ** 2 / rt * np.cos(Th) - T * np.sin(Th)))
    dK1 = K0 * (np.sin(Th) / rt ** 3
                + 1j * (T * (T ** 2 - 2) / rt ** 3 * np.cos(Th) + np.sin(Th)))
    W = np.broadcast_to(w[None, None, :], K1_SD.shape)
    return Struct(SD=K1_SD - M * k0 / (1 - M ** 2), t=t, dt=dK1, wdt=W * dK1, w=w,
                  X1=X1, X2=X2, Theta=Theta, k0=k0, K0=K0, beta=beta, M=M)


def far_field_saddle(theta, omega, Uinf, cinf, k3=0.0, lower=False) -> Struct:
    """Saddle points for observers at polar angles ``theta`` on the unit circle.

    Replaces the idiom used by every ``FarField_*`` script::

        [k1] = Steepest_Descent_Contour(cos(theta), sin(theta), ..., 1);
        Theta = diag(k1.Theta).';
        Saddle = k1.K0.*cos(Theta) - k1.M.*k1.k0/(1-k1.M^2);

    which built an ``N_theta x N_theta`` meshgrid only to take its diagonal.
    ``lower=True`` reflects ``x2 -> -x2`` (used for the ``x2 < 0`` side).
    """
    theta = np.asarray(theta, dtype=float)
    x1, x2 = np.cos(theta), np.sin(theta)
    if lower:
        x2 = -x2
    M = Uinf / cinf
    k0 = omega / cinf
    beta = np.sqrt(1 - M ** 2)
    Theta = np.arctan2(beta * x2, x1)
    K0 = np.real_if_close(k0 / (1 - M ** 2) * np.sqrt(complex(1 - beta ** 2 * k3 ** 2 / k0 ** 2)))
    return Struct(Theta=Theta, saddle=K0 * np.cos(Theta) - M * k0 / (1 - M ** 2),
                  K0=K0, k0=k0, M=M, beta=beta)


def steepest_descent_contour_anon(omega, Uinf, cinf, k3):
    """``Steepest_Descent_Contour_Anon.m`` -- the same maths as callables.

    Marked "Not finished, essentially a check the mathematics is correct" in
    the original.  Returns ``(f, F, K1_SD, k1_SD)``: the phase in physical and
    shifted wavenumber, and the SD contour ``(t, Theta) -> k``.
    """
    M = Uinf / cinf
    k0 = omega / cinf
    beta = np.sqrt(1 - M ** 2)
    K0 = k0 / (1 - M ** 2) * np.sqrt(complex(1 - beta ** 2 * k3 ** 2 / k0 ** 2))
    if np.imag(K0) != 0:
        log.warning("Acoustic wavenumber is complex")
    gamma, _, _ = gamma_ff(k3, omega, cinf, Uinf)
    f = lambda k1, Th: 1j * k1 * np.cos(Th) + gamma(k1) * np.sin(Th) / beta
    F = lambda K1, Th: (1j * K1 * np.cos(Th)
                        + sqrt_a(K1 - K0, 3 * np.pi / 2) * sqrt_a(K1 + K0, np.pi / 2) * np.sin(Th))
    K1_SD = lambda t, Th: K0 * (np.cos(Th) + t / np.sqrt(1 - t ** 2) * np.sin(Th)
                                - 1j * (t ** 2 / np.sqrt(1 - t ** 2) * np.cos(Th) - t * np.sin(Th)))
    k1_SD = lambda t, Th: K1_SD(t, Th) - M * k0 / (1 - M ** 2)
    return f, F, K1_SD, k1_SD


def oct_18_ff_sd(theta, y2, k3, omega, U, c, delta, N_BL=1000) -> Struct:
    """``Oct_18_FF_SD.m`` -- far field of a vortex-sheet source at height
    ``y2`` above a hard wall (first-order saddle contribution only).

    .. warning:: **Bug fixed.**  The wall-to-edge path was
       ``z.f = delta*(1-t)`` with ``z.df = -delta*t``; the derivative is
       ``-delta``.
    """
    sd = far_field_saddle(theta, omega, U.f(delta), c.f(delta), k3)
    K1 = sd.saddle
    zh = Path(f=lambda t: t * y2, df=lambda t: y2 * np.ones(np.shape(t)))
    hard_y2 = comp_hard_aug18a_skip(K1, k3, omega, U, c, N_BL, zh)
    log.info("Vortex sheet contribution computed")
    z = Path(f=lambda t: delta * (1 - t), df=lambda t: -delta * np.ones(np.shape(t)))
    dec_0 = comp_dec_aug18a_skip(K1, k3, omega, U, c, delta, N_BL, z)
    Disp_hard = dec_0.v
    Wronskian = Disp_hard / (c.f(0) ** 2 * (1j * (omega - U.f(0) * K1)) ** 5)
    log.info("Hard-wall contribution computed")
    Cauchy = 1 / (2j * np.pi * (K1 - omega / U.f(y2)))
    Conv_rat = (omega - U.f(y2) * K1) ** 4 / (omega - U.f(0) * K1) ** 4
    out = Struct(Theta=sd.Theta, Disp_hard=Disp_hard, Cauchy=Cauchy,
                 hard_y2=hard_y2, Wronskian=Wronskian)
    for name in ("phi", "p", "v", "dphi"):
        out[name] = Cauchy * hard_y2[name] / Wronskian * Conv_rat
    return out


def sep18_integrate_sd(x1, x2, y2, k3, omega, U, c, delta, N_SD=100, N_BL=100):
    """``Sep18_Integrate_SD.m`` -- **UNFINISHED in the original**
    ("A HECK OF A LOT OF AWKWARD WORK IS NEEDED").

    Direct numerical steepest-descent integration of the point-source field
    above a hard wall.  Ported line-for-line, including its known gaps:

    * the ``Dec``/``Hard`` values at ``y2`` deliberately store ``phi`` in
      every field (the original comment explains this);
    * the Wronskian normalisation is flagged by the author as "not yet
      correct";
    * the hard-wall solution above the boundary layer is "FUTURE WORK".

    Returns ``(I, Integrand, Hard, Dec, Wron, k1)``.
    """
    log.info("Setting up integration")
    k1 = steepest_descent_contour(x1, x2, omega, U.f(delta), c.f(delta), k3, N_SD)
    X1b = np.broadcast_to(k1.X1[..., None], k1.SD.shape)
    X2b = np.broadcast_to(k1.X2[..., None], k1.SD.shape)
    K1 = k1.SD
    fields = ("phi", "dphi", "p", "v")
    Dec = Struct(**{f: np.zeros(K1.shape, dtype=complex) for f in fields})
    Hard = Struct(**{f: np.zeros(K1.shape, dtype=complex) for f in fields})

    log.info("Computing auxiliary functions within the boundary-layer")
    m_dec = (X2b > y2) & (X2b < delta)
    m_hard = (X2b < y2) & (X2b < delta)
    dec = dec_aug18a_gridx_skip(X2b[m_dec], K1[m_dec], k3, omega, U, c, delta, N_BL)
    hard = hard_aug18a_gridx_skip(X2b[m_hard], K1[m_hard], k3, omega, U, c, delta, N_BL)
    for f in fields:
        Dec[f][m_dec] = dec[f]
        Hard[f][m_hard] = hard[f]

    log.info("Computing auxiliary functions outside the boundary-layer")
    gammainf, _, _ = gamma_ff(k3, omega, c.f(delta), U.f(delta))
    out = X2b >= delta
    G = gammainf(K1[out])
    C_inf = 1j * (omega - U.f(delta) * K1[out])
    dC_inf = -1j * U.df(delta) * K1[out]
    Dec.phi[out] = np.exp(-G * X2b[out])
    Dec.dphi[out] = -G * np.exp(-G * X2b[out])
    Dec.p[out] = -C_inf ** 3 * Dec.phi[out]
    Dec.v[out] = -3 * dC_inf * C_inf * Dec.phi[out] - C_inf ** 2 * Dec.dphi[out]

    log.info("Computing auxiliary functions at y2")
    m_dy = (X2b <= y2) & (X2b < delta)
    m_hy = X2b >= y2
    dec_y2 = dec_aug18a_gridx_skip(y2 * np.ones(m_dy.sum()), K1[m_dy], k3, omega, U, c, delta, N_BL)
    hard_y2 = hard_aug18a_gridx_skip(y2 * np.ones(m_hy.sum()), K1[m_hy], k3, omega, U, c, delta, N_BL)
    for f in fields:  # sic: phi stored in every field
        Dec[f][m_dy] = dec_y2.phi
        Hard[f][m_hy] = hard_y2.phi

    log.info("Computing dispersion relationship")
    dec_0 = dec_aug18a_gridx_skip(np.zeros(K1.size), K1.ravel(), k3, omega, U, c, delta, N_BL)
    Disp_hard = dec_0.v.reshape(K1.shape)
    C_0 = 1j * (omega - U.f(0) * K1)
    C_y2 = 1j * (omega - U.f(y2) * K1)
    Wron = (c.f(0) ** 2 * C_0 ** 4) / (c.f(y2) ** 2 * C_y2 ** 4) * Disp_hard

    Integrand = Struct(Fourier=np.exp(-1j * K1 * X1b), Weight=k1.wdt,
                       Cauchy=1 / (2j * np.pi * (K1 - omega / U.f(y2))), Disp_rec=1 / Wron)
    common = Integrand.Fourier * Integrand.Weight * Integrand.Cauchy * Integrand.Disp_rec
    I = Struct()
    for f in fields:
        Integrand[f] = common * Dec[f] * Hard[f]
        I[f] = Integrand[f].sum(axis=2)
    return I, Integrand, Hard, Dec, Wron, k1
