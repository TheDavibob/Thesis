"""Distributed (continuous-in-height) vorticity sources.

Ports of the *functions* in ``GenericInversion/Scattering/Distributed_Source``:
``Fbar``, ``kappa1_Quad``, ``Integrand_Distributed``, ``Distributed_WallPressure``,
``WallPressure``, ``Upstream_Distributed_{A,B,Integrated}`` and ``Scrappy_I0``.

A vortex sheet at height ``y2`` convects at ``U(y2)`` and so forces at
``kappa1 = omega / U(y2)``; a distribution over ``y2`` therefore becomes an
integral over ``kappa1`` from ``omega/U(delta)`` to ``omega/U(0)``, i.e.
along the critical-layer branch cut.
"""
from __future__ import annotations

import numpy as np

from ._utils import FIELDS, Struct, log, pchip, row
from .contours import path_to, straight_path
from .inversion import y_continuous_scattering, y_greens_generic
from .pridmore_brown import (comp_dec_aug18a, comp_dec_aug18a_skip, comp_generic_aug18a,
                             comp_generic_aug18a_skip, gamma_ff)
from .quadrature import lgwt


def fbar(k1, kappa1, dkappa1, Int):
    """``Fbar.m`` -- ``F(k1) = sum_j dkappa1_j Int_j / (k1 - kappa1_j)``.

    (The original comment says "rows k1, columns kappa1"; the matrix is
    actually rows kappa1, columns k1 -- the result is the same.)
    """
    k1, kappa1 = row(k1), row(kappa1)
    Cauchy = 1 / (k1[None, :] - kappa1[:, None])
    return (row(dkappa1) * row(Int)) @ Cauchy


def kappa1_quad(omega, U, delta, N):
    """``kappa1_Quad.m`` -- Gauss-Legendre nodes ``kappa1`` (and weights) on
    ``[omega/U(delta), omega/U(0)]``, and the matching heights ``y2`` with
    ``U(y2) = omega / kappa1`` (by linear interpolation on a fine grid).

    Weights are negative (``a > b`` in :func:`lgwt`), exactly as in MATLAB.
    Requires ``U(0) != 0`` (slip), otherwise ``omega/U(0)`` is infinite.
    """
    a = omega / U.f(0)
    b = omega / U.f(delta)
    kappa1, dkappa1 = lgwt(N, a, b)
    kappa1, dkappa1 = kappa1[::-1], dkappa1[::-1]
    y2_prior = np.linspace(0, delta, 10 * N)
    kappa1_prior = omega / U.f(y2_prior)
    order = np.argsort(kappa1_prior)
    y2 = np.interp(kappa1, kappa1_prior[order], y2_prior[order])
    return kappa1, dkappa1, y2


def integrand_distributed(omega, U, c, I0, omega_distributed, Km_distributed, y2):
    """``Integrand_Distributed.m``.

    .. warning:: **Bug fixed.**  The MATLAB used ``y2`` without it being an
       argument (it would error).  ``y2`` is now the final argument.
       Assumes a constant speed of sound (original note).
    """
    return -omega * c.f(0) ** 2 * I0 * omega_distributed / (Km_distributed * U.f(y2) ** 2)


def distributed_wall_pressure(y2, kappa1, shift_y2, omega, k3, U, c, delta, N_bl, l1, l2, HP):
    """``Distributed_WallPressure.m`` -- wall pressure from a unit source at each
    ``y2[j]`` (forcing ``kappa1[j]``), wall BC ``l1 phi' + l2 phi = 0``.

    ``shift_y2`` evaluates just above the source to stay off the critical
    layer.  ``HP`` picks straight (0) or semicircular (+-1) paths.
    Slow in MATLAB ("we really only want to run it once"); vectorised here.
    """
    y2, kappa1 = row(y2), row(kappa1)
    ys = y2 + shift_y2
    dec_y2 = comp_dec_aug18a_skip(kappa1, k3, omega, U, c, delta, N_bl, path_to(delta, ys, HP))
    dec_0 = comp_dec_aug18a_skip(kappa1, k3, omega, U, c, delta, N_bl, path_to(delta, 0.0, HP))
    C_0 = 1j * (omega - U.f(0) * kappa1)
    P_0 = l1(kappa1)  # "Just cancelling out the C_0^3"
    Rat = 1 / (C_0 * c.f(0) ** 2)  # "Absorbing y2 into p, rather than phi"
    return Rat * P_0 * dec_y2.p / (l1(kappa1) * dec_0.dphi + l2(kappa1) * dec_0.phi)


def wall_pressure(y2, omega, k3, U, c, delta, l1, l2, N_bl, HP):
    """``WallPressure.m`` -- earlier attempt at the same quantity.

    .. warning:: **Suspect in the original.**  ``Q0`` contains
       ``(omega - U(y2) kappa1)^4`` with ``kappa1 = omega/U(y2)``, which is
       identically zero, so the result is always 0 (or NaN).  Ported
       literally; :func:`distributed_wall_pressure` supersedes it.
    """
    y2 = row(y2)
    kappa1 = omega / U.f(y2)
    Q0 = (c.f(y2) ** 2 * (omega - U.f(y2) * kappa1) ** 4) / (c.f(0) ** 2 * (omega - U.f(0) * kappa1) ** 4)
    phid = y_continuous_scattering([0, 0], kappa1, k3, omega, U, c, delta, N_bl, np.sign(HP),
                                   truncate=HP != 0)
    D_l = l1(kappa1) * phid.dphi[0] + l2(kappa1) * phid.phi[0]
    I0 = np.zeros(y2.size, dtype=complex)
    for j in range(y2.size):
        Y = y_greens_generic(0, y2[j], kappa1[j], k3, omega, U, c, l1, l2, delta, N_bl, HP)
        I0[j] = Q0[j] * Y.p[0, 0] / D_l[j]
    return I0


def upstream_distributed_a(y2, kappa1, shift_y2, omega, k3, U, c, delta, N_bl, l1, l2, HP) -> Struct:
    """``Upstream_Distributed_A.m`` -- the ``x2``-independent part of the
    upstream Green's function for every source height ``y2[j]``.

    Stores full profiles of the decaying (above) and wall (below) solutions
    so :func:`upstream_distributed_b` can interpolate them onto any grid.

    Additionally stores ``dec_y2.v`` / ``generic_y2.v`` (see
    :func:`upstream_distributed_b` for why).
    """
    y2, kappa1 = row(y2), row(kappa1)
    dec_0 = comp_dec_aug18a_skip(kappa1, k3, omega, U, c, delta, N_bl, path_to(delta, 0.0, HP))
    DZ = l1(kappa1) * dec_0.dphi + l2(kappa1) * dec_0.phi
    C_p = 1j * (omega - U.f(y2 + shift_y2) * kappa1)
    C_m = 1j * (omega - U.f(y2 - shift_y2) * kappa1)
    C_0 = 1j * (omega - U.f(0) * kappa1)
    base = c.f(0) ** 2 * C_0 ** 4 * DZ
    # Values at the source, just above/below it (always straight paths, as in
    # the original where the semicircle variant was commented out).
    dec_y2 = comp_dec_aug18a_skip(kappa1, k3, omega, U, c, delta, N_bl,
                                  straight_path(delta, y2 + shift_y2))
    gen_y2 = comp_generic_aug18a_skip(kappa1, k3, omega, U, c, l1, l2, N_bl,
                                      straight_path(0.0, y2 - shift_y2))
    dec = comp_dec_aug18a(kappa1, k3, omega, U, c, delta, N_bl, straight_path(delta, 0.0))
    generic = comp_generic_aug18a(kappa1, k3, omega, U, c, l1, l2, N_bl, straight_path(0.0, delta))
    return Struct(dec=dec, generic=generic, DZ=DZ, dec_y2=dec_y2, generic_y2=gen_y2,
                  Q0_p=-C_p ** 3 / base, Q0_m=-C_m ** 3 / base, Q0=1 / base,
                  y2=y2, kappa1=kappa1, omega=omega, U=U, c=c)


def upstream_distributed_b(x2, A: Struct, omega, k3, U, c, delta) -> Struct:
    """``Upstream_Distributed_B.m`` -- Green's function ``G[x2, y2]`` by
    (pchip) interpolation of the profiles from :func:`upstream_distributed_a`.

    Assumes ``x2`` ascending (as the MATLAB concatenation did).  At the end
    ``phi``/``dphi`` are re-derived from ``p``/``v`` (continuity is imposed on
    pressure, so this keeps the critical-layer singularity well understood).

    .. warning:: **Bugs fixed.**  ``numel(x2_y2 == 1)`` was meant to be
       ``numel(x2_y2) == 1``; and the ``x2 == y2`` branch read
       ``dec_y2.v``, which ``_A`` never stored (it would crash whenever a grid
       point hit a source height exactly).  ``_A`` now stores it.
    """
    x2 = row(x2).astype(float)
    y2 = A.y2
    G = Struct(**{f: np.zeros((x2.size, y2.size), dtype=complex) for f in FIELDS})
    gamma = gamma_ff(k3, omega, np.real(c.f(delta)), np.real(U.f(delta)))[0](A.kappa1)
    zg = np.real(A.generic.z)
    zd = np.real(A.dec.z)
    for j in range(y2.size):
        ma = x2 < y2[j]
        my = x2 == y2[j]
        mb = (x2 > y2[j]) & (x2 < delta)
        mc = x2 >= delta
        Qd = A.Q0[j] * A.dec_y2.p[j]
        Qg = A.Q0[j] * A.generic_y2.p[j]
        for f in FIELDS:
            G[f][ma, j] = Qd * pchip(zg, A.generic[f][j], x2[ma])
            G[f][mb, j] = Qg * pchip(zd, A.dec[f][j], x2[mb])
        if my.any():
            G.phi[my, j] = Qg * A.dec_y2.phi[j]
            G.dphi[my, j] = 0  # "Ill-defined, a pain."
            G.p[my, j] = Qg * A.dec_y2.p[j]
            G.v[my, j] = Qg * A.dec_y2.v[j]
        E = np.exp(-gamma[j] * x2[mc])
        G.phi[mc, j] = A.Q0_m[j] * A.generic_y2.phi[j] * E
        G.dphi[mc, j] = -gamma[j] * G.phi[mc, j]
        Om = omega - U.f(delta) * A.kappa1[j]
        G.p[mc, j] = 1j * Om ** 3 * G.phi[mc, j]
        G.v[mc, j] = Om ** 2 * G.dphi[mc, j]

    inner = x2 < delta
    C = np.empty((x2.size, y2.size), dtype=complex)
    dC = np.zeros((x2.size, y2.size), dtype=complex)
    C[inner] = 1j * (omega - np.outer(U.f(x2[inner]), A.kappa1))
    C[~inner] = 1j * (omega - U.f(delta) * A.kappa1)[None, :]
    dC[inner] = -1j * np.outer(U.df(x2[inner]), A.kappa1)
    G.phi = -G.p / C ** 3
    G.dphi = -(G.v + 3 * C * dC * G.phi) / C ** 2
    G.x2, G.y2, G.kappa1 = x2, y2, A.kappa1
    G.omega, G.c, G.U = A.omega, A.c, A.U
    return G


#: Allowed values of the ``weighting`` toggle (see :func:`distribution_weight`).
WEIGHTINGS = ("thesis", "jacobian")


def distribution_weight(omega, U, c, y2, weighting="thesis"):
    """Factor ``g(y2)`` converting ``omega_c(y2) d kappa1`` into source strength.

    The distributed source is ``int F(y2) (...) dy2`` with ``F = -c0^2
    omega_c / U'`` (thesis 6.3.9), rewritten as an integral over
    ``kappa1 = omega / U(y2)`` (``dy2 = -U^2 d kappa1 / (omega U')``).  Both
    the incident field and ``F_bar`` use ``g * omega_c * d kappa1``.

    * ``"thesis"`` -- ``g = omega c0^2 / U^2``: thesis eq. (6.3.10) and the
      original MATLAB (``Upstream_Distributed_Integrated.m``,
      ``Distributed_TrailingEdge.m``).  (The printed (6.3.10) also carries a
      minus sign; the MATLAB absorbs it in the sign of ``I0``.)
    * ``"jacobian"`` -- ``g = c0^2 U^2 / (omega U'^2)``: the literal change of
      variables.  It differs from ``"thesis"`` by ``(U^2 / (omega U'))^2 =
      (dy2/d kappa1)^2``, i.e. the thesis applies the Jacobian upside down.

    .. warning:: With ``"jacobian"`` the integrand, back in ``y2``, is
       ``c0^2 omega_c / U'``, which is not integrable where ``U' = 0`` -- at
       ``y2 = delta`` for the parabolic profiles -- unless ``omega_c``
       vanishes there at least linearly.  For ``omega_c = 1`` the integral
       diverges logarithmically and results grow with ``N_distribution``.
    """
    if weighting == "thesis":
        return omega * c.f(y2) ** 2 / U.f(y2) ** 2
    if weighting == "jacobian":
        return c.f(y2) ** 2 * U.f(y2) ** 2 / (omega * U.df(y2) ** 2)
    raise ValueError(f"weighting must be one of {WEIGHTINGS}, not {weighting!r}")


def upstream_distributed_integrated(x1, omega_c, G: Struct, dkappa1, weighting="thesis") -> Struct:
    """``Upstream_Distributed_Integrated.m`` -- integrate the Green's function
    against the vorticity distribution ``omega_c(y2)``.  Returns
    ``phi, dphi, p, v`` of shape ``(len(x1), len(x2))``.  (Vectorised.)

    ``weighting`` selects the ``y2 -> kappa1`` weight; see
    :func:`distribution_weight`.  ``"thesis"`` reproduces the MATLAB.
    """
    x1 = row(x1)
    w = row(omega_c) * row(dkappa1) * distribution_weight(G.omega, G.U, G.c, G.y2, weighting)
    fourier = np.exp(-1j * np.outer(G.kappa1, x1))  # (y2, x1)
    return Struct(**{f: (fourier.T * w[None, :]) @ G[f].T for f in FIELDS})


def scrappy_i0(k1, y2, I, U0=0.3, c0=5.0, delta=1.0, N_bl=6000, omega=1.0, k3=0.0):
    """``Scrappy_I0.m`` -- exploratory: ``I0(k1)`` for a source at ``y2`` over a
    hard wall, before setting ``k1 = kappa1``.  Parameters were hard-coded
    in the MATLAB (linear profile ``U0 -> 1``, ``c = 5``); exposed here."""
    from .profiles import c_const, u_linear

    U, c = u_linear(U0, 1, delta), c_const(c0)
    C0 = lambda k: 1j * (omega - U.f(0) * k)
    dC0 = lambda k: -1j * U.df(0) * k
    l1 = lambda k: -C0(k) ** 2
    l2 = lambda k: -3 * C0(k) * dC0(k)
    k1 = row(k1)
    dec_y2 = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N_bl, path_to(delta, y2, I))
    dec_0 = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N_bl, path_to(delta, 0.0, I))
    C_0 = 1j * (omega - U.f(0) * k1)
    C_y2 = 1j * (omega - U.f(y2) * k1)
    P_0 = -C_0 ** 3 * l1(k1)
    return C_y2 ** 3 / C_0 ** 4 * P_0 * dec_y2.phi / (l1(k1) * dec_0.dphi + l2(k1) * dec_0.phi)
