"""Fourier inversion building blocks (``GenericInversion``).

Ports of ``Fourier_Nov2018.m``, ``Integral_Nov2018.m``, everything in
``PointSource_Cts`` that is a function, the ``Y_*`` functions in
``Scattering`` and ``Dirty_Dispersion.m``.

Every physical field is written as a contour integral

    q(x1, x2) = int_C  K(k) X(x1, k) Y(x2, k) dk,

discretised with points ``C.k`` and weights ``C.wk``:

* ``X = exp(-i k x1) / 2 pi``  (:func:`fourier_nov2018`),
* ``Y`` -- the wall-normal structure from the Pridmore-Brown solver
  (``Y_*`` functions; a :class:`Struct` with ``phi, dphi, p, v``, each of
  shape ``(len(x2), len(k))``),
* ``K`` -- everything independent of position (Cauchy factor, dispersion
  relation, Wiener-Hopf factors ...).

Vectorisation
-------------
The MATLAB ``Y_*`` functions looped over ``x2``, calling the boundary-layer
solver once per row.  Here all ``(x2, k1)`` pairs within the boundary layer
are integrated in a single RK4 sweep (each pair with its own path).  The
arithmetic is identical; ``max_batch`` bounds memory use.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from ._utils import FIELDS, Struct, log, row
from .contours import path_to, straight_path
from .pridmore_brown import (comp_dec_aug18a_skip, comp_generic_aug18a_skip,
                             comp_hard_aug18a_skip, gamma_ff)

#: Default upper bound on (x2 points) x (k points) per vectorised RK4 sweep.
MAX_BATCH = 400_000


# ==========================================================================
# Assembling integrals
# ==========================================================================
def fourier_nov2018(x, k):
    """``Fourier_Nov2018.m`` -- ``exp(-i x k) / 2 pi`` as a ``(len(x), len(k))`` array."""
    return np.exp(-1j * np.outer(row(x), row(k))) / (2 * np.pi)


def integral_nov2018(wk, K, X, Y):
    """``Integral_Nov2018.m`` -- ``I[x1, x2] = sum_j wk_j K_j X[x1, j] Y[x2, j]``.

    (Vectorised; the original looped over ``k`` with a waitbar.)
    """
    wk, K = row(wk), row(K)
    if wk.size != K.size:
        raise ValueError("Number of weights differs from number of points integrated over")
    X = np.reshape(X, (-1, K.size))
    Y = np.reshape(Y, (-1, K.size))
    return (X * (wk * K)[None, :]) @ Y.T


def integrate_fields(wk, K, X, Y: Struct) -> Struct:
    """:func:`integral_nov2018` applied to each of ``phi, dphi, p, v``.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    return Struct(**{f: integral_nov2018(wk, K, X, Y[f]) for f in FIELDS})


# ==========================================================================
# Vectorised boundary-layer sweeps
# ==========================================================================
def _sweep(solver, x2_targets, k1, max_batch=MAX_BATCH) -> Struct:
    """Solve for every ``(x2, k1)`` pair.  ``solver(K, Xs)`` receives flat
    arrays (``K`` k-values, ``Xs`` matching target heights) and returns a
    Struct of flat fields.  Output fields have shape ``(len(x2), len(k1))``."""
    x2_targets, k1 = row(x2_targets), row(k1)
    n, m = x2_targets.size, k1.size
    out = Struct(**{f: np.zeros((n, m), dtype=complex) for f in FIELDS})
    if n == 0:
        return out
    step = max(1, max_batch // max(m, 1))
    for i0 in range(0, n, step):
        xs = x2_targets[i0:i0 + step]
        K = np.tile(k1, xs.size)
        Xs = np.repeat(xs, m)
        res = solver(K, Xs)
        for f in FIELDS:
            out[f][i0:i0 + xs.size] = np.reshape(res[f], (xs.size, m))
    return out


def _outer(Y, x2, k1, k3, omega, U, c, delta, amp):
    """Fill rows with ``x2 >= delta`` by the free-stream solution ``amp exp(-gamma x2)``.

    As in the original, ``v = -C^2 phi'`` *omits* the ``-3 C C' phi`` term,
    i.e. assumes ``U'(delta) = 0`` (true for the parabolic profiles, **not**
    for ``u_linear``).
    """
    out = np.real(x2) >= delta
    if not out.any():
        return
    gammainf, _, _ = gamma_ff(k3, omega, np.real(c.f(delta)), np.real(U.f(delta)))
    G = gammainf(k1)
    C_inf = 1j * (omega - U.f(delta) * k1)
    E = np.exp(-np.outer(x2[out], G))
    Y.phi[out] = amp * E
    Y.dphi[out] = -amp * G * E
    Y.p[out] = -C_inf ** 3 * Y.phi[out]
    Y.v[out] = -C_inf ** 2 * Y.dphi[out]


def _truncate(Y, tol=5 * np.finfo(float).eps):
    for f in FIELDS:
        Y[f][np.abs(Y[f]) < tol] = 0


# ==========================================================================
# Scattering: decaying solution only (Y_Continuous_Scattering_*.m)
# ==========================================================================
def y_continuous_scattering(x2, k1, k3, omega, U, c, delta, N_bl, HP=0,
                            truncate=False, max_batch=MAX_BATCH) -> Struct:
    """Decaying solution ``phi_d(x2)`` for all ``x2``, ``k1``.

    Merges the four MATLAB variants:

    ==========================================  ======  ==========
    MATLAB                                      ``HP``  truncate
    ==========================================  ======  ==========
    ``Y_Continuous_Scattering_Nov2018``         0       False
    ``Y_Continuous_Scattering_Nov2018_Mod``     0       False [#]_
    ``Y_Continuous_Scattering_CL_Nov2018``      +-1     False
    ``Y_Continuous_Scattering_CL_Nov2018_Mod``  +-1     True
    ==========================================  ======  ==========

    .. [#] The ``_Mod`` straight-path version's header says it truncates tiny
       values, but the truncation lines are commented out in the MATLAB.

    ``HP = +-1`` integrates from ``delta`` to ``x2`` along a semicircle in the
    upper/lower half-plane (for ``k1`` near the critical layer).
    """
    x2, k1 = row(x2), row(k1)
    inbl = np.real(x2) < delta

    def solver(K, Xs):
        return comp_dec_aug18a_skip(K, k3, omega, U, c, delta, N_bl, path_to(delta, Xs, HP))

    Y = Struct(**{f: np.zeros((x2.size, k1.size), dtype=complex) for f in FIELDS})
    sub = _sweep(solver, x2[inbl], k1, max_batch)
    for f in FIELDS:
        Y[f][inbl] = sub[f]
    if truncate:
        _truncate(Y)
    _outer(Y, x2, k1, k3, omega, U, c, delta, 1.0)
    return Y


def y_continuous_scattering_nov2018(x2, k1, k3, omega, U, c, delta, N_bl):
    return y_continuous_scattering(x2, k1, k3, omega, U, c, delta, N_bl, 0, False)


y_continuous_scattering_nov2018_mod = y_continuous_scattering_nov2018


def y_continuous_scattering_cl_nov2018(x2, k1, k3, omega, U, c, delta, N_bl, I):
    """``y_continuous_scattering_cl_nov2018`` (thin wrapper).

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    return y_continuous_scattering(x2, k1, k3, omega, U, c, delta, N_bl, I, False)


def y_continuous_scattering_cl_nov2018_mod(x2, k1, k3, omega, U, c, delta, N_bl, I):
    """``y_continuous_scattering_cl_nov2018_mod`` (thin wrapper).

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    return y_continuous_scattering(x2, k1, k3, omega, U, c, delta, N_bl, I, True)


def dirty_dispersion(k1, k3, omega, U, c, delta, N_bl):
    """``Dirty_Dispersion.m`` -- ``D_p = p(0)`` of the decaying solution (straight path).

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    return comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N_bl, straight_path(delta, 0.0)).p


# ==========================================================================
# Point source: dispersion relations and the k-dependent factor
# ==========================================================================
def dispersion_generic_nov2018(k1, k3, omega, U, c, l1, l2, delta, N_bl, HP):
    """``Dispersion_Generic_Nov2018.m`` -- ``D_l = l1 phi_d'(0) + l2 phi_d(0)``.

    ``HP`` picks the path from ``delta`` to 0 (0: straight, +-1 semicircle).
    Returns an array shaped like ``k1``.
    """
    dec = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N_bl, path_to(delta, 0.0, HP))
    K = row(k1)
    return np.reshape(l1(K) * dec.dphi + l2(K) * dec.phi, np.shape(k1))


def dispersion_hard_nov2018(k1, k3, omega, U, c, delta, N_bl, HP):
    """``Dispersion_Hard_Nov2018.m`` -- ``D_h = v_d(0)``."""
    dec = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N_bl, path_to(delta, 0.0, HP))
    return np.reshape(dec.v, np.shape(k1))


def _k_common(y2, k1, omega, U, c, D, phi_scale):
    k1 = row(k1)
    kappa1 = omega / U.f(y2)
    K_Cauchy = 1 / (2j * np.pi) / (k1 - kappa1)
    C = lambda x2: 1j * (omega - U.f(x2) * k1)
    W_y2 = c.f(y2) ** 2 * C(y2) ** 4
    W_0 = c.f(0) ** 2 * C(0) ** 4
    return W_y2 * phi_scale(C) / (W_0 * D) * K_Cauchy


def k_generic(y2, k1, k3, omega, U, c, l1, l2, delta, N_bl, HP):
    """``K_Generic.m`` -- position-independent factor for a point source at
    height ``y2`` over a wall with BC ``l1 phi' + l2 phi = 0``:
    Cauchy term ``1/(2 pi i (k1 - kappa1))`` x Wronskian scaling / ``D_l``."""
    D = row(dispersion_generic_nov2018(row(k1), k3, omega, U, c, l1, l2, delta, N_bl, HP))
    return _k_common(y2, k1, omega, U, c, D, lambda C: 1.0)


def k_hard(y2, k1, k3, omega, U, c, delta, N_bl, HP):
    """``K_Hard.m`` -- as :func:`k_generic` for a hard wall (``D_h = v(0)``),
    including the ``-C(0)^5`` normalisation of the hard-wall solution."""
    D = row(dispersion_hard_nov2018(row(k1), k3, omega, U, c, delta, N_bl, HP))
    return _k_common(y2, k1, omega, U, c, D, lambda C: -C(0) ** 5)


def wake_dec_phi(k1, k3, omega, Um, cm, deltam, N_bl, z_below):
    """``WakeDec_phi.m`` -- ``phi`` of the lower-side decaying solution at the wake."""
    return comp_dec_aug18a_skip(k1, k3, omega, Um, cm, deltam, N_bl, z_below).phi


def wake_dec_dphi(k1, k3, omega, Um, cm, deltam, N_bl, z_below):
    """``WakeDec_dphi.m`` -- ``phi'`` of the lower-side decaying solution at the wake."""
    return comp_dec_aug18a_skip(k1, k3, omega, Um, cm, deltam, N_bl, z_below).dphi


# ==========================================================================
# Point source: Green's-function x2-structure (Y_Greens*.m)
# ==========================================================================
def _greens(x2, y2, k1, k3, omega, U, c, delta, N_bl, HP, lower_solver, max_batch):
    """Shared body of :func:`y_greens_hard` / :func:`y_greens_generic`.

    ``lower_solver(K, path)`` gives the wall-side solution (hard or generic)
    along a path starting at the wall.  Returns ``(Y, lower_at_y2)``.
    """
    x2, k1 = row(x2), row(k1)
    low_y2 = lower_solver(k1, path_to(0.0, y2, HP)).phi
    d = comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N_bl, path_to(delta, y2, HP)).phi

    Y = Struct(**{f: np.zeros((x2.size, k1.size), dtype=complex) for f in FIELDS})
    rx = np.real(x2)
    above = (rx < delta) & (rx > y2)
    below = (rx < delta) & (rx <= y2)
    dec = _sweep(lambda K, Xs: comp_dec_aug18a_skip(K, k3, omega, U, c, delta, N_bl,
                                                   path_to(delta, Xs, HP)),
                 x2[above], k1, max_batch)
    low = _sweep(lambda K, Xs: lower_solver(K, path_to(0.0, Xs, HP)), x2[below], k1, max_batch)
    for f in FIELDS:
        Y[f][above] = low_y2 * dec[f]
        Y[f][below] = d * low[f]
    _truncate(Y)
    _outer(Y, x2, k1, k3, omega, U, c, delta, low_y2)
    return Y


def y_greens_hard(x2, y2, k1, k3, omega, U, c, delta, N_bl, HP, max_batch=MAX_BATCH):
    """``Y_GreensHard.m`` -- ``phi_>(x2) phi_<(y2)`` for a hard wall: the
    decaying solution above the source, hard-wall solution below.  The jump
    conditions are absorbed into :func:`k_hard`."""
    solver = lambda K, z: comp_hard_aug18a_skip(K, k3, omega, U, c, N_bl, z)
    return _greens(x2, y2, k1, k3, omega, U, c, delta, N_bl, HP, solver, max_batch)


def y_greens_generic(x2, y2, k1, k3, omega, U, c, l1, l2, delta, N_bl, HP, max_batch=MAX_BATCH):
    """``Y_GreensGeneric.m`` -- as :func:`y_greens_hard` with wall BC
    ``l1 phi' + l2 phi = 0``."""
    solver = lambda K, z: comp_generic_aug18a_skip(K, k3, omega, U, c, l1, l2, N_bl, z)
    return _greens(x2, y2, k1, k3, omega, U, c, delta, N_bl, HP, solver, max_batch)


def source_pressure_limit(k1, k3, omega, U, c, delta, y2, side, l1=None, l2=None,
                          rel_eps=0.02, N_min=1000):
    """One-sided pressure limit ``p(y2) = lim -C^3 phi`` at the critical level.

    ``side = -1``: the wall solution (BC ``l1 phi' + l2 phi = 0``), approached
    from below along ``[0, y2)``.  ``side = +1``: the decaying solution,
    approached from above along ``(y2, delta]``.  Neither path crosses the
    critical level, so no branch choice is involved.

    ``p`` has Frobenius exponents 0 and 3 at the critical level with no linear
    term, so ``p(y2 -+ eps) = p(y2) + O(eps^2)``; one Richardson step on
    ``eps`` and ``eps/2`` gives ``O(eps^3)``.  ``eps = rel_eps * min(y2,
    delta - y2)``, with at least 40 RK steps per ``eps``.
    """
    k = np.atleast_1d(np.asarray(k1, dtype=complex))
    eps = rel_eps * min(y2, delta - y2)
    length = y2 if side < 0 else delta - y2
    vals = []
    for e in (eps, eps / 2):
        N = max(N_min, int(np.ceil(40 * length / e)))
        if side < 0:
            r = comp_generic_aug18a_skip(k, k3, omega, U, c, l1, l2, N, straight_path(0.0, y2 - e))
        else:
            r = comp_dec_aug18a_skip(k, k3, omega, U, c, delta, N, straight_path(delta, y2 + e))
        vals.append(r.p)
    return (4 * vals[1] - vals[0]) / 3


def y_greens_convected(x2, y2, kappa1, k3, omega, U, c, l1, l2, delta, N_bl, HP=1,
                       rel_eps=0.02, max_batch=MAX_BATCH) -> Struct:
    """Regularised Green's function for a vortex sheet at ``y2`` convecting with
    the flow, ``kappa1 = omega / U(y2)`` -- a replacement for evaluating
    ``Y_Greens*`` at ``kappa1``.

    Solves ``D_R^dagger g = exp(-i kappa1 x1) delta(x2 - y2)``.  Since
    ``D_R^dagger = -c0^2 C L`` with ``L`` the ODE (4.2.4), and ``L``'s
    Green's function carries ``c0^2 C^4|_{y2} / c0^2 C^4|_0`` (Abel),

        g = -C(y2)^3 phi_<(x_<) phi_>(x_>) / (c0(0)^2 C(0)^4 D)

    (thesis 2.20 / 6.3.11 lack the minus sign; 6.3.4 and 6.3.12 have it).
    At the source ``C(y2) = 0`` and ``phi(y2) = inf``, so this is written in
    the regular form

        g = p_<(y2) phi_>(x2) / (c0(0)^2 C(0)^4 D)    x2 > y2
        g = phi_<(x2) p_>(y2) / (c0(0)^2 C(0)^4 D)    x2 < y2

    with ``p = lim -C^3 phi`` taken one-sidedly (:func:`source_pressure_limit`).
    This limit equals the causal one (``kappa1 -> kappa1 + i0``, equivalently
    ``Im omega < 0``), which tests verify.  ``phi_<`` satisfies the wall BC
    with ``phi_<(0) = l1``, ``phi_<'(0) = -l2``, and ``D = W(0) = l1 phi_d'(0) +
    l2 phi_d(0)`` (``Dispersion_Generic``).

    The Wronskian ``D(kappa1)`` needs ``phi_d`` continued *through* the source's
    own critical level.  ``HP = +1`` (passing above it) is the causal branch;
    ``HP = -1`` gives the complex conjugate, as some original scripts used.

    Returns a Struct with fields ``phi, dphi, p, v`` (length ``len(x2)``; the
    exponential ``exp(-i kappa1 x1)`` is *not* included), plus ``D``,
    ``p_below = p_<(y2)``, ``p_above = p_>(y2)`` and the wall values
    ``phi0, dphi0, p0, v0`` (``p0 = I0``, the unit-source wall pressure, 6.3.5).
    At ``x2 == y2`` exactly, ``p`` is the (continuous) limit and ``phi, dphi,
    v`` are NaN.
    """
    if HP == 0:
        raise ValueError("HP must be +1 or -1: D(kappa1) needs a path round the critical level")
    x2 = row(x2)
    k = np.atleast_1d(complex(kappa1))
    C = lambda x: 1j * (omega - U.f(x) * k)
    p_lo = source_pressure_limit(k, k3, omega, U, c, delta, y2, -1, l1, l2, rel_eps)[0]
    p_hi = source_pressure_limit(k, k3, omega, U, c, delta, y2, +1, rel_eps=rel_eps)[0]
    dec0 = comp_dec_aug18a_skip(k, k3, omega, U, c, delta, N_bl, path_to(delta, 0.0, HP))
    D = (l1(k) * dec0.dphi + l2(k) * dec0.phi)[0]
    norm = 1 / (c.f(0) ** 2 * C(0)[0] ** 4 * D)

    rx = np.real(x2)
    above = (rx > y2) & (rx < delta)
    below = rx < y2
    at = rx == y2
    Y = Struct(**{f: np.zeros(x2.size, dtype=complex) for f in FIELDS})
    dec = _sweep(lambda K, Xs: comp_dec_aug18a_skip(K, k3, omega, U, c, delta, N_bl,
                                                   straight_path(delta, Xs)), x2[above], k, max_batch)
    low = _sweep(lambda K, Xs: comp_generic_aug18a_skip(K, k3, omega, U, c, l1, l2, N_bl,
                                                       straight_path(0.0, Xs)), x2[below], k, max_batch)
    for f in FIELDS:
        Y[f][above] = p_lo * dec[f][:, 0]
        Y[f][below] = p_hi * low[f][:, 0]
    Y2 = Struct(**{f: Y[f][:, None] for f in FIELDS})
    _outer(Y2, x2, k, k3, omega, U, c, delta, p_lo)
    out = Struct(**{f: norm * Y2[f][:, 0] for f in FIELDS})
    for f in ("phi", "dphi", "v"):
        out[f][at] = np.nan
    out.p[at] = norm * p_lo * p_hi
    C0, dC0 = C(0)[0], (-1j * U.df(0) * k)[0]
    phi0, dphi0 = norm * p_hi * l1(k)[0], norm * p_hi * -l2(k)[0]
    out.update(D=D, p_below=p_lo, p_above=p_hi, norm=norm, phi0=phi0, dphi0=dphi0,
               p0=-C0 ** 3 * phi0, v0=-3 * dC0 * C0 * phi0 - C0 ** 2 * dphi0)
    return out


def y_greens_wake(x2, y2, k1, k3, omega, Up, cp, Um, cm, l1, l2, deltap, deltam, N_bl, HP,
                  max_batch=MAX_BATCH) -> Struct:
    """``Y_GreensWake.m`` -- point source above a wake (``x2 = 0`` is the
    wake centreline, flow ``Up`` above and ``Um`` below as functions of ``-x2``).

    ``l1``/``l2`` encode the matching to the lower-side solution (see
    :func:`.scripts.point_source.point_source_wake`).

    .. warning:: **Bugs fixed.**

       * The loop ran over ``jmin + 1 + (1:N_ybl)``, which is only the right
         set of rows if ``x2`` contains ``-deltam`` exactly (it did for the
         default grid); masks are used instead.
       * The lower free-stream ``phi`` used ``x2(-x2 >= deltap)`` where
         ``deltam`` was meant (a size mismatch whenever ``deltap != deltam``).
    """
    x2, k1 = row(x2), row(k1)
    rx = np.real(x2)
    Ell = comp_generic_aug18a_skip(k1, k3, omega, Up, cp, l1, l2, N_bl, path_to(0.0, y2, HP)).phi
    d = comp_dec_aug18a_skip(k1, k3, omega, Up, cp, deltap, N_bl, path_to(deltap, y2, HP)).phi

    Y = Struct(**{f: np.zeros((x2.size, k1.size), dtype=complex) for f in FIELDS})
    m_above = (rx > y2) & (rx < deltap)
    m_mid = (rx >= 0) & (rx <= y2) & (rx < deltap)
    m_below = (rx < 0) & (-rx < deltam)
    dec_p = _sweep(lambda K, Xs: comp_dec_aug18a_skip(K, k3, omega, Up, cp, deltap, N_bl,
                                                     path_to(deltap, Xs, HP)),
                   x2[m_above], k1, max_batch)
    phil = _sweep(lambda K, Xs: comp_generic_aug18a_skip(K, k3, omega, Up, cp, l1, l2, N_bl,
                                                        path_to(0.0, Xs, HP)),
                  x2[m_mid], k1, max_batch)
    dec_m = _sweep(lambda K, Xs: comp_dec_aug18a_skip(K, k3, omega, Um, cm, deltam, N_bl,
                                                     path_to(deltam, Xs, HP)),
                   -x2[m_below], k1, max_batch)
    sign_below = {"phi": 1, "dphi": -1, "p": 1, "v": -1}
    for f in FIELDS:
        Y[f][m_above] = Ell * dec_p[f]
        Y[f][m_mid] = d * phil[f]
        Y[f][m_below] = sign_below[f] * d * dec_m[f]

    _outer(Y, x2, k1, k3, omega, Up, cp, deltap, Ell)

    out = -rx >= deltam
    if out.any():
        gm, _, _ = gamma_ff(k3, omega, np.real(cm.f(deltam)), np.real(Um.f(deltam)))
        Gm = gm(k1)
        Cm = 1j * (omega - Um.f(deltam) * k1)
        E = np.exp(np.outer(x2[out], Gm))
        Y.phi[out] = d * E
        Y.dphi[out] = d * Gm * E
        Y.p[out] = -Cm ** 3 * Y.phi[out]
        Y.v[out] = -Cm ** 2 * Y.dphi[out]
    return Y


def y_boundary_generic(y2, k1, k3, omega, U, c, l1, l2, delta, N_bl, HP) -> Struct:
    """``FarField/Y_BoundaryGeneric.m`` -- ``phi_l(y2)`` with the far-field
    ``x2``-dependence, for saddle-point evaluation (source inside the BL)."""
    if y2 > delta:
        raise ValueError("This code is calibrated only for sources within the boundary layer")
    k1 = row(k1)
    Ell = comp_generic_aug18a_skip(k1, k3, omega, U, c, l1, l2, N_bl, path_to(0.0, y2, HP)).phi
    gammainf, _, _ = gamma_ff(k3, omega, np.real(c.f(delta)), np.real(U.f(delta)))
    G = gammainf(k1)
    C_inf = 1j * (omega - U.f(delta) * k1)
    Y = Struct(phi=Ell, dphi=-Ell * G)
    Y.p = -C_inf ** 3 * Y.phi
    Y.v = -C_inf ** 2 * Y.dphi
    return Y
