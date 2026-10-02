"""The Pridmore-Brown equation across the boundary layer.

Ports of ``FixedStep/Common_Functions/Gamma_FF.m`` and everything in
``FixedStep/PridmoreBrown`` (except the steepest-descent files, which are in
:mod:`.steepest_descent`).

Formulation
-----------
With ``C(x2) = i(omega - U(x2) k1)`` the code works with an auxiliary
potential ``phi`` defined so that

    p = -C^3 phi,        v = -C^2 phi' - 3 C C' phi,

which turns Pridmore-Brown into ``phi'' + P phi' + Q phi = 0`` with

    P = 2 c'/c + 4 C'/C,
    Q = 3 (C'' c^2 - 2 C' c c') / (c^2 C) - (k1^2 + k3^2) - C^2/c^2.

(``ddCc2`` in the MATLAB is ``C'' c^2 - 2 C' c c'``.)  The state vector is
``[phi, phi']`` and every solver returns a :class:`Struct` with ``phi``,
``dphi``, ``p`` and ``v``.

Solutions
---------
* **decaying** (``dec``): matches ``exp(-gamma x2)`` at the boundary-layer
  edge ``x2 = delta`` and is integrated *down* to the wall;
* **hard** (``hard``): ``p(0) = 1, v(0) = 0`` at the wall, integrated up;
* **generic** (``phil``): satisfies ``l1 phi' + l2 phi = 0`` at the wall.

The ``comp_*`` ("complex contour") versions integrate along a path
``x2 = z.f(s)``, ``s in [0, 1]`` (see :mod:`.contours`), so the path can
be bent into the complex plane round the critical layer.  They assume
``z.f(0) = delta`` (decaying) or ``z.f(0) = 0`` (hard/generic).  Paths may be
*vectorised*: ``z.f(s)`` may return one point per ``k1`` element, which lets
many ``(x2, k1)`` combinations be solved in one sweep.
"""
from __future__ import annotations

import numpy as np

from ._utils import Path, Struct, log, row
from .complex_funcs import sqrt_a
from .contours import indented_path, semiellipse_path
from .ivp import (fixed_step_ivp, fixed_step_ivp_gridx, fixed_step_ivp_gridx_skip,
                  fixed_step_ivp_skip)


# ==========================================================================
# Free-stream (outer) solution
# ==========================================================================
def gamma_ff(k3, omega, cinf, Uinf):
    """``Gamma_FF.m`` -- outer decay rate and its Wiener-Hopf factors.

    ``gamma^2 = k1^2 + k3^2 - (omega - Uinf k1)^2 / cinf^2``; the solution
    beyond the boundary layer is ``exp(-gamma x2)``.  Returns callables
    ``(gamma, gamma_plus, gamma_minus)`` with ``gamma = gamma_plus *
    gamma_minus``, ``gamma_plus`` analytic in the upper half-plane (branch
    point ``Cm``, cut at 3pi/2) and ``gamma_minus`` in the lower (cut at pi/2).
    """
    if np.imag(Uinf) != 0:
        raise ValueError("Uinf must be real")
    if np.imag(cinf) != 0:
        raise ValueError("cinf must be real")
    if abs(Uinf) > cinf:
        raise ValueError("Flow must be subsonic")
    M = Uinf / cinf
    k0 = omega / cinf
    C2 = complex(k0 ** 2 / (1 - M ** 2) ** 2 - k3 ** 2 / (1 - M ** 2))
    r = np.sqrt(C2)
    Cpm = (r, -r)
    if np.imag(Cpm[0]) == 0:
        Cp, Cm = min(Cpm, key=np.real), max(Cpm, key=np.real)
    elif np.imag(Cpm[0]) > 0:
        Cp, Cm = Cpm
    else:
        Cm, Cp = Cpm
    beta = np.sqrt(1 - M ** 2)
    shift = M * k0 / (1 - M ** 2)
    gammap = lambda k1: np.sqrt(beta) * sqrt_a((k1 + shift) - Cm, 3 * np.pi / 2)
    gammam = lambda k1: np.sqrt(beta) * sqrt_a((k1 + shift) - Cp, np.pi / 2)
    gamma = lambda k1: gammap(k1) * gammam(k1)
    return gamma, gammap, gammam


def branch_points(omega, Uinf, cinf, k3):
    """Acoustic branch points ``(kb1, kb2)`` of ``gamma`` (``kb1 > kb2``).

    Inlined in every inversion script.  **Bug fixed:** the scripts used
    ``sqrt(k0^2 + k3^2 (1 - M^2))``; the zeros of ``gamma`` (cf.
    :func:`gamma_ff`) need ``sqrt(k0^2 - k3^2 (1 - M^2))``.  Identical when
    ``k3 = 0`` (every script in the repo).
    """
    M = Uinf / cinf
    k0 = omega / cinf
    root = np.sqrt(complex(k0 ** 2 - k3 ** 2 * (1 - M ** 2))) / (1 - M ** 2)
    centre = -M / (1 - M ** 2) * k0
    return np.real_if_close(centre + root), np.real_if_close(centre - root)


# ==========================================================================
# ODE right-hand sides and initial conditions
# ==========================================================================
def _pq(x2, k1, omega, c, U, k3):
    c0 = c.f(x2)
    dc0 = c.df(x2)
    C = 1j * (omega - U.f(x2) * k1)
    dC = -1j * U.df(x2) * k1
    ddCc2 = -1j * U.d2f(x2) * c0 ** 2 * k1 - 2 * dC * c0 * dc0
    P = 2 * dc0 / c0 + 4 * dC / C
    Q = 3 * ddCc2 / (c0 ** 2 * C) - (k1 ** 2 + k3 ** 2) - C ** 2 / c0 ** 2
    return P, Q


def ode_fun_pb(x2, Phi, k1, omega, c, U, k3):
    """``ODEFun_PB_Aug18A.m`` / ``ODEFun_PB_Aug18A_GridX.m`` (merged).

    ``x2`` may be a scalar or an array broadcastable against ``k1``.
    """
    P, Q = _pq(x2, k1, omega, c, U, k3)
    return np.stack([Phi[1] * np.ones_like(Q), -Q * Phi[0] - P * Phi[1]])


def comp_ode_fun_pb(s, Phi, k1, omega, c, U, k3, z: Path):
    """``Comp_ODEFun_PB_Aug18A.m`` -- the ODE in the path parameter ``s``."""
    x2 = z.f(s)
    dz = z.df(s)
    P, Q = _pq(x2, k1, omega, c, U, k3)
    return np.stack([dz * Phi[1] * np.ones_like(Q), -dz * Q * Phi[0] - dz * P * Phi[1]])


def ic_decaying(k1, k3, omega, U, c, delta):
    """``IC_Decaying_PB_Aug18A.m`` (= ``Comp_IC_Decaying_PB_Aug18A.m``).

    ``phi = exp(-delta gamma)`` at ``x2 = delta`` with ``phi'`` chosen so ``p``
    and ``v`` are continuous with the outer solution ``exp(-gamma x2)``.
    """
    Cinf = 1j * (omega - U.f(delta) * k1)
    dCinf = -1j * U.df(delta) * k1
    gamma, _, _ = gamma_ff(k3, omega, np.real(c.f(delta)), np.real(U.f(delta)))
    g = gamma(k1)
    e = np.exp(-delta * g)
    return np.stack([e, (-g - 3 * dCinf / Cinf) * e])


def ic_hard(k1, omega, U, c):
    """``IC_Hard_PB_Aug18A.m``: ``p(0) = 1``, ``v(0) = 0``."""
    C0 = 1j * (omega - U.f(0) * k1)
    dC0 = -1j * U.df(0) * k1
    return np.stack([-1 / C0 ** 3, 3 * dC0 / C0 ** 4])


def ic_generic(k1, l1, l2):
    """``IC_Generic_PB_Aug18A.m``: ``l1 phi' + l2 phi = 0`` at the wall,
    normalised as ``phi(0) = l1``, ``phi'(0) = -l2``."""
    return np.stack([np.asarray(l1(k1)) * np.ones(np.shape(k1)),
                     -np.asarray(l2(k1)) * np.ones(np.shape(k1))])


def _pv(C, dC, phi, dphi):
    return -C ** 3 * phi, -3 * dC * C * phi - C ** 2 * dphi


# ==========================================================================
# Real-axis solvers (FixedStep/PridmoreBrown/*.m)
# ==========================================================================
def dec_aug18a(k1, k3, omega, U, c, delta, N):
    """``Dec_Aug18A.m`` -- decaying solution on ``x2 in [0, delta]`` (real axis).

    Returns ``dec`` with ``x`` (length N) and ``phi``/``dphi``/``p``/``v`` of
    shape ``(k_N, N)``; column ``i`` corresponds to ``x2 = dec.x[i]``.
    """
    k1 = row(k1)
    f = lambda t, phi, k: -ode_fun_pb(delta - t, phi, k, omega, c, U, k3)
    phi1, x1 = fixed_step_ivp(f, lambda k: ic_decaying(k, k3, omega, U, c, delta),
                              delta, N, k1)
    dec = Struct(x=x1, phi=phi1[::-1, 0, :].T, dphi=phi1[::-1, 1, :].T)
    C = 1j * (omega - np.outer(k1, U.f(dec.x)))
    dC = -1j * np.outer(k1, U.df(dec.x))
    dec.p, dec.v = _pv(C, dC, dec.phi, dec.dphi)
    return dec


def hard_aug18a(k1, k3, omega, U, c, delta, N):
    """``Hard_Aug18A.m`` -- hard-wall solution on ``[0, delta]`` (real axis)."""
    k1 = row(k1)
    f = lambda x2, phi, k: ode_fun_pb(x2, phi, k, omega, c, U, k3)
    phi2, x2 = fixed_step_ivp(f, lambda k: ic_hard(k, omega, U, c), delta, N, k1)
    hard = Struct(x=x2, phi=phi2[:, 0, :].T, dphi=phi2[:, 1, :].T)
    C = 1j * (omega - np.outer(k1, U.f(hard.x)))
    dC = -1j * np.outer(k1, U.df(hard.x))
    hard.p, hard.v = _pv(C, dC, hard.phi, hard.dphi)
    return hard


def dec_aug18a_gridx(x2_min, k1, k3, omega, U, c, delta, N):
    """``Dec_Aug18A_GridX.m`` -- decaying solution down to ``x2_min[k]`` per k1.

    .. warning:: **Bug fixed.**  The MATLAB reversed ``phi`` along the grid
       but not ``dec.x`` (and ``C``, ``dC``, ``p``, ``v`` were evaluated on the
       un-reversed ``x``), so ``phi`` and ``x`` were paired back-to-front
       unless the grid happened to be symmetric.  Here ``dec.x`` is reversed
       too, so ``dec.phi[k, i]`` is at ``x2 = dec.x[k, i]`` (ascending).
    """
    k1 = row(k1)
    f = lambda t, phi, k: -ode_fun_pb(delta - t, phi, k, omega, c, U, k3)
    phi1, X, phi1_end = fixed_step_ivp_gridx(
        f, lambda k: ic_decaying(k, k3, omega, U, c, delta), delta - np.ravel(x2_min), N, k1)
    dec = Struct(x=(delta - X)[:, ::-1], phi=phi1[::-1, 0, :].T, dphi=phi1[::-1, 1, :].T)
    C = 1j * (omega - k1[:, None] * U.f(dec.x))
    dC = -1j * k1[:, None] * U.df(dec.x)
    dec.p, dec.v = _pv(C, dC, dec.phi, dec.dphi)
    dec.phi_end, dec.dphi_end = phi1_end[0], phi1_end[1]
    return dec


def dec_aug18a_gridx_skip(x2_min, k1, k3, omega, U, c, delta, N):
    """``Dec_Aug18A_GridX_skip.m`` -- decaying solution evaluated at ``x2_min``.

    Output fields have the shape of ``k1``.
    """
    shape = np.shape(k1)
    K1 = row(k1)
    X2_min = row(x2_min)
    f = lambda t, phi, k: -ode_fun_pb(delta - t, phi, k, omega, c, U, k3)
    phi = fixed_step_ivp_gridx_skip(f, lambda k: ic_decaying(k, k3, omega, U, c, delta),
                                    delta - X2_min, N, K1)
    C = 1j * (omega - K1 * U.f(X2_min))
    dC = -1j * K1 * U.df(X2_min)
    p, v = _pv(C, dC, phi[0], phi[1])
    return Struct(phi=phi[0].reshape(shape), dphi=phi[1].reshape(shape),
                  p=p.reshape(shape), v=v.reshape(shape))


def hard_aug18a_gridx(x2_max, k1, k3, omega, U, c, N):
    """``Hard_Aug18A_GridX.m`` -- hard-wall solution up to ``x2_max[k]`` per k1."""
    k1 = row(k1)
    x2_max = row(x2_max)
    f = lambda x2, phi, k: ode_fun_pb(x2, phi, k, omega, c, U, k3)
    phi2, X, phi2_end = fixed_step_ivp_gridx(f, lambda k: ic_hard(k, omega, U, c), x2_max, N, k1)
    hard = Struct(x=X, phi=phi2[:, 0, :].T, dphi=phi2[:, 1, :].T)
    C = 1j * (omega - k1[:, None] * U.f(X))
    dC = -1j * k1[:, None] * U.df(X)
    hard.p, hard.v = _pv(C, dC, hard.phi, hard.dphi)
    hard.phi_end, hard.dphi_end = phi2_end[0], phi2_end[1]
    C_end = 1j * (omega - k1 * U.f(x2_max))
    dC_end = -1j * k1 * U.df(x2_max)
    hard.p_end, hard.v_end = _pv(C_end, dC_end, phi2_end[0], phi2_end[1])
    return hard


def hard_aug18a_gridx_skip(x2_max, k1, k3, omega, U, c, delta, N):
    """``Hard_Aug18A_GridX_skip.m`` -- hard-wall solution evaluated at ``x2_max``.

    ``delta`` is unused (it was ``~`` in MATLAB) but kept for call compatibility.
    """
    shape = np.shape(k1)
    K1 = row(k1)
    X2_max = row(x2_max)
    f = lambda x2, phi, k: ode_fun_pb(x2, phi, k, omega, c, U, k3)
    phi = fixed_step_ivp_gridx_skip(f, lambda k: ic_hard(k, omega, U, c), X2_max, N, K1)
    C = 1j * (omega - K1 * U.f(X2_max))
    dC = -1j * K1 * U.df(X2_max)
    p, v = _pv(C, dC, phi[0], phi[1])
    return Struct(phi=phi[0].reshape(shape), dphi=phi[1].reshape(shape),
                  p=p.reshape(shape), v=v.reshape(shape))


# ==========================================================================
# Complex-path solvers (FixedStep/PridmoreBrown/Complex_Contour/*.m)
# ==========================================================================
def critical_point(U, omega, k1, lo, hi, n_grid=2001, newton=4):
    """Complex critical point ``x_c`` with ``U(x_c) = omega / k1`` near
    ``[lo, hi]`` (one per ``k1``; NaN where ``Re(omega/k1)`` is outside the
    range of ``U``).  Assumes ``U`` monotonic on the interval: a real estimate
    by interpolation, then Newton steps in the complex plane."""
    k1 = np.atleast_1d(np.asarray(k1, dtype=complex))
    target = omega / k1
    xg = np.linspace(min(lo, hi), max(lo, hi), n_grid)
    Ug = np.real(U.f(xg))
    order = np.argsort(Ug)
    inside = (np.real(target) > Ug.min()) & (np.real(target) < Ug.max())
    x = np.interp(np.real(target), Ug[order], xg[order]).astype(complex)
    for _ in range(newton):
        x = x - (U.f(x) - target) / U.df(x)
    return np.where(inside, x, np.nan + 0j)


def _indent_for_nonanalytic(U, z: Path, k1, omega, radius=0.05, clearance=0.05):
    """Replace a complex (semicircular) path by a real-axis path between the
    same end points, lifted smoothly over each element's critical point on the
    same side as the original path (:func:`.contours.indented_path`).

    The lift is only applied where the critical point is within ``clearance``
    of the real segment or on the wrong side of it; its radius grows if the critical point is
    further off the axis, but is capped to stay inside one analytic piece
    (``U.breaks``) and away from the end points.
    """
    ends = [np.real(z.f(si)) for si in (0.0, 1.0)]
    hp = np.sign(np.imag(z.f(0.5)))
    hp = np.where(hp == 0, 1.0, hp)
    a, b = np.broadcast_arrays(ends[0], ends[1], np.empty(np.size(k1)))[:2]
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    xc = critical_point(U, omega, k1, float(np.min(lo)), float(np.max(hi)))
    xr, xi = np.real(xc), np.imag(xc)
    # lift only where the straight segment would pass within `clearance` of the
    # critical point, or on the wrong side of it (squeezing a bump against an end
    # point for a critical point already safely off the axis loses accuracy)
    near = np.isfinite(xr) & (xr > lo) & (xr < hi) & (hp * xi > -clearance)
    gaps = [np.abs(xr - lo), np.abs(hi - xr)] + [np.abs(xr - br) for br in getattr(U, "breaks", ())]
    r_cap = 0.45 * np.min(np.vstack(gaps), axis=0)
    r = np.minimum(np.maximum(radius, 1.25 * np.maximum(hp * np.nan_to_num(xi), 0)), r_cap)
    bad = near & (r <= hp * xi)
    if bad.any():
        log.warning("%d critical point(s) cannot be passed on the requested side within one "
                    "analytic piece; results there are unreliable", int(bad.sum()))
    r = np.where(near, r, 0.0)
    return indented_path(a, b, np.nan_to_num(xr), hp, r)


def spurious_critical_points(U, omega, k1, a, b, half_plane, height, eps=0.02, nb=200):
    """Number of roots of ``U(z) = omega/k1`` strictly inside the region swept by
    the half-ellipse from ``a`` to ``b`` of the given ``height`` (excluding a
    strip ``|Im z| < eps`` round the real axis, where the physical critical
    point lives).  Argument principle on the region boundary; vectorised over
    elements (``k1``, ``a``, ``b``, ``height`` broadcast together)."""
    k1, a, b, height = np.broadcast_arrays(np.asarray(k1, dtype=complex), *(np.asarray(v, dtype=float)
                                                                            for v in (a, b, height)))
    I = 1.0 if half_plane > 0 else -1.0
    m, r = (a + b) / 2, np.abs(a - b) / 2
    ok = height > eps
    yb = np.where(ok, eps, 0.5 * height)                    # inner boundary height
    th0 = np.arcsin(np.clip(yb / np.where(height > 0, height, 1), 0, 1))
    t = np.linspace(0, 1, nb)[:, None]
    th = th0 + (np.pi - 2 * th0) * t                         # arc, right to left
    arc = m + r * np.cos(th) + 1j * I * height * np.sin(th)
    darc = (np.pi - 2 * th0) * (-r * np.sin(th) + 1j * I * height * np.cos(th))
    xl, xr = m - r * np.cos(th0), m + r * np.cos(th0)
    seg = xl + (xr - xl) * t + 1j * I * yb                   # inner segment, left to right
    dseg = (xr - xl) * np.ones_like(t)
    u = omega / k1
    g = lambda z: U.df(z) / (U.f(z) - u)
    tw = np.full(nb, 1.0 / (nb - 1)); tw[[0, -1]] *= 0.5     # trapezium weights
    total = (tw[:, None] * (g(seg) * dseg + g(arc) * darc)).sum(axis=0) / (2j * np.pi)
    total = np.where(np.isfinite(total), total, 0)          # k1 = 0: no critical point at all
    return np.where(ok, np.round(I * total.real).astype(int), 0)   # orientation: CCW for I = +1


def _safe_ellipse(U, z: Path, k1, omega, fractions=(1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.12)):
    """Semi-ellipse with the same end points and side as the semicircle ``z``,
    as tall as possible (up to the semicircle itself) while staying clear of
    every spurious critical point (none within a 25 %-taller ellipse).  Identical to ``z`` for profiles without such
    points (parabolic, cubic, power law)."""
    a, b = (np.real(z.f(si)) for si in (0.0, 1.0))
    I = 1.0 if np.all(np.imag(z.f(0.5)) >= 0) else -1.0
    a, b, _ = np.broadcast_arrays(a, b, np.empty(np.size(k1)))
    full = np.abs(a - b) / 2
    H = np.zeros(np.size(k1))
    todo = np.ones(np.size(k1), bool)
    for frac in fractions:
        if not todo.any():
            break
        # safety margin: no spurious root within a 25 % taller (+0.03) ellipse, so
        # the chosen path also stays clear of roots lying just outside it
        n = spurious_critical_points(U, omega, np.ravel(k1)[todo], a[todo], b[todo], I,
                                     1.25 * frac * full[todo] + 0.03)
        good = n == 0
        idx = np.where(todo)[0][good]
        H[idx] = (frac * full)[idx]
        todo[idx] = False
    if todo.any():
        log.warning("%d path(s): spurious critical points even inside a flat ellipse; using it anyway",
                    int(todo.sum()))
        H[todo] = fractions[-1] * full[todo]
    if np.all(H == full):
        return z                                             # nothing to change
    return semiellipse_path(a, b, I, H)


def _is_complex_path(z: Path) -> bool:
    """Whether the path leaves the real axis (scalar samples, so per-element
    vectorised paths work too)."""
    return any(np.any(np.abs(np.imag(z.f(si))) > 1e-12) for si in (0.25, 0.5, 0.75))


def _prepare_path(U, z: Path, k1, omega) -> Path:
    """Make a complex (semicircular) integration path safe for profile ``U``.

    The original code deforms the boundary-layer integration onto a semicircle
    of radius delta/2 to step round the critical layer.  That is only
    equivalent to the real path if no *other* root of ``U(z) = omega/k1``
    lies between them:

    * analytic profiles -- true for parabolic, cubic and power-law profiles,
      but not for higher-degree polynomials (a degree-8 fit to the paper's
      profiles has spurious critical points near 0.39 +- 0.43i, giving O(1)
      errors).  The semicircle is flattened into a semi-ellipse just enough to
      exclude them (:func:`_safe_ellipse`); otherwise it is left unchanged.
    * piecewise profiles -- never valid; replaced by a real-axis path lifted
      locally over the critical point (:func:`_indent_for_nonanalytic`).
    """
    if not _is_complex_path(z):
        return z
    if getattr(U, "analytic", True):
        return _safe_ellipse(U, z, k1, omega)
    return _indent_for_nonanalytic(U, z, k1, omega)


def _comp_full(k1, k3, omega, U, c, N, z, ic):
    if not getattr(U, "analytic", True) and _is_complex_path(z):
        raise NotImplementedError("full-profile output along a complex path needs an analytic "
                                  "profile; use the _skip solvers or profiles.fit_analytic")
    k1 = row(k1)
    f = lambda s, phi, k: comp_ode_fun_pb(s, phi, k, omega, c, U, k3, z)
    phi1, s = fixed_step_ivp(f, ic, 1.0, N, k1)
    out = Struct(s=s, z=z.f(s), phi=phi1[:, 0, :].T, dphi=phi1[:, 1, :].T)
    C = 1j * (omega - np.outer(k1, U.f(out.z)))
    dC = -1j * np.outer(k1, U.df(out.z))
    out.p, out.v = _pv(C, dC, out.phi, out.dphi)
    return out


def _comp_skip(k1, k3, omega, U, c, N, z, ic):
    k1 = row(k1)
    z = _prepare_path(U, z, k1, omega)
    f = lambda s, phi, k: comp_ode_fun_pb(s, phi, k, omega, c, U, k3, z)
    phi1 = fixed_step_ivp_skip(f, ic, 1.0, N, k1)
    z1 = z.f(1.0)
    C = 1j * (omega - k1 * U.f(z1))
    dC = -1j * k1 * U.df(z1)
    p, v = _pv(C, dC, phi1[0], phi1[1])
    return Struct(phi=phi1[0], dphi=phi1[1], p=p, v=v)


def comp_dec_aug18a(k1, k3, omega, U, c, delta, N, z: Path):
    """``Comp_Dec_Aug18A.m`` -- decaying solution along path ``z`` (``z.f(0)=delta``).

    Returns ``s``, ``z`` and fields of shape ``(k_N, N)``.
    """
    return _comp_full(k1, k3, omega, U, c, N, z,
                      lambda k: ic_decaying(k, k3, omega, U, c, delta))


def comp_dec_aug18a_skip(k1, k3, omega, U, c, delta, N, z: Path):
    """``Comp_Dec_Aug18A_skip.m`` -- decaying solution at the path end ``z.f(1)``.

    This is *the* workhorse of the inversion code: with ``z`` running from
    ``delta`` to ``0`` it gives the dispersion relations ``D_h = v(0)`` and
    ``D_p = p(0)``.
    """
    return _comp_skip(k1, k3, omega, U, c, N, z,
                      lambda k: ic_decaying(k, k3, omega, U, c, delta))


def comp_generic_aug18a(k1, k3, omega, U, c, l1, l2, N, z: Path):
    """``Comp_Generic_Aug18A.m`` -- wall-BC solution along ``z`` (``z.f(0)=0``)."""
    return _comp_full(k1, k3, omega, U, c, N, z, lambda k: ic_generic(k, l1, l2))


def comp_generic_aug18a_skip(k1, k3, omega, U, c, l1, l2, N, z: Path):
    """``Comp_Generic_Aug18A_skip.m`` -- as above, end point only."""
    return _comp_skip(k1, k3, omega, U, c, N, z, lambda k: ic_generic(k, l1, l2))


def comp_hard_aug18a_skip(k1, k3, omega, U, c, N, z: Path):
    """``Comp_Hard_Aug18A_skip.m`` -- hard-wall solution at ``z.f(1)``."""
    return _comp_skip(k1, k3, omega, U, c, N, z, lambda k: ic_hard(k, omega, U, c))


# ==========================================================================
# NumericalAccuracy.m
# ==========================================================================
def numerical_accuracy(k1=(-1, 0, 1j, 3, 10, 1.5 + 0.1j, 1.5 + 0.05j, 1.5 + 0.01j),
                       k3=0.0, omega=1.0, delta=1.0, c=None, U=None,
                       N_count=100, N_logmin=0.3, N_logmax=4, plot=False):
    """``NumericalAccuracy.m`` -- convergence of ``D_h = v(0)`` with step count.

    Returns ``(N, err)`` where ``err[i, j] = |D_h(k1_i, N_j) - D_h(k1_i, N_inf)|``
    and ``N_inf = 10^(N_logmax + 1)``.  (Slow at the default settings.)
    """
    from .profiles import c_const, u_parabolic

    c = c or c_const(2)
    U = U or u_parabolic(0.5, 1, delta)
    k1 = np.asarray(k1, dtype=complex)
    N = np.round(np.logspace(N_logmin, N_logmax, N_count)).astype(int)
    N = np.maximum(N, 2)
    D = np.zeros((k1.size, N.size), dtype=complex)
    D_max = np.zeros(k1.size, dtype=complex)
    for i, k in enumerate(k1):
        for j, n in enumerate(N):
            D[i, j] = dec_aug18a(k, k3, omega, U, c, delta, n).v[0, 0]
        D_max[i] = dec_aug18a(k, k3, omega, U, c, delta, int(10 ** (N_logmax + 1))).v[0, 0]
    err = np.abs(D - D_max[:, None])
    if plot:
        import matplotlib.pyplot as plt

        for i, k in enumerate(k1):
            plt.loglog(N, err[i], lw=2, label=f"$k_1 = {k}$")
        plt.xlabel(r"Integration points, $N = 1/(\delta x_2)$")
        plt.ylabel(r"$|D_h(k_1,N) - D_h(k_1,N_\infty)|$")
        plt.legend()
    return N, err
