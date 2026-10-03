"""Integration contours in the complex wavenumber (k1) and wall-normal (x2) planes.

Three families from the MATLAB tree live here:

* ``CommonTools/Integration/*Contour*.m`` -- infinite contours ``z = C(t)``,
  ``t`` real, used as the Wiener-Hopf factorisation contour (``contour.Cp``,
  ``contour.Cm``).  Each returns the pair ``(C, dC)`` of callables.
* ``GenericInversion/Contours_Weights/*.m`` -- *discretised* contours for the
  Fourier inversion, returned as a :class:`Struct` with ``k`` (points) and
  ``wk`` (quadrature weight x dk/dt).
* Paths in the ``x2`` plane along which the Pridmore-Brown equation is
  integrated (``z.f``/``z.df`` structs) -- :func:`straight_path`,
  :func:`semicircle_path`, :func:`contour_semicircle`.  Deforming these into
  the complex ``x2`` plane is how the code steps around the critical layer
  ``U(x2) = omega/k1``.
"""
from __future__ import annotations

import numpy as np

from ._utils import Path, Struct, heaviside
from .quadrature import lgwt


# ==========================================================================
# CommonTools/Integration -- infinite factorisation contours z = C(t)
# ==========================================================================
def _sech2(u):
    """``sech(u)^2`` without overflow for large ``|u|``."""
    e = np.exp(-2 * np.abs(u))
    return 4 * e / (1 + e) ** 2


def tanh_contour(a, b, tm, tp):
    """``TanhContour.m``: ``C(t) = t + i(a + b tanh(2(t - (tp+tm)/2)/(tp-tm)))``.

    Tends to ``a +- b`` (imaginary part) as ``t -> +-inf``.
    """
    s = 2 / (tp - tm)
    mid = (tp + tm) / 2
    C = lambda t: t + 1j * (a + b * np.tanh(s * (t - mid)))
    dC = lambda t: 1 + 1j * b * s * _sech2(s * (t - mid))
    return C, dC


def tanh_contour_b(a, bm, bp, tm, tp):
    """``TanhContourB.m``: tanh contour with independent halves either side of 0.

    .. warning:: **Bug fixed.**  The original derivative of the ``t < 0`` half
       was ``(bm + a) sech^2(t/tm)/tm``; differentiating
       ``fm = a - (bm - a) tanh(t/tm)`` gives ``-(bm - a) sech^2(t/tm)/tm``.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    fp = lambda t: a + (bp - a) * np.tanh(t / tp)
    fm = lambda t: a - (bm - a) * np.tanh(t / tm)
    dfp = lambda t: (bp - a) * _sech2(t / tp) / tp
    dfm = lambda t: -(bm - a) * _sech2(t / tm) / tm
    C = lambda t: t + 1j * (heaviside(t) * fp(t) + heaviside(-t) * fm(t))
    dC = lambda t: 1 + 1j * (heaviside(t) * dfp(t) + heaviside(-t) * dfm(t))
    return C, dC


def gauss_contour(Cinf, H, W):
    """``GaussContour.m``: Gaussian bump through 0 and ``W``, max height ``H``,
    tending to ``Cinf < 0`` (imaginary part) at both ends.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    sigma2 = np.log((Cinf - H) / Cinf) / (W / 2) ** 2
    t0 = W / 2
    g = lambda t: Cinf + (H - Cinf) * np.exp(-sigma2 * (t - t0) ** 2)
    dg = lambda t: -2 * sigma2 * (H - Cinf) * (t - t0) * np.exp(-sigma2 * (t - t0) ** 2)
    return (lambda t: t + 1j * g(t)), (lambda t: 1 + 1j * dg(t))


def ellipse_contour(H, W):
    """``EllipseContour.m``: half-ellipse of height ``H`` over ``[0, W]``, else flat.

    .. warning:: **Bug fixed.**  The original had ``dCp = 1 + dg`` (missing
       the factor ``1i``) although ``Cp = t + 1i*g``.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    b = 2 * H / W
    ge = lambda t: np.sqrt(H ** 2 - b ** 2 * (t - W / 2) ** 2 + 0j)
    dge = lambda t: -b ** 2 * (t - W / 2) / ge(t)
    g = lambda t: ge(t) * heaviside(W - t)
    dg = lambda t: dge(t) * heaviside(W - t)
    C = lambda t: (t + 1j * g(t)) * heaviside(t) + t * heaviside(-t)
    dC = lambda t: (1 + 1j * dg(t)) * heaviside(t) + 1.0 * heaviside(-t)
    return C, dC


def veitch_contour(c, D, W):
    """``VeitchContour.m``: contour of Veitch & Peake (2008) crossing the real
    axis at ``c``; ``D`` sets the height, ``W`` the width."""
    B = 3 * W ** 4
    A = 4j * D * W ** 3
    C = lambda t: t + c - A * t / (B + t ** 4)
    dC = lambda t: 1 - A * (B - 3 * t ** 4) / (B + t ** 4) ** 2
    return C, dC


def veitch_contour_b(Dp, Wp, Dm, Wm):
    """``VeitchContourB.m``: Veitch-Peake contour through 0 with different
    ``(D, W)`` on either side of the imaginary axis.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    Dp = -Dp  # "so both heights above real axis" (original comment)
    Bp, Ap = 3 * Wp ** 4, 4j * Dp * Wp ** 3
    Bm, Am = 3 * Wm ** 4, 4j * Dm * Wm ** 3
    Cp = lambda t: t - Ap * t / (Bp + t ** 4)
    dCp = lambda t: 1 - Ap * (Bp - 3 * t ** 4) / (Bp + t ** 4) ** 2
    Cm = lambda t: t - Am * t / (Bm + t ** 4)
    dCm = lambda t: 1 - Am * (Bm - 3 * t ** 4) / (Bm + t ** 4) ** 2
    C = lambda t: Cp(t) * heaviside(t) + Cm(t) * heaviside(-t)
    dC = lambda t: dCp(t) * heaviside(t) + dCm(t) * heaviside(-t)
    return C, dC


def truncate(L, x, N):
    """``Truncate.m``: ``L(x)`` for ``|x| < N``, zero elsewhere (avoids NaNs)."""
    x = np.asarray(x)
    out = np.zeros(x.shape, dtype=complex)
    m = np.abs(x) < N
    out[m] = L(x[m])
    return out


def default_factorisation_contour():
    """``Factorisation/Contour_Setup.m``: identical tanh contours ``C+ = C-``."""
    Cp, dCp = tanh_contour(0, 1, -1, 1)
    Cm, dCm = tanh_contour(0, 1, -1, 1)
    return Struct(Cp=Cp, dCp=dCp, Cm=Cm, dCm=dCm)


# ==========================================================================
# GenericInversion/Contours_Weights -- discretised inversion contours
# ==========================================================================
def inf_to_finite(N):
    """``Inf_to_Finite.m``: map ``(-inf, inf)`` to ``(-1, 1)`` via
    ``k = s / sqrt(1 - s^2)`` and return Gauss-Legendre points/weights.

    Returns a Struct with ``s``, ``ws`` (GL rule), ``k`` and ``dk`` (the
    weights including the Jacobian).
    """
    s, w_s = lgwt(N, -1, 1)
    s = s[::-1]  # weights are symmetric, so not reversing w_s is harmless
    return Struct(s=s, ws=w_s, k=s / np.sqrt(1 - s ** 2),
                  dk=w_s / np.sqrt(1 - s ** 2) ** 3)


def parabola_3points(zm, z0, zp):
    """``Parabola_3Points.m``: complex parabola with ``z(-1)=zm, z(0)=z0, z(1)=zp``.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing calls it except :func:`finite_parabola_3points`, itself unused.
       Exercised only by ``tests/test_legacy.py``.
    """
    A = 0.5 * (zm + zp) - z0
    B = 0.5 * (zp - zm)
    return Struct(f=lambda t: A * t ** 2 + B * t + z0,
                  df=lambda t: 2 * A * t + B,
                  d2f=lambda t: 2 * A * np.ones(np.shape(t)))


def symmetric_hyperbola(z0r, apex, curvature, angle):
    """``Symmetric_Hyperbola.m``: complex hyperbola centred at real ``z0r``.

    ``apex`` is the extreme imaginary part (at ``Re z = z0r``); positive apex
    curves down into the lower half-plane and vice versa.  ``angle`` is the
    asymptotic argument; ``curvature`` controls how sharp the turn is.
    """
    C2 = curvature ** 2
    A2 = np.tan(angle) ** 2
    if apex >= 0:
        zi = lambda t: (apex + np.sqrt(C2)) - np.sqrt(C2 + A2 * (t - z0r) ** 2)
        dzi = lambda t: -A2 * (t - z0r) / np.sqrt(C2 + A2 * (t - z0r) ** 2)
    else:
        zi = lambda t: (apex - np.sqrt(C2)) + np.sqrt(C2 + A2 * (t - z0r) ** 2)
        dzi = lambda t: A2 * (t - z0r) / np.sqrt(C2 + A2 * (t - z0r) ** 2)
    return Struct(f=lambda t: (t + z0r) + 1j * zi(t + z0r),
                  df=lambda t: 1 + 1j * dzi(t + z0r))


def _finite(z, N):
    C_int = inf_to_finite(N)
    return Struct(k=z.f(C_int.k), wk=z.df(C_int.k) * C_int.dk)


def finite_parabola_3points(N, zm, z0, zp):
    """``Finite_Parabola_3Points.m``: discretised :func:`parabola_3points`.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    return _finite(parabola_3points(zm, z0, zp), N)


def finite_symmetric_hyperbola(N, z0r, apex, curvature, angle):
    """``Finite_Symmetric_Hyperbola.m``: discretised :func:`symmetric_hyperbola`.

    This is the "steepest-descent-like" contour hugging each acoustic branch
    cut (``Cp`` round ``kb2`` with ``-apex``, ``Cm`` round ``kb1`` with ``+apex``).
    """
    return _finite(symmetric_hyperbola(z0r, apex, curvature, angle), N)


def finite_cl_ellipse(N, a, b, shift1, shift2):
    """``Finite_CL_Ellipse.m``: two half-ellipses wrapping the critical-layer
    branch cut ``[a, b]`` (``a = omega/U_inf``, ``b = omega/U(0)``).

    ``shift1`` extends the ellipse past each end, ``shift2`` is its height.
    Returns ``(C_CLp, C_CLm)`` each with ``k``/``wk``.
    """
    A = 0.5 * (a + b)
    B = 0.5 * (b - a) + shift1
    C = 1j * shift2
    s, ws = lgwt(N, -1, 1)
    s = s[::-1]
    arg = np.pi * s / 2
    C_CLp = Struct(k=A + B * np.sin(arg) + C * np.cos(arg),
                   wk=np.pi / 2 * ws * (B * np.cos(arg) - C * np.sin(arg)))
    C_CLm = Struct(k=A - B * np.sin(arg) - C * np.cos(arg),
                   wk=np.pi / 2 * ws * (-B * np.cos(arg) + C * np.sin(arg)))
    return C_CLp, C_CLm


def pole_circle(z0, radius=0.01, step=0.01):
    """Small clockwise circle round a known pole ``z0`` (inlined in several
    MATLAB scripts as ``C_pole.k`` / ``C_pole.wk``).  Trapezium weights."""
    th = np.arange(0, 2 * np.pi + step / 2, step)
    th = th[th <= 2 * np.pi]
    return Struct(k=z0 + radius * np.exp(-1j * th),
                  wk=-step * radius * 1j * np.exp(-1j * th))


def append_contour(C: Struct, extra: Struct) -> Struct:
    """Concatenate two discretised contours (``[C.k, extra.k]`` etc.)."""
    return Struct(k=np.concatenate([C.k, extra.k]),
                  wk=np.concatenate([C.wk, extra.wk]))


def cl_ellipse_contour(kinf, k0, a, b):
    """``FixedStep/Common_Functions/CL_Ellipse_Contour.m``.

    Elliptical loop round the critical layer ``[kinf, k0]``, parametrised by
    ``t in [-1, 1]``; ``C.U``/``C.dU`` upper half, ``C.L``/``C.dL`` lower.

    .. warning:: **Unused legacy port.**  Kept for completeness of the MATLAB
       port; nothing in the package, the example scripts or the paper
       figures calls it.  Exercised only by ``tests/test_legacy.py``.
    """
    if kinf > k0:
        raise ValueError("End of critical layer must be larger than start")
    dk = 0.5 * (k0 - kinf)
    sk = 0.5 * (k0 + kinf)
    scale = 1 / np.sqrt(1 - dk ** 2 / (dk + a) ** 2)
    xf = lambda t: (a + dk) * t + sk
    xdf = lambda t: (a + dk) * np.ones(np.shape(t))
    root = lambda t: np.sqrt(1 - (xf(t) - sk) ** 2 / (dk + a) ** 2 + 0j)
    droot = lambda t: -(xf(t) - sk) / (dk + a) ** 2 / root(t)
    return Struct(U=lambda t: xf(t) + 1j * b * scale * root(t),
                  dU=lambda t: xdf(t) + 1j * b * scale * xdf(t) * droot(t),
                  L=lambda t: xf(t) - 1j * b * scale * root(t),
                  dL=lambda t: xdf(t) - 1j * b * scale * xdf(t) * droot(t))


# ==========================================================================
# Paths in the x2 plane for the Pridmore-Brown integration
# ==========================================================================
def straight_path(a, b) -> Path:
    """``z(s) = a (1 - s) + b s``: straight from ``a`` (s=0) to ``b`` (s=1).

    ``a``/``b`` may be arrays (one path per element), which lets many
    integrations run simultaneously in the vectorised solvers.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    return Path(f=lambda s: a * (1 - s) + b * s,
                df=lambda s: (b - a) * np.ones(np.shape(s)))


def semicircle_path(a, b, half_plane) -> Path:
    """Semicircle from ``a`` (s=0) to ``b`` (s=1) through the upper
    (``half_plane > 0``) or lower (``< 0``) half of the complex x2-plane.

    Generalises the several hand-written variants in the MATLAB
    (``zdy``, ``zly``, ``Contour_Semicircle`` ...), which are all instances.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    I = 1 if half_plane > 0 else -1
    m = (a + b) / 2
    r = (a - b) / 2
    sig = I * np.where(r >= 0, 1, -1)
    return Path(f=lambda s: m + r * np.exp(sig * 1j * np.pi * s),
                df=lambda s: r * sig * 1j * np.pi * np.exp(sig * 1j * np.pi * s))


def graded_path(z: Path) -> Path:
    """``z`` reparametrised by ``s -> (1 - cos(pi s))/2``: the same path, with
    fixed RK steps clustered quadratically at both ends.

    Critical points approach the ends of a boundary-layer path (the wall, or
    the edge ``x2 = delta``) as ``k1`` approaches the ends of the
    critical-layer cut; a uniform step there is too coarse long before the
    path itself becomes unsuitable.
    """
    g = lambda s: 0.5 * (1 - np.cos(np.pi * s))
    dg = lambda s: 0.5 * np.pi * np.sin(np.pi * s)
    return Path(f=lambda s: z.f(g(s)), df=lambda s: z.df(g(s)) * dg(s))


def semiellipse_path(a, b, half_plane, height) -> Path:
    """Half-ellipse from ``a`` (s=0) to ``b`` (s=1) bulging to ``height`` into
    the upper (``half_plane > 0``) or lower half-plane.  ``height = |a-b|/2``
    is :func:`semicircle_path`.  Arrays give one path per element; like the
    semicircle it leaves and meets the axis vertically.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    I = np.where(np.asarray(half_plane) > 0, 1.0, -1.0)
    m, r, H = (a + b) / 2, (a - b) / 2, np.asarray(height, dtype=float)
    return Path(f=lambda s: m + r * np.cos(np.pi * s) + 1j * I * H * np.sin(np.pi * s),
                df=lambda s: np.pi * (-r * np.sin(np.pi * s) + 1j * I * H * np.cos(np.pi * s)))


def contour_semicircle(delta, HP) -> Path:
    """``Contour_Semicircle.m``: from ``x2 = delta`` to ``0`` via a semicircle
    in the upper (``HP > 0``) or lower (``HP < 0``) half-plane."""
    if HP == 0:
        raise ValueError("Half-plane of integration unclear")
    return semicircle_path(delta, 0.0, HP)


def path_to(a, b, HP) -> Path:
    """Straight path if ``HP == 0`` else a semicircle in half-plane ``HP``.

    This is the ``if HP == 0 ... else ... end`` block repeated throughout
    ``GenericInversion``.
    """
    return straight_path(a, b) if HP == 0 else semicircle_path(a, b, HP)


def indented_path(a, b, xc, half_plane, radius) -> Path:
    """Real-axis path from ``a`` to ``b`` lifted smoothly over the point ``xc``
    into the ``half_plane`` (+1 above, -1 below):

        z(x) = x + i * half_plane * radius * B((x - xc) / radius),
        B(u) = e * exp(-1 / (1 - u^2))  for |u| < 1, else 0,

    so the path passes ``xc`` at height ``radius`` and is real outside
    ``|x - xc| < radius``.  ``B`` is C-infinity: a semicircular indentation has
    corners where it meets the axis, and fixed-step RK4 loses its order at each
    corner (errors of several % near the critical layer).

    All arguments may be arrays (one path per element); ``radius = 0`` gives a
    straight path.  Used for non-analytic (piecewise) profiles, for which the
    large semicircles of :func:`semicircle_path` are invalid; the bump stays
    inside one analytic piece as long as ``radius`` is smaller than the
    distance from ``xc`` to the nearest join.
    """
    a, b, xc, r, hp = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in
                                            (a, b, xc, radius, half_plane)))
    safe_r = np.where(r > 0, r, 1.0)

    def _bump(x):
        u = np.where(r > 0, (x - xc) / safe_r, 2.0)
        inside = np.abs(u) < 1
        w = np.where(inside, 1 - u ** 2, 1.0)
        B = np.where(inside, np.e * np.exp(-1 / w), 0.0)
        dB = np.where(inside, B * (-2 * u / w ** 2), 0.0)   # dB/du
        return B, dB

    def f(s):
        x = a + (b - a) * s
        B, _ = _bump(x)
        return x + 1j * hp * r * B

    def df(s):
        x = a + (b - a) * s
        _, dB = _bump(x)
        return (b - a) * (1 + 1j * hp * dB)

    return Path(f=f, df=df)
