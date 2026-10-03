"""Numerical Wiener-Hopf factorisation.

Ports of ``Factorisation/Fact_*.m``, ``Prec_Setup.m`` and the legacy
``CommonTools/Integration/AddFactB.m``.

Conventions
-----------
``contour`` is a :class:`Struct` with ``Cp, dCp, Cm, dCm`` callables
(see :func:`.contours.default_factorisation_contour`): ``F+`` is the Cauchy
integral along ``Cp`` (analytic *above* it) and ``F-`` along ``Cm``
(analytic *below*).  ``prec`` has ``N`` (number of quadrature points) and, for
the non-``Finite`` versions, ``min``/``max`` (truncation of the real
parameter ``t``).

Every routine assumes vector input ``k`` so thousands of evaluation points are
handled as one matrix product.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from ._utils import Struct, heaviside, log
from .complex_funcs import log_a, log_curved_branch
from .quadrature import lgwt


def prec_setup(N=1000, t_max=20, t_min=-20) -> Struct:
    """``Prec_Setup.m`` default precision struct."""
    return Struct(N=N, max=t_max, min=t_min)


def _cauchy_sum(Fz, dz, z, k):
    """``(1/2 pi i) sum_j F(z_j) dz_j / (z_j - k)``, returned with the shape of ``k``."""
    k = np.asarray(k)
    K = np.ravel(k)
    out = (Fz * dz) @ (1 / (2j * np.pi) / (z[:, None] - K[None, :]))
    return out.reshape(k.shape)


def fact_additive(F: Callable, contour: Struct, prec: Struct):
    """``Fact_Additive.m`` -- ``F = F+ + F-`` for ``F`` decaying at infinity.

    Gauss-Legendre on the truncated interval ``t in [prec.min, prec.max]``.
    Returns callables ``(Fp, Fm)``.
    """
    t, wt = lgwt(prec.N, prec.min, prec.max)
    t, wt = t[::-1], wt[::-1]
    zp, dzp = contour.Cp(t), wt * contour.dCp(t)
    zm, dzm = contour.Cm(t), wt * contour.dCm(t)
    Fzp, Fzm = F(zp), F(zm)  # evaluated once, as MATLAB's closures re-evaluate
    Fp = lambda k: _cauchy_sum(Fzp, dzp, zp, k)
    Fm = lambda k: -_cauchy_sum(Fzm, dzm, zm, k)
    return Fp, Fm


def fact_additive_finite(F: Callable, contour: Struct, prec: Struct):
    """``Fact_Additive_Finite.m`` -- as :func:`fact_additive` over the whole
    line, via ``t = s / (1 - s^2)``, ``s in (-1, 1)``.  NaNs in ``F`` are
    zeroed.  Only ``prec.N`` is used.
    """
    s, ws = lgwt(prec.N, -1, 1)
    s, ws = s[::-1], ws[::-1]
    t = s / (1 - s ** 2)
    wt = ws * (1 + s ** 2) / (1 - s ** 2) ** 2
    zp, dzp = contour.Cp(t), wt * contour.dCp(t)
    zm, dzm = contour.Cm(t), wt * contour.dCm(t)

    def no_nan(z):
        v = np.array(F(z), dtype=complex)
        v[np.isnan(v)] = 0
        return v

    Fzp, Fzm = no_nan(zp), no_nan(zm)
    Fp = lambda k: _cauchy_sum(Fzp, dzp, zp, k)
    Fm = lambda k: -_cauchy_sum(Fzm, dzm, zm, k)
    return Fp, Fm


def _log(G, branch_cut):
    # Curved-cut support was added to the Finite version only in Nov 2018
    # ("This functionality should be added to the others"); here all get it.
    if callable(branch_cut):
        return lambda z: log_curved_branch(G(z), branch_cut)
    return lambda z: log_a(G(z), branch_cut)


def fact_mult(G: Callable, contour: Struct, prec: Struct, branch_cut):
    """``Fact_Mult.m`` -- ``G = G+ G-`` for ``G -> 1`` at infinity, via
    ``log G``.  Assumes ``log G`` does not wind round 0 or cross the cut.
    """
    Fp, Fm = fact_additive(_log(G, branch_cut), contour, prec)
    return (lambda k: np.exp(Fp(k))), (lambda k: np.exp(Fm(k)))


def fact_mult_finite(G: Callable, contour: Struct, prec: Struct, branch_cut):
    """``Fact_Mult_Finite.m`` -- :func:`fact_mult` using :func:`fact_additive_finite`.

    ``branch_cut`` may be an angle or a callable ``theta = f(|z|)`` (curved cut).
    """
    Fp, Fm = fact_additive_finite(_log(G, branch_cut), contour, prec)
    return (lambda k: np.exp(Fp(k))), (lambda k: np.exp(Fm(k)))


def _continue(G, Gp, Gm, contour, shift_p, shift_m):
    """Analytic continuation across the factorisation contours.

    ``G+`` is used above ``Cp`` and ``G/G-`` below it; ``G-`` below ``Cm``
    and ``G/G+`` above.  Requires ``Re C(t) = t + const``.

    Unlike MATLAB (which evaluated both branches everywhere and blended with
    0/1 masks) only the needed branch is evaluated.  Besides being faster --
    ``G`` is the expensive kernel -- this avoids ``NaN * 0 = NaN`` leaking
    from the unused branch.
    """
    def plus(k1):
        k1 = np.asarray(k1, dtype=complex)
        above = heaviside(np.imag(k1) - np.imag(contour.Cp(np.real(k1) - shift_p))).astype(bool)
        out = np.empty(k1.shape, dtype=complex)
        if above.any():
            out[above] = Gp(k1[above])
        if (~above).any():
            kb = k1[~above]
            out[~above] = G(kb) / Gm(kb)
        return out

    def minus(k1):
        k1 = np.asarray(k1, dtype=complex)
        above = heaviside(np.imag(k1) - np.imag(contour.Cm(np.real(k1) - shift_m))).astype(bool)
        out = np.empty(k1.shape, dtype=complex)
        if above.any():
            ka = k1[above]
            out[above] = G(ka) / Gp(ka)
        if (~above).any():
            out[~above] = Gm(k1[~above])
        return out

    return plus, minus


def fact_mult_cont(G, contour, prec, branch_cut):
    """``Fact_Mult_Cont.m`` -- :func:`fact_mult` + analytic continuation
    (requires ``Re C(t) = t``)."""
    Gp, Gm = fact_mult(G, contour, prec, branch_cut)
    return _continue(G, Gp, Gm, contour, 0.0, 0.0)


def fact_mult_cont_finite(G, contour, prec, branch_cut):
    """``Fact_Mult_Cont_Finite.m`` -- as above using the finite map; allows
    ``Re C(t) = t + Re C(0)``."""
    shift_p = np.real(contour.Cp(0.0))
    shift_m = np.real(contour.Cm(0.0))
    Gp, Gm = fact_mult_finite(G, contour, prec, branch_cut)
    return _continue(G, Gp, Gm, contour, shift_p, shift_m)


def _scaled(G, Gff, contour, prec, branch_cut, cont):
    H = lambda k: G(k) / Gff.f(k)  # -> 1 at infinity
    Hplus, Hminus = cont(H, contour, prec, branch_cut)
    return Struct(f=G, p=lambda k: Hplus(k) * Gff.p(k), m=lambda k: Hminus(k) * Gff.m(k))


def fact_mult_cont_scaled(G, Gff, contour, prec, branch_cut) -> Struct:
    """``Fact_Mult_Cont_Scaled.m`` -- factorise ``G`` whose far-field form
    ``Gff.f = Gff.p * Gff.m`` is known analytically.  Returns ``g`` with
    ``g.f``, ``g.p``, ``g.m``."""
    return _scaled(G, Gff, contour, prec, branch_cut, fact_mult_cont)


def fact_mult_cont_scaled_finite(G, Gff, contour, prec, branch_cut) -> Struct:
    """``Fact_Mult_Cont_Scaled_Finite.m`` -- the version every scattering
    script uses."""
    return _scaled(G, Gff, contour, prec, branch_cut, fact_mult_cont_finite)


#: Samples used by :func:`winding_number` (the kernel is evaluated at each).
WINDING_SAMPLES = 8001


def winding_number(G: Callable, C: Callable, t_max=200.0, n=None) -> int:
    """Winding number of ``G`` round 0 along the contour ``z = C(t)``,
    ``t in [-t_max, t_max]`` (unwrapped change of ``arg G`` over ``2 pi``).

    Multiplicative factorisation of ``G -> 1`` needs this to be zero: otherwise
    ``log G`` cannot be continuous along the contour and the Cauchy split is
    wrong *and* contour-dependent.  The kernels' far-field factors take an
    integer ``W`` to cancel it (see :func:`.kernels.ff_wake_kernel`).
    """
    n = WINDING_SAMPLES if n is None else n
    t = np.sinh(np.linspace(-np.arcsinh(t_max), np.arcsinh(t_max), n))  # dense near 0
    g = G(C(t))
    g = g[np.isfinite(g)]  # the kernel can overflow at very large |k1|
    ph = np.unwrap(np.angle(g))
    return int(np.round((ph[-1] - ph[0]) / (2 * np.pi)))


def choose_tanh_height(G: Callable, Gff0: Struct, heights=(1, 2, 3, 4)):
    """Height of a tanh factorisation contour ``t + i h tanh(t)`` lying above
    every zero of ``G`` in the upper half-plane (within ``max(heights)``).

    Zeros of a Wiener-Hopf kernel that is a ratio of dispersion functions
    (downstream / upstream) are downstream modes, which by Briggs-Bers must lie
    below the factorisation contour; unstable ones rise with frequency and
    shear and can cross a fixed height-1 contour.  Raising the contour changes
    the winding of ``G/Gff`` once per zero crossed, so the lowest height whose
    winding matches the tallest has all those zeros below it.  Returns
    ``(height, winding of G/Gff0 along it)``.
    """
    from .contours import tanh_contour

    wind = {}

    def w(h):
        if h not in wind:
            wind[h] = winding_number(lambda k: G(k) / Gff0.f(k), tanh_contour(0, h, -1, 1)[0])
        return wind[h]

    target = w(heights[-1])
    b = heights[0] if w(heights[0]) == target else next(h for h in heights[1:] if w(h) == target)
    if b != heights[0]:
        log.info("kernel zero(s) above the height-%g contour; using height %g", heights[0], b)
    return b, w(b)


def fact_mult_scaled_auto_winding(G: Callable, Gff_of_W: Callable, contour: Struct, prec: Struct,
                                  branch_cut, W="auto", max_W=3) -> Struct:
    """:func:`fact_mult_cont_scaled_finite` with the far-field winding chosen
    automatically.

    ``Gff_of_W(W)`` builds the analytic far-field factorisation with integer
    winding correction ``W`` (every ``FFKernel_*`` / ``FFWake_Kernel`` takes
    one).  The original scripts always used ``W = 0``.  When ``G/Gff`` winds
    round 0 along the contour, that makes ``log(G/Gff)`` discontinuous and the
    result wrong and contour-dependent.  Here ``W`` is chosen, by measurement
    along ``contour.Cp``, so that ``G/Gff`` has zero winding; an error is raised
    if no ``|W| <= max_W`` achieves this.  Returns the factorised Struct with
    the chosen ``W`` attached.
    """
    if W == "auto":
        wind = lambda w: winding_number(lambda k: G(k) / Gff_of_W(w).f(k), contour.Cp)
        n = wind(0)
        W = 0
        if n:
            for cand in sorted({-n, n, *range(-max_W, max_W + 1)}, key=abs):
                if wind(cand) == 0:
                    W = cand
                    break
            else:
                raise RuntimeError(f"cannot cancel winding number {n} of G/Gff with |W| <= {max_W}")
            log.info("G/Gff winds %d times along the factorisation contour; using W = %d", n, W)
    out = fact_mult_cont_scaled_finite(G, Gff_of_W(W), contour, prec, branch_cut)
    out.W = W
    return out


def add_fact_b(F, C, dC, N1, N2, N3, t1, t2):
    """``AddFactB.m`` (2017, legacy) -- **cannot be ported: dependencies missing.**

    Calls ``FullIntervalIntA`` and ``CauchyKernel``, neither of which is in the
    repository.  The routine is superseded by :func:`fact_additive` (April
    2018), which performs the same Cauchy-integral split; use that instead.
    """
    raise NotImplementedError(
        "AddFactB depends on FullIntervalIntA.m and CauchyKernel.m, which are "
        "missing from the repository. Use fact_additive / fact_additive_finite.")
