"""Known-answer tests for the numerical building blocks."""
import numpy as np
import pytest

from shear_flow_scattering import (complex_funcs as cf, contours, factorisation as fa, frobenius,
                                   inversion, ivp, profiles, quadrature as q, zeros)
from shear_flow_scattering._utils import Struct, pchip
from shear_flow_scattering.pridmore_brown import (comp_dec_aug18a_skip, comp_hard_aug18a_skip,
                                                  dec_aug18a, dec_aug18a_gridx, gamma_ff,
                                                  hard_aug18a, _comp_full, ic_hard)


# --------------------------------------------------------------------------- quadrature
def test_lgwt_matches_matlab_ordering_and_is_exact():
    x, w = q.lgwt(5, 0, 2)
    assert np.all(np.diff(x) < 0)  # MATLAB lgwt runs from b down to a
    assert np.isclose(np.sum(w * x ** 9), 2 ** 10 / 10)


@pytest.mark.parametrize("rule, f, exact", [
    (lambda: q.gauss_legendre(8), lambda x: x ** 4, 2 / 5),
    (lambda: q.gauss_hermite(8), lambda x: x ** 2, np.sqrt(np.pi) / 2),
    (lambda: q.gauss_laguerre(8, 0.0), lambda x: x ** 3, 6.0),
])
def test_eigen_gauss_rules(rule, f, exact):
    x, w = rule()
    assert np.isclose(np.sum(w * f(x)), exact)


def test_adaptive_integral_finite_infinite_complex():
    assert np.isclose(q.adaptive_integral(np.exp, 0, 1), np.e - 1)
    assert np.isclose(q.adaptive_integral(np.exp, -np.inf, 0), 1)
    assert np.isclose(q.adaptive_integral(lambda x: np.exp(-x ** 2), -np.inf, np.inf), np.sqrt(np.pi))
    # residue theorem round a box
    assert np.isclose(q.polyline_integral(lambda z: 1 / z, q.box_contour([-1, 1], [-1, 1])), 2j * np.pi)


# --------------------------------------------------------------------------- complex functions
def test_branch_cut_functions():
    z = np.array([1 + 1j, -2 + 0.5j, 0.3 - 4j])
    assert np.allclose(cf.log_a(z, np.pi), np.log(z))
    for cut in (np.pi / 2, 3 * np.pi / 2):
        assert np.allclose(cf.sqrt_a(z, cut) ** 2, z)
        assert np.allclose(np.exp(cf.log_a(z, cut)), z)
    # cut along +i axis: arg in (-3pi/2, pi/2]
    assert np.isclose(cf.log_a(-1 + 0j, np.pi / 2).imag, -np.pi)
    assert np.allclose(np.exp(cf.log_curved_branch(z, lambda r: 0 * r + 0.3)), z)


def test_pchip_matches_linear_data_and_complex():
    x = np.linspace(0, 1, 7)
    y = (2 + 1j) * x + 1
    xi = np.linspace(0, 1, 23)
    assert np.allclose(pchip(x[::-1], y[::-1], xi), (2 + 1j) * xi + 1)


# --------------------------------------------------------------------------- RK4
def test_rk4_damped_oscillator():
    k = np.array([0.1, 0.5])
    phi, x = ivp.fixed_step_ivp_test(5.0, 2001, k)
    wd = np.sqrt(1 - k ** 2)
    exact = np.exp(-k * x[:, None]) * (np.cos(wd * x[:, None]) + k / wd * np.sin(wd * x[:, None]))
    assert np.allclose(phi[:, 0, :], exact, atol=1e-10)
    end = ivp.fixed_step_ivp_test_gridx_skip(np.array([5.0, 5.0]), 2001, k)
    assert np.allclose(end[0], exact[-1])
    _, _, pe = ivp.fixed_step_ivp_test_gridx(np.array([2.5, 5.0]), 2001, k)
    assert np.isclose(pe[0, 1], exact[-1, 1])


# --------------------------------------------------------------------------- profiles
@pytest.mark.parametrize("U", [
    profiles.u_linear(0.3, 1, 1), profiles.u_parabolic(0.3, 1, 1), profiles.u_cubic(0.2, 1.0),
    profiles.u_cubic_b([0.2, 0.7, 1], [0, 0.4, 1]), profiles.u_quartic([0.2, 0.7, 1], [0, 0.4, 1]),
    profiles.u_parabolic_b([0.2, 0.7, 1], [0, 0.4, 1]),
    profiles.u_pw_parabolic([0.6, 0.7, 1], [0, 0.1, 1]),
    profiles.u_power(0.3, 1, 1, 3), profiles.u_power(0.3, 1, 1, 8),
])
def test_profile_derivatives(U):
    x, h = np.array([0.05, 0.37, 0.8]), 1e-6
    assert np.allclose((U.f(x + h) - U.f(x - h)) / (2 * h), U.df(x), atol=1e-6)
    assert np.allclose((U.df(x + h) - U.df(x - h)) / (2 * h), U.d2f(x), atol=1e-5)


def test_parabolic_profile_edge_values():
    U = profiles.u_parabolic(0.3, 1, 1)
    assert np.isclose(U.f(0), 0.3) and np.isclose(U.f(1), 1) and np.isclose(U.df(1), 0)


# --------------------------------------------------------------------------- Pridmore-Brown
def test_gamma_ff_dispersion_and_factors():
    gamma, gp, gm = gamma_ff(0.3, 1.0, 5.0, 1.0)
    k = np.array([-3.0, 0.1 + 0.2j, 2 - 1j])
    assert np.allclose(gamma(k) ** 2, k ** 2 + 0.3 ** 2 - (1 - k) ** 2 / 25)
    assert np.allclose(gamma(k), gp(k) * gm(k))
    assert np.all(gamma(np.array([-50.0, 50.0])).real > 0)  # decays for large |k1|


def test_uniform_flow_decaying_solution_is_exponential():
    """For U, c constant, Pridmore-Brown reduces to phi'' = gamma^2 phi."""
    U, c, omega, delta = profiles.u_const(0.4), profiles.c_const(3.0), 1.0, 1.0
    k1 = np.array([-2.0, 0.7 + 0.3j, 3.0])
    g = gamma_ff(0.0, omega, 3.0, 0.4)[0](k1)
    dec = dec_aug18a(k1, 0.0, omega, U, c, delta, 400)
    assert np.allclose(dec.phi, np.exp(-np.outer(g, dec.x)), rtol=1e-8)
    for HP in (0, 1, -1):  # analytic in x2 -> path-independent end value
        d = comp_dec_aug18a_skip(k1, 0.0, omega, U, c, delta, 400, contours.path_to(delta, 0.0, HP))
        assert np.allclose(d.phi, 1.0, rtol=1e-7)


def test_path_independence_away_from_critical_layer():
    U, c = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5)
    k1 = np.array([-1.5, -0.3 + 0.2j])  # omega - U k1 never vanishes on [0, 1]
    a = comp_dec_aug18a_skip(k1, 0.0, 1.0, U, c, 1.0, 2000, contours.straight_path(1.0, 0.0))
    b = comp_dec_aug18a_skip(k1, 0.0, 1.0, U, c, 1.0, 2000, contours.contour_semicircle(1.0, 1))
    assert np.allclose(a.v, b.v, rtol=1e-6)


def test_real_and_complex_solvers_agree_and_hard_ic():
    U, c, om = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5), 1.0
    k1 = np.array([-1.0, 0.5 + 0.5j])
    d_real = dec_aug18a(k1, 0.0, om, U, c, 1.0, 801)
    d_comp = comp_dec_aug18a_skip(k1, 0.0, om, U, c, 1.0, 801, contours.straight_path(1.0, 0.0))
    assert np.allclose(d_real.v[:, 0], d_comp.v)
    h = hard_aug18a(k1, 0.0, om, U, c, 1.0, 801)
    assert np.allclose(h.p[:, 0], 1) and np.allclose(h.v[:, 0], 0, atol=1e-12)
    hs = comp_hard_aug18a_skip(k1, 0.0, om, U, c, 801, contours.straight_path(0.0, 1.0))
    assert np.allclose(h.p[:, -1], hs.p)


def test_gridx_decaying_solution_pairs_x_with_phi():
    U, c = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5)
    k1 = np.array([-1.0, -0.5])
    x2_min = np.array([0.2, 0.5])
    d = dec_aug18a_gridx(x2_min, k1, 0.0, 1.0, U, c, 1.0, 301)
    assert np.allclose(d.x[:, 0], x2_min) and np.allclose(d.x[:, -1], 1.0)
    # value at x2_min must equal the skip-solver value there
    for j in range(2):
        s = comp_dec_aug18a_skip(k1[j:j + 1], 0.0, 1.0, U, c, 1.0, 301,
                                 contours.straight_path(1.0, x2_min[j]))
        assert np.isclose(d.phi[j, 0], s.phi[0])


# --------------------------------------------------------------------------- factorisation
def _straight():
    return Struct(Cp=lambda t: t + 0j, dCp=lambda t: np.ones_like(t) + 0j,
                  Cm=lambda t: t + 0j, dCm=lambda t: np.ones_like(t) + 0j)


def test_additive_factorisation_known_split():
    F = lambda k: 1 / (k ** 2 + 1)
    Fp, Fm = fa.fact_additive_finite(F, _straight(), Struct(N=4000))
    assert np.isclose(Fp(0.5j), -1 / (2j * (0.5j + 1j)), atol=1e-6)
    assert np.isclose(Fm(-0.5j), 1 / (2j * (-0.5j - 1j)), atol=1e-6)


def test_multiplicative_factorisation_and_continuation():
    G = lambda k: (k ** 2 + 4) / (k ** 2 + 1)
    Gp, Gm = fa.fact_mult_cont_finite(G, contours.default_factorisation_contour(), Struct(N=4000), np.pi)
    for k in (0.5j, 3 + 0.2j, -0.5j):  # last one uses the analytic continuation G/G-
        assert np.isclose(Gp(k), (k + 2j) / (k + 1j), rtol=1e-5)
        assert np.isclose(Gm(k), (k - 2j) / (k - 1j), rtol=1e-5)


# --------------------------------------------------------------------------- zeros
def test_newton_and_argument_principle():
    f = lambda z: (z - 0.3) * (z + 0.2j)
    df = lambda z: 2 * z - 0.3 + 0.2j
    assert np.isclose(zeros.newton_raphson(f, df, 0.5 + 0.1j), 0.3)
    assert np.isclose(zeros.zeros_poles(f, df, [-1, 1], [-1, 1]), 2)
    zz = np.sort_complex(zeros.simple_zero_finder(f, df, [-1, 1], [-1, 1]))
    assert np.allclose(zz, np.sort_complex([0.3, -0.2j]))
    many = lambda z: np.sin(3 * z)
    found = np.sort(zeros.complete_zero_finder(many, lambda z: 3 * np.cos(3 * z), [-3.3, 3.2], [-1, 0.9]).real)  # avoid zeros on a bisection line
    assert np.allclose(found, np.pi / 3 * np.arange(-3, 4), atol=1e-6)


def test_zero_finder_recursive_finds_lower_half_plane_zeros():
    roots = np.array([1.0, 2 + 1j, 1 - 1j])
    f = lambda z: np.prod([z - r for r in roots], axis=0)
    out = zeros.zero_finder_recursive(f, [1.1, 1.9 + 0.9j, 1.1 - 0.9j])
    key = lambda a: sorted(a, key=lambda z: (round(z.real, 6), round(z.imag, 6)))
    assert np.allclose(key(out.z), key(roots), atol=1e-8)


def test_zeros_box_counts():
    f = lambda z: (z - 0.5) * (z - 0.2j) * (z + 0.3 + 0.4j)
    Z = zeros.zeros_box(f, [-1, 1], [-1, 1])
    assert (Z.z, Z.r, Z.U, Z.L) == (3, 1, 1, 1)


# --------------------------------------------------------------------------- Frobenius
def test_frobenius_bessel_j0():
    y = frobenius.two_solutions_complex(lambda z: 1 + 0 * z, lambda z: z ** 2)
    assert y.I == 2  # repeated indicial root
    assert np.isclose(y.y1(0.5 + 0j), 0.9384698072408129)
    for fn, dfn in ((y.y1, y.dy1), (y.y2, y.dy2)):  # check ODE residual
        z, h = 0.5 + 0.1j, 1e-4
        d2 = (dfn(z + h) - dfn(z - h)) / (2 * h)
        assert abs(z ** 2 * d2 + z * dfn(z) + z ** 2 * fn(z)) < 1e-6


def test_frobenius_integer_separation():
    # z^2 y'' - 2 y + z^2 y = 0 : indicial roots 2 and -1
    y = frobenius.two_solutions_complex(lambda z: 0 * z, lambda z: -2 + z ** 2)
    assert y.I == 3 and np.isclose(y.s1, 2) and np.isclose(y.s2, -1)
    for fn, dfn in ((y.y1, y.dy1), (y.y2, y.dy2)):
        z, h = 0.4 + 0.2j, 1e-4
        d2 = (dfn(z + h) - dfn(z - h)) / (2 * h)
        assert abs(z ** 2 * d2 + (-2 + z ** 2) * fn(z)) < 1e-5


# --------------------------------------------------------------------------- inversion
def test_integral_nov2018_matches_loop():
    rng = np.random.default_rng(0)
    wk, K = rng.normal(size=4) + 0j, rng.normal(size=4) + 1j
    X, Y = rng.normal(size=(3, 4)), rng.normal(size=(5, 4))
    loop = sum(wk[j] * K[j] * np.outer(X[:, j], Y[:, j]) for j in range(4))
    assert np.allclose(inversion.integral_nov2018(wk, K, X, Y), loop)


def test_cl_ellipse_encloses_cut():
    Cp, Cm = contours.finite_cl_ellipse(400, 1.0, 3.0, 0.2, 0.3)
    loop = np.sum(Cp.wk / (Cp.k - 2.0)) + np.sum(Cm.wk / (Cm.k - 2.0))
    assert np.isclose(abs(loop), 2 * np.pi, rtol=1e-6)


def test_semicircle_paths():
    for a, b, I in ((1.0, 0.0, 1), (0.0, 0.4, 1), (1.0, 0.3, -1)):
        z = contours.semicircle_path(a, b, I)
        assert np.isclose(z.f(0.0), a) and np.isclose(z.f(1.0), b)
        assert np.sign(z.f(0.5).imag) == I
        s, h = 0.3, 1e-6
        assert np.isclose((z.f(s + h) - z.f(s - h)) / (2 * h), z.df(s))


def test_vectorised_y_matches_row_by_row():
    U, c = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5)
    x2 = np.array([0.0, 0.3, 0.7, 1.2])
    k = np.array([-1.0, 0.4 + 0.4j, 2.0 - 0.3j])
    Y = inversion.y_continuous_scattering(x2, k, 0.0, 1.0, U, c, 1.0, 300, HP=1)
    for i in range(3):
        d = comp_dec_aug18a_skip(k, 0.0, 1.0, U, c, 1.0, 300, contours.semicircle_path(1.0, x2[i], 1))
        assert np.allclose(Y.p[i], d.p)
    g = gamma_ff(0.0, 1.0, 5.0, 1.0)[0](k)
    assert np.allclose(Y.phi[3], np.exp(-1.2 * g))


def test_greens_function_continuous_at_source():
    U, c = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5)
    k = np.array([-1.0, 0.3 + 0.5j])
    y2 = 0.4
    Y = inversion.y_greens_hard(np.array([y2 - 1e-9, y2 + 1e-9]), y2, k, 0.0, 1.0, U, c, 1.0, 800, 0)
    assert np.allclose(Y.phi[0], Y.phi[1], rtol=1e-6)


@pytest.mark.parametrize("U0", [0.3, 0.0])
def test_hardsoft_kernel_factorisation(U0):
    """K+ K- = K between two separated tanh contours (Contour_Setup.m's
    'pair of separated Tanh contours').  U0 = 0 exercises the no-slip far
    field, which returned NaN before the complex-log fix."""
    from shear_flow_scattering import kernels
    U, c, om = profiles.u_parabolic(U0, 1, 1), profiles.c_const(5), 1.0
    zk = contours.contour_semicircle(1.0, 1)
    Kern = lambda k: kernels.kernel_hardsoft(k, 0.0, om, U, c, 1.0, 400, zk)
    Kff = (kernels.ff_kernel_hardsoft_slip(om, U.f(0), c.f(0), 0.0, 0) if U0 else
           kernels.ff_kernel_hardsoft_noslip(om, U.df(0), c.f(0), 0.0, 0))
    Cp, dCp = contours.tanh_contour(-0.15, 1, -1, 1)
    Cm, dCm = contours.tanh_contour(0.15, 1, -1, 1)
    g = fa.fact_mult_cont_scaled_finite(Kern, Kff, Struct(Cp=Cp, dCp=dCp, Cm=Cm, dCm=dCm),
                                        Struct(N=1000), np.pi)
    k = contours.tanh_contour(0, 1, -1, 1)[0](np.array([-4, -1.5, -0.3, 0.2, 0.8, 2.0, 4.0]))
    assert np.allclose(g.p(k) * g.m(k), Kern(k), rtol=1e-3)


def test_distribution_weight_jacobian_ratio():
    """'jacobian' differs from 'thesis' by exactly (dy2/dkappa1)^2."""
    from shear_flow_scattering.distributed import distribution_weight
    U, c, om = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5), 2.0
    y = np.array([0.0, 0.3, 0.8])
    dy_dk = -U.f(y) ** 2 / (om * U.df(y))
    assert np.allclose(distribution_weight(om, U, c, y, "jacobian"),
                       distribution_weight(om, U, c, y, "thesis") * dy_dk ** 2)
    with pytest.raises(ValueError):
        distribution_weight(om, U, c, y, "other")


def test_jacobian_weighting_equals_y2_integral():
    """With 'jacobian', sum_j g w_j omega_c = int_0^delta c^2 omega_c / U' dy2 (the
    literal source integral, up to the sign convention of dkappa1)."""
    from shear_flow_scattering.distributed import distribution_weight, kappa1_quad
    U, c, om = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5), 1.0
    k, w, y = kappa1_quad(om, U, 1.0, 400)
    oc = lambda y: (1 - y) ** 2
    lhs = np.sum(distribution_weight(om, U, c, y, "jacobian") * oc(y) * w)
    yy, wy = q.lgwt(400, 0, 1)
    rhs = np.sum(25 * oc(yy) / U.df(yy) * wy)
    assert np.isclose(abs(lhs), abs(rhs), rtol=1e-4)


# --------------------------------------------------------------------------- regularised Green's function
def _hard_ops(om, U):
    C0 = lambda k: 1j * (om - U.f(0) * k)
    dC0 = lambda k: -1j * U.df(0) * k
    return (lambda k: -C0(k) ** 2), (lambda k: -3 * C0(k) * dC0(k)), C0


def test_convected_greens_matches_causal_limit():
    """The regular form equals lim_{eta -> 0+} of the unregularised formula at
    kappa1 + i eta (finite there; straight paths are legitimate)."""
    from shear_flow_scattering.inversion import y_greens_convected
    from shear_flow_scattering.pridmore_brown import comp_generic_aug18a_skip as gen
    U, c, om, d = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5), 1.0, 1.0
    l1, l2, C0 = _hard_ops(om, U)
    C = lambda x, k: 1j * (om - U.f(x) * k)
    y2, x2 = 0.4, np.array([0.0, 0.15, 0.6, 0.9])
    kap = om / U.f(y2)
    ref = y_greens_convected(x2, y2, kap, 0.0, om, U, c, l1, l2, d, 2000).p
    errs = []
    for eta in (1e-2, 1e-3):
        k = np.atleast_1d(kap + 1j * eta)
        N = 20000
        s0 = comp_dec_aug18a_skip(k, 0.0, om, U, c, d, N, contours.straight_path(d, 0.0))
        D = (l1(k) * s0.dphi + l2(k) * s0.phi)[0]
        phl = lambda x: gen(k, 0.0, om, U, c, l1, l2, N, contours.straight_path(0.0, x)).phi[0]
        phd = lambda x: comp_dec_aug18a_skip(k, 0.0, om, U, c, d, N, contours.straight_path(d, x)).phi[0]
        g = [-C(y2, k)[0] ** 3 * (phl(y2) * phd(xx) if xx > y2 else phd(y2) * phl(xx)) for xx in x2]
        pg = np.array([-C(xx, k)[0] ** 3 * gi for xx, gi in zip(x2, g)]) / (c.f(0) ** 2 * C0(k)[0] ** 4 * D)
        errs.append(np.max(np.abs(pg - ref)) / np.max(np.abs(ref)))
    assert errs[1] < 0.02 and errs[1] < errs[0] / 5  # converging linearly in eta


@pytest.mark.parametrize("y2", [0.1, 0.4])
def test_convected_greens_continuous_and_resolution_independent(y2):
    from shear_flow_scattering.inversion import y_greens_convected
    U, c, om = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5), 1.0
    l1, l2, _ = _hard_ops(om, U)
    x = np.array([y2 - 0.01, y2 + 0.01])
    a = y_greens_convected(x, y2, om / U.f(y2), 0.0, om, U, c, l1, l2, 1.0, 1000)
    b = y_greens_convected(x, y2, om / U.f(y2), 0.0, om, U, c, l1, l2, 1.0, 4000)
    assert np.isclose(a.p[0], a.p[1], rtol=1e-3)   # original jumps ~((1-y2)/y2)^3
    assert np.isclose(a.p0, b.p0, rtol=1e-6)        # original grows like N^3
    # the anti-causal branch is the complex conjugate (real-coefficient problem)
    m = y_greens_convected(x, y2, om / U.f(y2), 0.0, om, U, c, l1, l2, 1.0, 1000, HP=-1)
    assert np.isclose(m.D, np.conj(a.D), rtol=1e-8)
    with pytest.raises(ValueError):
        y_greens_convected(x, y2, om / U.f(y2), 0.0, om, U, c, l1, l2, 1.0, 1000, HP=0)


# --------------------------------------------------------------------------- factorisation winding
def test_wake_factorisation_winding_and_contour_independence():
    """At omega = 5 the wake kernel winds once round 0 along the factorisation
    contour; with the automatic W correction K+ / K- no longer depend on the
    contour (with the original W = 0 they differed by factors of ~2)."""
    from shear_flow_scattering.scripts.trailing_edge import _wake_factorisation
    U, c, om = profiles.u_parabolic(0.3, 1, 1), profiles.c_const(5.0), 5.0
    kp, km = np.array([1 + 4j]), np.array([2 - 4j])
    vals = []
    for fc in ((0, 1, -1, 1), (0, 2, -0.5, 0.5)):
        C, dC = contours.tanh_contour(*fc)
        K = _wake_factorisation(0.0, om, U, c, 1.0, U, c, 1.0, np.inf, 300,
                                Struct(Cp=C, dCp=dC, Cm=C, dCm=dC), 1500, 3 * np.pi / 2)
        assert K.W == 1
        vals.append((K.p(kp)[0], K.m(km)[0]))
    assert np.allclose(vals[0], vals[1], rtol=1e-3)


# --------------------------------------------------------------------------- non-analytic / fitted profiles
def test_fit_analytic_matches_afshari_profiles():
    for T in profiles.afshari_profiles():
        F = profiles.fit_analytic(T, degree=8)
        assert F.analytic and F.monotonic and F.inflection_free
        assert F.fit_max < 4e-3
        assert np.isclose(F.f(0.0), T.f(0.0)) and np.isclose(F.f(1.0), 1.0) and abs(F.df(1.0)) < 1e-12


def test_indented_path_matches_semicircle_for_analytic_profile():
    """A parabola flagged non-analytic is integrated on the indented real-axis
    path; it must reproduce the (valid) semicircle result, including on the
    critical-layer cut."""
    P = profiles.u_parabolic(0.3, 1, 1)
    Pna = profiles.Profile(P.f, P.df, P.d2f, analytic=False)
    c, om = profiles.c_const(5.0), 5.0
    k = np.array([6.0 + 0j, 9.0, 12.0 + 0.01j, 7.0 - 0.02j, 3 + 1j])
    for HP in (1, -1):
        ref = comp_dec_aug18a_skip(k, 0.0, om, P, c, 1.0, 3000, contours.contour_semicircle(1.0, HP)).v
        ind = comp_dec_aug18a_skip(k, 0.0, om, Pna, c, 1.0, 3000, contours.contour_semicircle(1.0, HP)).v
        assert np.allclose(ind, ref, rtol=1e-4)


def test_two_piece_on_cut_matches_approach_from_above():
    """For the piecewise profile, the indented on-cut value is the limit of the
    straight-path value at kappa1 + i eta (eta -> 0+)."""
    sheared = profiles.afshari_profiles()[1]
    c, om = profiles.c_const(35.0), 5.0
    kap = om / np.array([0.6, 0.85])
    on = comp_dec_aug18a_skip(kap, 0.0, om, sheared, c, 1.0, 6000, contours.contour_semicircle(1.0, 1)).v
    ab = comp_dec_aug18a_skip(kap + 0.01j, 0.0, om, sheared, c, 1.0, 20000, contours.straight_path(1.0, 0.0)).v
    assert np.allclose(on, ab, rtol=0.06)   # O(eta) offset at eta = 0.01


def test_spurious_critical_points_and_safe_ellipse():
    """A degree-8 fit has spurious roots of U(z) = omega/k1 inside the radius-1/2
    semicircle; the path is flattened to avoid them, which restores agreement
    with the straight path (the semicircle itself gives O(1) errors)."""
    from shear_flow_scattering.pridmore_brown import spurious_critical_points
    om = 5.0
    F = profiles.fit_analytic(profiles.afshari_profiles()[1])
    P = profiles.u_parabolic(0.3, 1, 1)
    k = om / np.array([0.6, 0.85]) + 0j
    assert np.all(spurious_critical_points(P, om, k, 1.0, 0.0, 1, 0.5) == 0)
    assert np.all(spurious_critical_points(F, om, k, 1.0, 0.0, 1, 0.5) == 1)
    kk = np.array([5.5 + 1j, 9.0 + 1j])            # critical points off the axis: straight path is exact
    c = profiles.c_const(5.0)
    ref = comp_dec_aug18a_skip(kk, 0.0, om, F, c, 1.0, 20000, contours.straight_path(1.0, 0.0)).v
    v = comp_dec_aug18a_skip(kk, 0.0, om, F, c, 1.0, 3000, contours.contour_semicircle(1.0, 1)).v
    assert np.allclose(v, ref, rtol=1e-5)


def test_tikz_export(tmp_path):
    """matlab2tikz-style export: lines (log axis, NaN gaps), pcolormesh raster
    with colorbar, polar axes."""
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    from shear_flow_scattering.tikz import save_figure_tikz

    fig = plt.figure()
    a = fig.add_subplot(1, 3, 1)
    f = np.logspace(2, 4, 5)
    y = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    a.plot(f, y, "C0--", label=r"$y_2 = 0.1\delta$")
    a.plot(f, y + 1, "o", color="C1")
    a.set_xscale("log"); a.set_xlabel("$f$"); a.legend()
    b = fig.add_subplot(1, 3, 2)
    x1, x2 = np.linspace(-1, 2, 7), np.linspace(0, 1, 4)
    data = np.add.outer(x2, x1)                       # increases up and right
    m = b.pcolormesh(x1, x2, data, shading="auto", cmap="viridis")
    fig.colorbar(m, ax=b)
    c = fig.add_subplot(1, 3, 3, projection="polar")
    c.plot(np.linspace(0, 2 * np.pi, 50), np.ones(50), label="unit"); c.legend()

    paths = save_figure_tikz(fig, str(tmp_path / "fig"))
    assert [p[-5:] for p in paths] == ["a.tex", "b.tex", "c.tex"]
    ta, tb, tc = (open(p).read() for p in paths)
    assert "xmode=log" in ta and "unbounded coords=jump" in ta and "dashed" in ta
    assert r"\addlegendentry{$y_2 = 0.1\delta$}" in ta and ta.count(r"\addlegendentry") == 1
    assert "nan" in ta and "mark=*" in ta
    assert "graphics[" in tb and "colorbar" in tb and "fig_b-1.png" in tb
    img = imread(str(tmp_path / "fig_b-1.png"))
    assert img.shape[:2] == (4, 7)
    # top-right pixel is the largest value -> the yellow end of viridis
    assert img[0, -1, 0] > 0.9 and img[-1, 0, 0] < 0.3
    assert "axis equal image" in tc and r"\addlegendentry{unit}" in tc
    for t in (ta, tb, tc):
        assert t.count(r"\begin{axis}") == t.count(r"\end{axis}") == 1
        assert t.count("{") == t.count("}")
    plt.close(fig)
