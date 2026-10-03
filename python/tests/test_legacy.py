"""Basic functionality of the unused legacy ports.

These functions were ported from the MATLAB for completeness but nothing in
the package, the example scripts or the paper figures calls them (each
carries an "Unused legacy port" warning in its docstring).  The tests check
that each does what its docstring says, with known answers where possible.
"""
import numpy as np
import pytest

from shear_flow_scattering import complex_funcs as cf, contours, profiles, quadrature as q, zeros
from shear_flow_scattering._utils import Struct

PARABOLA = profiles.u_parabolic(0.3, 1, 1)
UNIFORM = profiles.u_const(0.4)
C5 = profiles.c_const(5.0)
C3 = profiles.c_const(3.0)


def _fd(f, t, h=1e-6):
    return (f(t + h) - f(t - h)) / (2 * h)


def _a(*t):
    """Evaluation points as an array (the heaviside-based contours need arrays)."""
    return np.array(t, dtype=float)


def _hard_ops(omega, U):
    C0 = lambda k: 1j * (omega - U.f(0) * k)
    dC0 = lambda k: -1j * U.df(0) * k
    return (lambda k: -C0(k) ** 2), (lambda k: -3 * C0(k) * dC0(k))


# --------------------------------------------------------------------------- _utils, complex_funcs
class TestLegacyUtilities:
    def test_progress_yields_every_item_in_order(self):
        from shear_flow_scattering._utils import progress
        assert list(progress(range(7), "x", every=3)) == list(range(7))

    def test_power_b_with_constant_exponent(self):
        z = np.array([1 + 1j, -2 + 0.5j, 0.3 - 4j])
        assert np.allclose(cf.power_b(z, lambda z: 2 + 0 * z, np.pi), z ** 2)
        assert np.allclose(cf.power_b(z, lambda z: 0.5 + 0 * z, np.pi), np.sqrt(z))

    def test_nan_to_zero(self):
        out = cf.nan_to_zero(lambda x: np.array([1.0, np.nan, 2.0 + 1j]), None)
        assert np.allclose(out, [1, 0, 2 + 1j])

    def test_log_cts_continues_past_the_principal_branch(self):
        C = lambda t: np.exp(1j * t)               # C(0) = 1, winds round 0
        dC = lambda t: 1j * np.exp(1j * t)
        assert np.isclose(cf.log_cts(3 * np.pi, C, dC), 3j * np.pi, atol=1e-8)

    def test_log_cts_f_from_minus_infinity(self):
        a, b = 1j, 2j                               # F -> 1 at -inf along the real line
        F = lambda z: (z - a) / (z - b)
        dF = lambda z: (a - b) / (z - b) ** 2
        val = cf.log_cts_f(0.3, F, dF, lambda e: e, lambda e: 1.0 + 0 * e)
        assert np.isclose(val, np.log(F(0.3)), atol=1e-7)


# --------------------------------------------------------------------------- contours
class TestLegacyContours:
    def test_tanh_contour_b_limits_and_derivative(self):
        a, bm, bp = 0.2, -0.7, 1.3
        C, dC = contours.tanh_contour_b(a, bm, bp, 1.5, 0.8)
        assert np.allclose(C(_a(1e-9, -1e-9)), 1j * a)          # continuous through 0 (see below for t = 0)
        assert np.allclose(C(_a(60.0, -60.0)).imag, [bp, bm])
        t = _a(-1.3, -0.4, 0.5, 2.0)                # both halves (sign bug in the original)
        assert np.allclose(dC(t), _fd(C, t))

    @pytest.mark.xfail(strict=True, reason="KNOWN FAILURE (inherited from TanhContourB.m): at t = 0 both "
                       "halves are switched on, since heaviside(0) = heaviside(-0) = 1, so C(0) = 2ia instead of ia "
                       "and dC(0) is the sum of both slopes.")
    def test_tanh_contour_b_at_zero(self):
        C, _ = contours.tanh_contour_b(0.2, -0.7, 1.3, 1.5, 0.8)
        assert np.isclose(C(_a(0.0))[0], 0.2j)

    def test_gauss_contour(self):
        Cinf, H, W = -0.5, 1.2, 3.0
        C, dC = contours.gauss_contour(Cinf, H, W)
        assert np.isclose(C(0.0).imag, 0, atol=1e-12) and np.isclose(C(W).imag, 0, atol=1e-12)
        assert np.isclose(C(W / 2).imag, H)
        assert np.isclose(C(50.0).imag, Cinf) and np.isclose(C(-50.0).imag, Cinf)
        for t in (-1.0, 0.7, 2.1):
            assert np.isclose(dC(t), _fd(C, t))

    def test_ellipse_contour(self):
        H, W = 0.8, 2.0
        C, dC = contours.ellipse_contour(H, W)
        assert np.isclose(C(_a(W / 2))[0], W / 2 + 1j * H)
        assert np.allclose(C(_a(-1.0, 3.0)), [-1.0, 3.0])
        t = _a(0.3, 1.4)                             # the original lacked the factor 1i
        assert np.allclose(dC(t), _fd(C, t))

    def test_veitch_contour_b_heights(self):
        Dp, Wp, Dm, Wm = 0.6, 1.0, 0.9, 2.0
        C, dC = contours.veitch_contour_b(Dp, Wp, Dm, Wm)
        assert np.isclose(C(_a(0.0))[0], 0)
        assert np.allclose(C(_a(Wp, -Wm)).imag, [Dp, Dm])                  # peaks above the axis
        assert np.all(np.abs(C(_a(200.0, -200.0)).imag) < 1e-5)
        t = _a(-1.5, 0.6)
        assert np.allclose(dC(t), _fd(C, t))

    def test_parabola_through_three_points(self):
        zm, z0, zp = -2 - 1j, 0.5j, 3 - 0.5j
        P = contours.parabola_3points(zm, z0, zp)
        assert np.allclose([P.f(-1.0), P.f(0.0), P.f(1.0)], [zm, z0, zp])
        assert np.isclose(P.df(0.3), _fd(P.f, 0.3)) and np.isclose(P.d2f(0.3), _fd(P.df, 0.3))

    def test_finite_parabola_integrates_decaying_function(self):
        # infinite parabola through -1-i, 0.5i, 1-i: the pole at 2i lies above it
        C = contours.finite_parabola_3points(4000, -1 - 1j, 0.5j, 1 - 1j)
        val = np.sum(C.wk / (C.k - 2j) ** 2)          # decays like 1/z^2, no pole crossed
        assert abs(val) < 1e-6

    def test_cl_ellipse_contour_encloses_the_cut(self):
        kinf, k0, a, b = 1.0, 3.0, 0.2, 0.3
        C = contours.cl_ellipse_contour(kinf, k0, a, b)
        t, w = q.gauss_legendre(400)
        assert np.isclose(C.U(-1.0).imag, 0, atol=1e-12) and np.isclose(C.U(1.0).imag, 0, atol=1e-12)
        assert np.allclose(C.L(t), np.conj(C.U(t)))
        assert np.isclose(C.dU(0.3), _fd(C.U, 0.3))
        loop = np.sum(w * C.dU(t) / (C.U(t) - 2.0)) - np.sum(w * C.dL(t) / (C.L(t) - 2.0))
        assert np.isclose(abs(loop), 2 * np.pi, rtol=2e-3)   # dU ~ 1/sqrt at the ends limits Gauss accuracy
        with pytest.raises(ValueError):
            contours.cl_ellipse_contour(3.0, 1.0, a, b)


# --------------------------------------------------------------------------- distributed
class TestLegacyDistributed:
    def test_integrand_distributed_formula(self):
        from shear_flow_scattering.distributed import integrand_distributed
        val = integrand_distributed(2.0, PARABOLA, C5, 0.7 + 0.1j, 1.5, 0.3 - 0.2j, 0.4)
        assert np.isclose(val, -2.0 * 25 * (0.7 + 0.1j) * 1.5 / ((0.3 - 0.2j) * PARABOLA.f(0.4) ** 2))

    def test_wall_pressure_is_identically_zero(self):
        """As documented: Q0 contains (omega - U(y2) kappa1)^4 = 0."""
        from shear_flow_scattering.distributed import wall_pressure
        l1, l2 = _hard_ops(1.0, PARABOLA)
        out = wall_pressure(np.array([0.3, 0.6]), 1.0, 0.0, PARABOLA, C5, 1.0, l1, l2, 300, 1)
        assert np.all((out == 0) | np.isnan(out))

    def test_scrappy_i0_path_independent_off_the_critical_layer(self):
        from shear_flow_scattering.distributed import scrappy_i0
        k = np.array([-1.0 + 0.3j, 0.5 + 0.8j])
        a = scrappy_i0(k, 0.4, 0, N_bl=800)
        b = scrappy_i0(k, 0.4, 1, N_bl=800)
        assert np.all(np.isfinite(a)) and np.allclose(a, b, rtol=1e-6)


# --------------------------------------------------------------------------- factorisation
class TestLegacyFactorisation:
    def test_prec_setup_defaults(self):
        from shear_flow_scattering.factorisation import prec_setup
        p = prec_setup()
        assert (p.N, p.max, p.min) == (1000, 20, -20)

    def test_fact_mult_cont_scaled_known_split(self):
        """G = Gff (k^2+4)/(k^2+1) with Gff = (k+3i)(k-3i): G+ = Gff+ (k+2i)/(k+i)."""
        from shear_flow_scattering import factorisation as fa
        Gff = Struct(p=lambda k: k + 3j, m=lambda k: k - 3j)
        Gff.f = lambda k: Gff.p(k) * Gff.m(k)
        G = lambda k: Gff.f(k) * (k ** 2 + 4) / (k ** 2 + 1)
        g = fa.fact_mult_cont_scaled(G, Gff, contours.default_factorisation_contour(), fa.prec_setup(N=4000), np.pi)
        kp = np.array([0.5j, 3 + 0.2j, -2 + 1j])                 # above the contour
        km = np.conj(kp)
        Gp = lambda k: (k + 3j) * (k + 2j) / (k + 1j)
        Gm = lambda k: (k - 3j) * (k - 2j) / (k - 1j)
        # constant-free checks (a split is only defined up to K+ -> A K+, K- -> K-/A) ...
        assert np.allclose(g.p(kp) * g.m(kp), G(kp), rtol=1e-6)
        assert np.allclose(g.p(kp) / g.p(kp[0]), Gp(kp) / Gp(kp[0]), rtol=2e-4)
        assert np.allclose(g.m(km) / g.m(km[0]), Gm(km) / Gm(km[0]), rtol=2e-4)
        # ... and the scaling by Gff pins the constant (G/Gff -> 1 at infinity)
        assert np.allclose(g.p(kp), Gp(kp), rtol=2e-4)

    def test_add_fact_b_not_portable(self):
        from shear_flow_scattering.factorisation import add_fact_b
        with pytest.raises(NotImplementedError):
            add_fact_b(None, None, None, 1, 1, 1, 0, 1)


# --------------------------------------------------------------------------- inversion, pridmore_brown
class TestLegacyInversionWrappers:
    def test_integrate_fields_applies_inversion_sum_to_each_field(self):
        from shear_flow_scattering import FIELDS, inversion
        rng = np.random.default_rng(1)
        wk, K, X = rng.normal(size=4) + 0j, rng.normal(size=4) + 1j, rng.normal(size=(3, 4))
        Y = Struct(**{f: rng.normal(size=(5, 4)) + 0j for f in FIELDS})
        out = inversion.integrate_fields(wk, K, X, Y)
        for f in FIELDS:
            assert np.allclose(out[f], inversion.integral_nov2018(wk, K, X, Y[f]))

    @pytest.mark.parametrize("name, truncate", [("y_continuous_scattering_cl_nov2018", False),
                                                ("y_continuous_scattering_cl_nov2018_mod", True)])
    def test_cl_wrappers(self, name, truncate):
        from shear_flow_scattering import inversion
        x2, k = np.array([0.2, 0.8]), np.array([1.5 + 0.2j, 2.0 - 0.3j])
        a = getattr(inversion, name)(x2, k, 0.0, 1.0, PARABOLA, C5, 1.0, 200, 1)
        b = inversion.y_continuous_scattering(x2, k, 0.0, 1.0, PARABOLA, C5, 1.0, 200, 1, truncate)
        assert np.allclose(a.p, b.p) and np.allclose(a.phi, b.phi)

    def test_dirty_dispersion_uniform_flow(self):
        """Uniform flow: phi_d = exp(-gamma x2), so D_p = p(0) = -C(0)^3."""
        from shear_flow_scattering.inversion import dirty_dispersion
        k = np.array([-2.0, 0.7 + 0.3j, 3.0])
        assert np.allclose(dirty_dispersion(k, 0.0, 1.0, UNIFORM, C3, 1.0, 400), -(1j * (1.0 - 0.4 * k)) ** 3)


class TestLegacyGridSolvers:
    """Per-wavenumber end points, against the uniform-flow exact solutions
    phi_d = exp(-gamma x2) (decaying) and p = cosh(gamma x2) (hard wall)."""
    k = np.array([-2.0, 0.7 + 0.3j, 3.0])

    def _gamma(self):
        from shear_flow_scattering.pridmore_brown import gamma_ff
        return gamma_ff(0.0, 1.0, 3.0, 0.4)[0](self.k)

    def test_hard_aug18a_gridx(self):
        from shear_flow_scattering.pridmore_brown import hard_aug18a_gridx
        x2max = np.array([0.3, 0.6, 0.9])
        h = hard_aug18a_gridx(x2max, self.k, 0.0, 1.0, UNIFORM, C3, 400)
        assert np.allclose(h.p_end, np.cosh(self._gamma() * x2max), rtol=1e-8)

    def test_hard_aug18a_gridx_skip(self):
        from shear_flow_scattering.pridmore_brown import hard_aug18a_gridx_skip
        x2max = np.array([0.3, 0.6, 0.9])
        h = hard_aug18a_gridx_skip(x2max, self.k, 0.0, 1.0, UNIFORM, C3, 1.0, 400)
        assert np.allclose(h.p, np.cosh(self._gamma() * x2max), rtol=1e-8)

    def test_dec_aug18a_gridx_skip(self):
        from shear_flow_scattering.pridmore_brown import dec_aug18a_gridx_skip
        x2min = np.array([0.1, 0.4, 0.8])
        d = dec_aug18a_gridx_skip(x2min, self.k, 0.0, 1.0, UNIFORM, C3, 1.0, 400)
        assert np.allclose(d.phi, np.exp(-self._gamma() * x2min), rtol=1e-8)


# --------------------------------------------------------------------------- steepest descent
class TestLegacySteepestDescent:
    def test_contour_passes_through_saddle_with_constant_phase(self):
        from shear_flow_scattering.pridmore_brown import gamma_ff
        from shear_flow_scattering.steepest_descent import far_field_saddle, steepest_descent_contour
        omega, U, c = 5.0, 1.0, 5.0
        x1, x2 = np.array([-3.0, 2.0]), np.array([4.0])
        sd = steepest_descent_contour(x1, x2, omega, U, c, 0.0, 41)     # odd: t = 0 is a node
        gamma = gamma_ff(0.0, omega, c, U)[0]
        for i in range(x1.size):
            K = sd.SD[0, i]
            th = np.arctan2(x2[0], x1[i])
            assert np.isclose(K[20], far_field_saddle([th], omega, U, c).saddle[0])
            E = -1j * K * x1[i] - gamma(K) * x2[0]                       # exponent of the integrand
            assert np.allclose(E.imag, E[20].imag, atol=1e-8)           # steepest descent: constant phase
            assert np.all(E.real <= E[20].real + 1e-9)                   # maximum at the saddle

    def test_anonymous_contour_matches_discretised(self):
        from shear_flow_scattering.steepest_descent import steepest_descent_contour, steepest_descent_contour_anon
        omega, U, c = 5.0, 1.0, 5.0
        sd = steepest_descent_contour(np.array([-3.0]), np.array([4.0]), omega, U, c, 0.0, 21)
        _, _, K1_SD, k1_SD = steepest_descent_contour_anon(omega, U, c, 0.0)
        assert np.allclose(k1_SD(sd.t, sd.Theta[0, 0]), sd.SD[0, 0])
        assert np.isclose(K1_SD(0.0, 0.7), sd.K0 * np.cos(0.7))

    def test_oct_18_far_field_finite_and_consistent(self):
        from shear_flow_scattering.steepest_descent import oct_18_ff_sd
        th = np.linspace(0.3, np.pi - 0.3, 7)
        r = oct_18_ff_sd(th, 0.5, 0.0, 1.0, PARABOLA, C5, 1.0, N_BL=600)
        assert np.all(np.isfinite(r.p))
        assert np.allclose(r.p, r.Cauchy * r.hard_y2.p / r.Wronskian
                           * (1.0 - PARABOLA.f(0.5) * np.real_if_close(r.Cauchy * 0 + 1)) ** 0
                           * ((1.0 - PARABOLA.f(0.5) * (1 / (2j * np.pi * r.Cauchy) + 1.0 / PARABOLA.f(0.5)))
                              / (1.0 - PARABOLA.f(0) * (1 / (2j * np.pi * r.Cauchy) + 1.0 / PARABOLA.f(0.5)))) ** 4)

    def test_sep18_unfinished_integration_runs(self):
        """Unfinished in the original (see docstring); only check it runs and
        returns fields of the right shape."""
        from shear_flow_scattering.steepest_descent import sep18_integrate_sd
        x1, x2 = np.array([-2.0, 1.0, 3.0]), np.array([0.3, 1.5])
        I, *_ = sep18_integrate_sd(x1, x2, 0.5, 0.0, 1.0, PARABOLA, C5, 1.0, N_SD=21, N_BL=200)
        assert I.p.shape == (x2.size, x1.size)


# --------------------------------------------------------------------------- zeros
class TestLegacyZeroTools:
    f = staticmethod(lambda z: (z - 0.3) * (z + 0.2j))
    df = staticmethod(lambda z: 2 * z - 0.3 + 0.2j)

    def test_complex_integral_line_and_box(self):
        assert np.isclose(zeros.complex_integral_line(lambda z: z ** 2, 0, 1 + 1j, 20), (1 + 1j) ** 3 / 3)
        assert np.isclose(zeros.complex_integral_box(lambda z: 1 / z, [-1, 1], [-1, 1], 200, 200), 2j * np.pi,
                          rtol=1e-6)

    def test_number_zeros_box(self):
        assert zeros.number_zeros_box(Struct(f=self.f, df=self.df), [-1, 1], [-1, 1], 200, 200, 1e-3) == 2

    def test_number_zeros_bisection_locates_zeros(self):
        r = [0.33 + 0.11j, -0.12 - 0.21j]                 # off the 2^3 bisection grid lines
        f = lambda z: (z - r[0]) * (z - r[1])
        df = lambda z: 2 * z - r[0] - r[1]
        N, Z, sz, X, Y = zeros.number_zeros_bisection(Struct(f=f, df=df), [-1, 1], [-1, 1], 100, 100, 3, 1e-3)
        assert N.sum() == 2
        assert np.allclose(np.sort_complex(sz[N == 1]), np.sort_complex(r), atol=1e-8)

    def test_simple_zero_tracker_follows_a_root(self):
        p = np.linspace(1, 4, 7)
        z = zeros.simple_zero_tracker(lambda z, p: z ** 2 - p, p, 1.0 + 0j)
        assert np.allclose(z, np.sqrt(p))

    def test_pointwise_tools(self):
        x = np.linspace(0, 1, 2001)
        assert np.isclose(zeros.pointwise_trapezium(lambda t: t ** 2, x), 1 / 3, atol=1e-6)
        s = np.linspace(-1, 1, 2001)
        assert np.isclose(zeros.pointwise_trapezium_box(lambda z: 1 / z, s, s), 2j * np.pi, rtol=1e-5)
        assert np.isclose(zeros.pointwise_zp(self.f, s, s), 2, atol=1e-3)

    def test_pv_trapz(self):
        assert abs(zeros.pv_trapz(lambda x: 1 / x, -1, 1, [0.0], 1e-3, N=20001)) < 1e-6
        # exact for any symmetric excision: log(eps) + (log 2 - log eps)
        assert np.isclose(zeros.pv_trapz(lambda x: 1 / x, -1, 2, [0.0], 1e-2, N=200001), np.log(2), atol=1e-5)

    def test_branch_cut_integral_detects_a_cut(self):
        B = np.linspace(-2, -0.5, 200) + 0j                  # along the cut of the principal sqrt
        d = 1e-3
        I1, I2 = zeros.branch_cut_integral(lambda z: z ** 2, B, d)
        assert abs(I1 + I2) < 20 * d                           # analytic: only O(delta) end effects
        J1, J2 = zeros.branch_cut_integral(np.sqrt, B, 1e-3)
        assert abs(J1 + J2) > 0.5                              # sqrt jumps across the negative real axis

    def test_zeros_poles_branch_cuts_not_portable(self):
        with pytest.raises(NotImplementedError):
            zeros.zeros_poles_branch_cuts(None, None, [0, 1], [0, 1])


# --------------------------------------------------------------------------- trailing edge
@pytest.mark.slow
def test_wall_pressure_far_field_routine_loops_the_single_frequency_version():
    from shear_flow_scattering.scripts import trailing_edge as te
    r = te.wall_pressure_point_source_far_field_routine(y2=[0.3], omega=[2.0], N_bl=200)
    single = te.wall_pressure_point_source_far_field(2.0, te._study_profiles()[0], np.array([0.3]), N_bl=200)
    assert r.FF_Pressure_base.shape == (1, 1)
    assert np.allclose(r.FF_Pressure_base[0], single)


class TestHeaviside:
    """``_utils.heaviside`` -- used by the legacy contours above and by
    ``factorisation`` (always with arrays there)."""

    def test_arrays(self):
        from shear_flow_scattering._utils import heaviside
        assert np.array_equal(heaviside(np.array([-1.0, 0.0, 2.0 + 1j])), [0.0, 1.0, 1.0])

    @pytest.mark.xfail(strict=True, reason="KNOWN FAILURE: heaviside(x) = (np.real(x) >= 0).astype(float) "
                       "fails for a scalar x, where the comparison returns a Python bool (no .astype); so "
                       "tanh_contour_b, ellipse_contour and veitch_contour_b crash when called with a scalar t.")
    def test_scalar(self):
        from shear_flow_scattering._utils import heaviside
        assert heaviside(0.5) == 1.0 and heaviside(-0.5) == 0.0
