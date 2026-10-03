"""Wiener-Hopf factorisation (``factorisation``, ``kernels``): known additive and
multiplicative splits, the hard->soft kernel, and the wake kernel's winding."""
import numpy as np
import pytest

from shear_flow_scattering import contours, factorisation as fa, kernels, profiles
from shear_flow_scattering._utils import Struct


def _straight():
    return Struct(Cp=lambda t: t + 0j, dCp=lambda t: np.ones_like(t) + 0j,
                  Cm=lambda t: t + 0j, dCm=lambda t: np.ones_like(t) + 0j)


class TestAdditiveFactorisation:
    @pytest.mark.slow
    def test_known_split(self):
        F = lambda k: 1 / (k ** 2 + 1)
        Fp, Fm = fa.fact_additive_finite(F, _straight(), Struct(N=4000))
        assert np.isclose(Fp(0.5j), -1 / (2j * (0.5j + 1j)), atol=1e-6)
        assert np.isclose(Fm(-0.5j), 1 / (2j * (-0.5j - 1j)), atol=1e-6)


class TestMultiplicativeFactorisation:
    @pytest.mark.slow
    def test_known_split_and_analytic_continuation(self):
        G = lambda k: (k ** 2 + 4) / (k ** 2 + 1)
        Gp, Gm = fa.fact_mult_cont_finite(G, contours.default_factorisation_contour(), Struct(N=4000), np.pi)
        for k in (0.5j, 3 + 0.2j, -0.5j):  # last one uses the analytic continuation G/G-
            assert np.isclose(Gp(k), (k + 2j) / (k + 1j), rtol=1e-5)
            assert np.isclose(Gm(k), (k - 2j) / (k - 1j), rtol=1e-5)


class TestHardSoftKernel:
    @pytest.mark.parametrize("U0", [0.3, 0.0], ids=["slip", "no-slip"])
    def test_product_of_factors_recovers_kernel(self, U0):
        """K+ K- = K between two separated tanh contours.  U0 = 0 exercises the
        no-slip far field, which returned NaN before the complex-log fix."""
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


class TestWakeKernelWinding:
    @pytest.mark.slow
    def test_winding_cancelled_and_contour_independent(self):
        """At omega = 5 the wake kernel winds once round 0 along the factorisation
        contour; with the automatic W correction K+ and K- no longer depend on the
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


class TestFarFieldKernelFactors:
    """Closed-form far-field factors Kff used to scale the numerical
    factorisation: Kff.p * Kff.m = Kff.f, with the winding factor placed on
    the correct branch points (k3 != 0 included)."""

    @pytest.mark.parametrize("k3", [0.0, 0.4])
    @pytest.mark.parametrize("W", [0, 1, -1])
    def test_hardsoft_slip_product_and_winding_factor(self, k3, W):
        from shear_flow_scattering.pridmore_brown import branch_points, gamma_ff
        omega, U0, c0 = 5.0, 0.3, 5.0
        Kff = kernels.ff_kernel_hardsoft_slip(omega, U0, c0, k3, W)
        k = np.array([-3 + 0.5j, -0.2 + 2j, 0.4 - 1j, 3 - 0.5j])
        assert np.allclose(Kff.p(k) * Kff.m(k), Kff.f(k))
        gamma = gamma_ff(k3, omega, c0, U0)[0]
        kb1, kb2 = branch_points(omega, U0, c0, k3)
        C = 1j * (omega - U0 * k)
        assert np.allclose(Kff.f(k) * gamma(k) / -C, ((k - kb2) / (k - kb1)) ** W)

    @pytest.mark.parametrize("Z", [np.inf, 0.0, 2 - 2j], ids=["hard", "soft", "lined"])
    @pytest.mark.parametrize("U0", [0.3, 0.0], ids=["slip", "no-slip"])
    def test_wake_factors_multiply(self, Z, U0):
        omega, c0 = 5.0, 5.0
        Kff = kernels.ff_wake_kernel(0.0, omega, U0, c0, 1.4, -1.4, Z, 0)
        k = np.array([-3 + 0.5j, -0.2 + 2j, 0.4 - 1j, 3 - 0.5j])
        assert np.all(np.isfinite(Kff.f(k)))
        assert np.allclose(Kff.p(k) * Kff.m(k), Kff.f(k))
