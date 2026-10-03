"""Generic numerical building blocks: quadrature, branch cuts, interpolation, RK4.

Modules: ``quadrature``, ``complex_funcs``, ``_utils.pchip``, ``ivp``.
"""
import numpy as np
import pytest

from shear_flow_scattering import complex_funcs as cf, ivp, quadrature as q
from shear_flow_scattering._utils import pchip


class TestGaussRules:
    def test_lgwt_matches_matlab_ordering_and_is_exact(self):
        x, w = q.lgwt(5, 0, 2)
        assert np.all(np.diff(x) < 0)  # MATLAB lgwt runs from b down to a
        assert np.isclose(np.sum(w * x ** 9), 2 ** 10 / 10)

    @pytest.mark.parametrize("rule, f, exact", [
        (lambda: q.gauss_legendre(8), lambda x: x ** 4, 2 / 5),
        (lambda: q.gauss_hermite(8), lambda x: x ** 2, np.sqrt(np.pi) / 2),
        (lambda: q.gauss_laguerre(8, 0.0), lambda x: x ** 3, 6.0),
    ])
    def test_eigen_gauss_rules(self, rule, f, exact):
        x, w = rule()
        assert np.isclose(np.sum(w * f(x)), exact)


class TestAdaptiveIntegral:
    def test_finite_and_infinite_ranges(self):
        assert np.isclose(q.adaptive_integral(np.exp, 0, 1), np.e - 1)
        assert np.isclose(q.adaptive_integral(np.exp, -np.inf, 0), 1)
        assert np.isclose(q.adaptive_integral(lambda x: np.exp(-x ** 2), -np.inf, np.inf), np.sqrt(np.pi))

    def test_polyline_contour_residue(self):
        assert np.isclose(q.polyline_integral(lambda z: 1 / z, q.box_contour([-1, 1], [-1, 1])), 2j * np.pi)


class TestBranchCuts:
    z = np.array([1 + 1j, -2 + 0.5j, 0.3 - 4j])

    def test_principal_log(self):
        assert np.allclose(cf.log_a(self.z, np.pi), np.log(self.z))

    @pytest.mark.parametrize("cut", [np.pi / 2, 3 * np.pi / 2])
    def test_rotated_cuts_invert(self, cut):
        assert np.allclose(cf.sqrt_a(self.z, cut) ** 2, self.z)
        assert np.allclose(np.exp(cf.log_a(self.z, cut)), self.z)

    def test_cut_along_positive_imaginary_axis(self):
        # arg in (-3 pi/2, pi/2]
        assert np.isclose(cf.log_a(-1 + 0j, np.pi / 2).imag, -np.pi)

    def test_curved_branch_log_inverts(self):
        assert np.allclose(np.exp(cf.log_curved_branch(self.z, lambda r: 0 * r + 0.3)), self.z)


class TestPchip:
    def test_reproduces_linear_complex_data_from_descending_nodes(self):
        x = np.linspace(0, 1, 7)
        y = (2 + 1j) * x + 1
        xi = np.linspace(0, 1, 23)
        assert np.allclose(pchip(x[::-1], y[::-1], xi), (2 + 1j) * xi + 1)


class TestFixedStepRK4:
    k = np.array([0.1, 0.5])

    def _exact(self, x):
        wd = np.sqrt(1 - self.k ** 2)
        return np.exp(-self.k * x[:, None]) * (np.cos(wd * x[:, None]) + self.k / wd * np.sin(wd * x[:, None]))

    def test_damped_oscillator_full_solution(self):
        phi, x = ivp.fixed_step_ivp_test(5.0, 2001, self.k)
        assert np.allclose(phi[:, 0, :], self._exact(x), atol=1e-10)

    def test_end_point_and_per_k_grids(self):
        exact_end = self._exact(np.array([5.0]))[0]
        end = ivp.fixed_step_ivp_test_gridx_skip(np.array([5.0, 5.0]), 2001, self.k)
        assert np.allclose(end[0], exact_end)
        _, _, pe = ivp.fixed_step_ivp_test_gridx(np.array([2.5, 5.0]), 2001, self.k)
        assert np.isclose(pe[0, 1], exact_end[1])
