"""Integration paths and contours (``contours``)."""
import numpy as np
import pytest

from shear_flow_scattering import contours


class TestSemicirclePath:
    @pytest.mark.parametrize("a, b, I", [(1.0, 0.0, 1), (0.0, 0.4, 1), (1.0, 0.3, -1)])
    def test_end_points_side_and_derivative(self, a, b, I):
        z = contours.semicircle_path(a, b, I)
        assert np.isclose(z.f(0.0), a) and np.isclose(z.f(1.0), b)
        assert np.sign(z.f(0.5).imag) == I
        s, h = 0.3, 1e-6
        assert np.isclose((z.f(s + h) - z.f(s - h)) / (2 * h), z.df(s))


class TestCriticalLayerEllipse:
    def test_encloses_the_cut_once(self):
        Cp, Cm = contours.finite_cl_ellipse(400, 1.0, 3.0, 0.2, 0.3)
        loop = np.sum(Cp.wk / (Cp.k - 2.0)) + np.sum(Cm.wk / (Cm.k - 2.0))
        assert np.isclose(abs(loop), 2 * np.pi, rtol=1e-6)
