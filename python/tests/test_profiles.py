"""Mean-flow profiles (``profiles``): derivatives, edge values, analytic fits."""
import numpy as np
import pytest

from shear_flow_scattering import profiles


class TestProfileDerivatives:
    @pytest.mark.parametrize("U", [
        profiles.u_linear(0.3, 1, 1), profiles.u_parabolic(0.3, 1, 1), profiles.u_cubic(0.2, 1.0),
        profiles.u_cubic_b([0.2, 0.7, 1], [0, 0.4, 1]), profiles.u_quartic([0.2, 0.7, 1], [0, 0.4, 1]),
        profiles.u_parabolic_b([0.2, 0.7, 1], [0, 0.4, 1]),
        profiles.u_pw_parabolic([0.6, 0.7, 1], [0, 0.1, 1]),
        profiles.u_power(0.3, 1, 1, 3), profiles.u_power(0.3, 1, 1, 8),
    ])
    def test_first_and_second_derivatives_match_finite_differences(self, U):
        x, h = np.array([0.05, 0.37, 0.8]), 1e-6
        assert np.allclose((U.f(x + h) - U.f(x - h)) / (2 * h), U.df(x), atol=1e-6)
        assert np.allclose((U.df(x + h) - U.df(x - h)) / (2 * h), U.d2f(x), atol=1e-5)


class TestParabolicProfile:
    def test_edge_values(self):
        U = profiles.u_parabolic(0.3, 1, 1)
        assert np.isclose(U.f(0), 0.3) and np.isclose(U.f(1), 1) and np.isclose(U.df(1), 0)


class TestFitAnalytic:
    @pytest.mark.parametrize("which", [0, 1], ids=["base", "sheared"])
    def test_fits_to_afshari_profiles(self, which):
        T = profiles.afshari_profiles()[which]
        F = profiles.fit_analytic(T, degree=8)
        assert F.analytic and F.monotonic and F.inflection_free
        assert F.fit_max < 4e-3
        assert np.isclose(F.f(0.0), T.f(0.0)) and np.isclose(F.f(1.0), 1.0) and abs(F.df(1.0)) < 1e-12
