"""Inversion integrals and Green's functions (``inversion``), including the
regularised Green's function for a convected vortex sheet."""
import numpy as np
import pytest

from shear_flow_scattering import contours, inversion, profiles
from shear_flow_scattering.pridmore_brown import comp_dec_aug18a_skip, gamma_ff

from _helpers import hard_wall_operators

PARABOLA = profiles.u_parabolic(0.3, 1, 1)
C5 = profiles.c_const(5.0)


class TestIntegralNov2018:
    def test_matches_explicit_loop(self):
        rng = np.random.default_rng(0)
        wk, K = rng.normal(size=4) + 0j, rng.normal(size=4) + 1j
        X, Y = rng.normal(size=(3, 4)), rng.normal(size=(5, 4))
        loop = sum(wk[j] * K[j] * np.outer(X[:, j], Y[:, j]) for j in range(4))
        assert np.allclose(inversion.integral_nov2018(wk, K, X, Y), loop)


class TestScatteringY:
    def test_vectorised_matches_row_by_row_and_outer_region(self):
        x2 = np.array([0.0, 0.3, 0.7, 1.2])
        k = np.array([-1.0, 0.4 + 0.4j, 2.0 - 0.3j])
        Y = inversion.y_continuous_scattering(x2, k, 0.0, 1.0, PARABOLA, C5, 1.0, 300, HP=1)
        for i in range(3):
            d = comp_dec_aug18a_skip(k, 0.0, 1.0, PARABOLA, C5, 1.0, 300, contours.semicircle_path(1.0, x2[i], 1))
            assert np.allclose(Y.p[i], d.p)
        g = gamma_ff(0.0, 1.0, 5.0, 1.0)[0](k)
        assert np.allclose(Y.phi[3], np.exp(-1.2 * g))


class TestHardWallGreensFunction:
    def test_continuous_at_source(self):
        k = np.array([-1.0, 0.3 + 0.5j])
        y2 = 0.4
        Y = inversion.y_greens_hard(np.array([y2 - 1e-9, y2 + 1e-9]), y2, k, 0.0, 1.0, PARABOLA, C5, 1.0, 800, 0)
        assert np.allclose(Y.phi[0], Y.phi[1], rtol=1e-6)


class TestConvectedGreensFunction:
    """The regular form p<(y2) phi>(x2) / (c0(0)^2 C(0)^4 D) of the vortex-sheet
    Green's function at kappa1 = omega/U(y2)."""

    @pytest.mark.slow
    def test_matches_causal_limit(self):
        """Equals lim_{eta -> 0+} of the unregularised formula at kappa1 + i eta
        (finite there; straight paths are legitimate)."""
        from shear_flow_scattering.inversion import y_greens_convected
        from shear_flow_scattering.pridmore_brown import comp_generic_aug18a_skip as gen
        om, d = 1.0, 1.0
        l1, l2, C0 = hard_wall_operators(om, PARABOLA)
        C = lambda x, k: 1j * (om - PARABOLA.f(x) * k)
        y2, x2 = 0.4, np.array([0.0, 0.15, 0.6, 0.9])
        kap = om / PARABOLA.f(y2)
        ref = y_greens_convected(x2, y2, kap, 0.0, om, PARABOLA, C5, l1, l2, d, 2000).p
        errs = []
        for eta in (1e-2, 1e-3):
            k = np.atleast_1d(kap + 1j * eta)
            N = 20000
            s0 = comp_dec_aug18a_skip(k, 0.0, om, PARABOLA, C5, d, N, contours.straight_path(d, 0.0))
            D = (l1(k) * s0.dphi + l2(k) * s0.phi)[0]
            phl = lambda x: gen(k, 0.0, om, PARABOLA, C5, l1, l2, N, contours.straight_path(0.0, x)).phi[0]
            phd = lambda x: comp_dec_aug18a_skip(k, 0.0, om, PARABOLA, C5, d, N, contours.straight_path(d, x)).phi[0]
            g = [-C(y2, k)[0] ** 3 * (phl(y2) * phd(xx) if xx > y2 else phd(y2) * phl(xx)) for xx in x2]
            pg = np.array([-C(xx, k)[0] ** 3 * gi for xx, gi in zip(x2, g)]) / (C5.f(0) ** 2 * C0(k)[0] ** 4 * D)
            errs.append(np.max(np.abs(pg - ref)) / np.max(np.abs(ref)))
        assert errs[1] < 0.02 and errs[1] < errs[0] / 5  # converging linearly in eta

    @pytest.mark.slow
    @pytest.mark.parametrize("y2", [0.1, 0.4])
    def test_continuous_resolution_independent_and_causal(self, y2):
        from shear_flow_scattering.inversion import y_greens_convected
        om = 1.0
        l1, l2, _ = hard_wall_operators(om, PARABOLA)
        x = np.array([y2 - 0.01, y2 + 0.01])
        a = y_greens_convected(x, y2, om / PARABOLA.f(y2), 0.0, om, PARABOLA, C5, l1, l2, 1.0, 1000)
        b = y_greens_convected(x, y2, om / PARABOLA.f(y2), 0.0, om, PARABOLA, C5, l1, l2, 1.0, 4000)
        assert np.isclose(a.p[0], a.p[1], rtol=1e-3)   # original jumps ~((1-y2)/y2)^3
        assert np.isclose(a.p0, b.p0, rtol=1e-6)        # original grows like N^3
        # the anti-causal branch is the complex conjugate (real-coefficient problem)
        m = y_greens_convected(x, y2, om / PARABOLA.f(y2), 0.0, om, PARABOLA, C5, l1, l2, 1.0, 1000, HP=-1)
        assert np.isclose(m.D, np.conj(a.D), rtol=1e-8)

    def test_rejects_straight_path_through_source(self):
        from shear_flow_scattering.inversion import y_greens_convected
        l1, l2, _ = hard_wall_operators(1.0, PARABOLA)
        with pytest.raises(ValueError):
            y_greens_convected(np.array([0.3]), 0.4, 1.0 / PARABOLA.f(0.4), 0.0, 1.0, PARABOLA, C5,
                               l1, l2, 1.0, 100, HP=0)
