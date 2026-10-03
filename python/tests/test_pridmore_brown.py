"""Boundary-layer ODE solvers (``pridmore_brown``): far-field gamma, decaying and
hard-wall solutions, path independence, and the path-safety logic for the
critical layer (indented paths, spurious critical points)."""
import numpy as np
import pytest

from shear_flow_scattering import contours, profiles
from shear_flow_scattering.pridmore_brown import (comp_dec_aug18a_skip, comp_hard_aug18a_skip,
                                                  dec_aug18a, dec_aug18a_gridx, gamma_ff,
                                                  hard_aug18a, spurious_critical_points)

PARABOLA = profiles.u_parabolic(0.3, 1, 1)
C5 = profiles.c_const(5.0)


class TestGammaFarField:
    def test_dispersion_relation_and_factors(self):
        gamma, gp, gm = gamma_ff(0.3, 1.0, 5.0, 1.0)
        k = np.array([-3.0, 0.1 + 0.2j, 2 - 1j])
        assert np.allclose(gamma(k) ** 2, k ** 2 + 0.3 ** 2 - (1 - k) ** 2 / 25)
        assert np.allclose(gamma(k), gp(k) * gm(k))
        assert np.all(gamma(np.array([-50.0, 50.0])).real > 0)  # decays for large |k1|


class TestDecayingSolution:
    def test_uniform_flow_is_exponential_on_every_path(self):
        """For U, c constant, Pridmore-Brown reduces to phi'' = gamma^2 phi."""
        U, c, omega, delta = profiles.u_const(0.4), profiles.c_const(3.0), 1.0, 1.0
        k1 = np.array([-2.0, 0.7 + 0.3j, 3.0])
        g = gamma_ff(0.0, omega, 3.0, 0.4)[0](k1)
        dec = dec_aug18a(k1, 0.0, omega, U, c, delta, 400)
        assert np.allclose(dec.phi, np.exp(-np.outer(g, dec.x)), rtol=1e-8)
        for HP in (0, 1, -1):  # analytic in x2 -> path-independent end value
            d = comp_dec_aug18a_skip(k1, 0.0, omega, U, c, delta, 400, contours.path_to(delta, 0.0, HP))
            assert np.allclose(d.phi, 1.0, rtol=1e-7)

    def test_path_independent_away_from_critical_layer(self):
        k1 = np.array([-1.5, -0.3 + 0.2j])  # omega - U k1 never vanishes on [0, 1]
        a = comp_dec_aug18a_skip(k1, 0.0, 1.0, PARABOLA, C5, 1.0, 2000, contours.straight_path(1.0, 0.0))
        b = comp_dec_aug18a_skip(k1, 0.0, 1.0, PARABOLA, C5, 1.0, 2000, contours.contour_semicircle(1.0, 1))
        assert np.allclose(a.v, b.v, rtol=1e-6)

    def test_real_and_complex_solvers_agree(self):
        k1 = np.array([-1.0, 0.5 + 0.5j])
        d_real = dec_aug18a(k1, 0.0, 1.0, PARABOLA, C5, 1.0, 801)
        d_comp = comp_dec_aug18a_skip(k1, 0.0, 1.0, PARABOLA, C5, 1.0, 801, contours.straight_path(1.0, 0.0))
        assert np.allclose(d_real.v[:, 0], d_comp.v)

    def test_gridx_pairs_x_with_phi(self):
        k1 = np.array([-1.0, -0.5])
        x2_min = np.array([0.2, 0.5])
        d = dec_aug18a_gridx(x2_min, k1, 0.0, 1.0, PARABOLA, C5, 1.0, 301)
        assert np.allclose(d.x[:, 0], x2_min) and np.allclose(d.x[:, -1], 1.0)
        for j in range(2):  # value at x2_min equals the skip-solver value there
            s = comp_dec_aug18a_skip(k1[j:j + 1], 0.0, 1.0, PARABOLA, C5, 1.0, 301,
                                     contours.straight_path(1.0, x2_min[j]))
            assert np.isclose(d.phi[j, 0], s.phi[0])


class TestHardWallSolution:
    def test_initial_condition_and_skip_solver(self):
        k1 = np.array([-1.0, 0.5 + 0.5j])
        h = hard_aug18a(k1, 0.0, 1.0, PARABOLA, C5, 1.0, 801)
        assert np.allclose(h.p[:, 0], 1) and np.allclose(h.v[:, 0], 0, atol=1e-12)
        hs = comp_hard_aug18a_skip(k1, 0.0, 1.0, PARABOLA, C5, 801, contours.straight_path(0.0, 1.0))
        assert np.allclose(h.p[:, -1], hs.p)


class TestPathSafety:
    @pytest.mark.slow
    def test_indented_path_matches_semicircle_for_analytic_profile(self):
        """A parabola flagged non-analytic is integrated on the indented real-axis
        path; it must reproduce the (valid) semicircle result, including on the
        critical-layer cut."""
        Pna = profiles.Profile(PARABOLA.f, PARABOLA.df, PARABOLA.d2f, analytic=False)
        k = np.array([6.0 + 0j, 9.0, 12.0 + 0.01j, 7.0 - 0.02j, 3 + 1j])
        for HP in (1, -1):
            ref = comp_dec_aug18a_skip(k, 0.0, 5.0, PARABOLA, C5, 1.0, 3000, contours.contour_semicircle(1.0, HP)).v
            ind = comp_dec_aug18a_skip(k, 0.0, 5.0, Pna, C5, 1.0, 3000, contours.contour_semicircle(1.0, HP)).v
            assert np.allclose(ind, ref, rtol=1e-4)

    @pytest.mark.slow
    def test_two_piece_on_cut_matches_approach_from_above(self):
        """For the piecewise profile, the indented on-cut value is the limit of the
        straight-path value at kappa1 + i eta (eta -> 0+)."""
        sheared = profiles.afshari_profiles()[1]
        c, om = profiles.c_const(35.0), 5.0
        kap = om / np.array([0.6, 0.85])
        on = comp_dec_aug18a_skip(kap, 0.0, om, sheared, c, 1.0, 6000, contours.contour_semicircle(1.0, 1)).v
        ab = comp_dec_aug18a_skip(kap + 0.01j, 0.0, om, sheared, c, 1.0, 20000, contours.straight_path(1.0, 0.0)).v
        assert np.allclose(on, ab, rtol=0.06)   # O(eta) offset at eta = 0.01

    def test_spurious_critical_point_count(self):
        """A degree-8 fit has a spurious root of U(z) = omega/k1 inside the
        radius-1/2 semicircle; a parabola has none."""
        F = profiles.fit_analytic(profiles.afshari_profiles()[1])
        k = 5.0 / np.array([0.6, 0.85]) + 0j
        assert np.all(spurious_critical_points(PARABOLA, 5.0, k, 1.0, 0.0, 1, 0.5) == 0)
        assert np.all(spurious_critical_points(F, 5.0, k, 1.0, 0.0, 1, 0.5) == 1)

    @pytest.mark.slow
    def test_safe_ellipse_restores_straight_path_value(self):
        """The flattened path avoids the spurious root, restoring agreement with the
        straight path (the semicircle itself gives O(1) errors)."""
        F = profiles.fit_analytic(profiles.afshari_profiles()[1])
        kk = np.array([5.5 + 1j, 9.0 + 1j])  # critical points off the axis: straight path is exact
        ref = comp_dec_aug18a_skip(kk, 0.0, 5.0, F, C5, 1.0, 20000, contours.straight_path(1.0, 0.0)).v
        v = comp_dec_aug18a_skip(kk, 0.0, 5.0, F, C5, 1.0, 3000, contours.contour_semicircle(1.0, 1)).v
        assert np.allclose(v, ref, rtol=1e-5)


class TestBranchPoints:
    """Acoustic branch points: the zeros of gamma, including k3 != 0 (where the
    MATLAB scripts had the wrong sign of k3^2; every script uses k3 = 0)."""

    @pytest.mark.parametrize("k3", [0.0, 0.3, 0.9, 1.5], ids=["k3=0", "k3=0.3", "k3=0.9", "k3>k0"])
    @pytest.mark.parametrize("Uinf, cinf, omega", [(1.0, 5.0, 5.0), (0.3, 5.0, 5.0), (1.0, 35.0, 10.0)])
    def test_are_zeros_of_gamma(self, k3, Uinf, cinf, omega):
        from shear_flow_scattering.pridmore_brown import branch_points
        kb1, kb2 = branch_points(omega, Uinf, cinf, k3)
        gamma = gamma_ff(k3, omega, cinf, Uinf)[0]
        scale = omega / cinf + k3
        assert abs(gamma(np.asarray(kb1, complex))) < 1e-6 * scale
        assert abs(gamma(np.asarray(kb2, complex))) < 1e-6 * scale

    @pytest.mark.parametrize("k3", [0.0, 0.3, 0.9])
    def test_ordered_and_each_belongs_to_one_factor(self, k3):
        """kb1 > kb2 when real; gamma+ vanishes at kb1 (its branch point is in the
        lower half-plane under causality), gamma- at kb2."""
        from shear_flow_scattering.pridmore_brown import branch_points
        omega, Uinf, cinf = 5.0, 1.0, 5.0
        kb1, kb2 = branch_points(omega, Uinf, cinf, k3)
        assert np.isreal(kb1) and np.isreal(kb2) and kb1 > kb2
        _, gp, gm = gamma_ff(k3, omega, cinf, Uinf)
        assert abs(gp(np.asarray(kb1, complex))) < 1e-6
        assert abs(gm(np.asarray(kb2, complex))) < 1e-6
