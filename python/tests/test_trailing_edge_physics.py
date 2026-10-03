"""Physics checks on the trailing-edge solution (``scripts.trailing_edge``).

These test *what* the solver computes rather than that it runs:

* the near field satisfies the boundary conditions it was built from (hard
  plate, continuity across the wake) and the Kutta condition at the edge;
* in uniform flow the wake kernel and the far field match closed forms;
* the two far-field representations agree, the far field matches direct
  near-field inversion for asymmetric flow and below the plate, and the
  critical angle of an asymmetric wake sits where its formula says.

All are slow (each needs a factorisation of the wake kernel).
"""
import numpy as np
import pytest

from shear_flow_scattering.profiles import c_const, u_parabolic
from shear_flow_scattering.pridmore_brown import gamma_ff
from shear_flow_scattering.scripts import trailing_edge
from shear_flow_scattering.steepest_descent import far_field_saddle

pytestmark = pytest.mark.slow

C5 = c_const(5.0)
OMEGA = 5.0
KAPPA1 = OMEGA / 0.8                 # unit wall-pressure gust convecting at 0.8 U_inf


def physical_amplitude(FF, theta, Uinf, c0, lower=False):
    """lim r^(1/2) |p'| from the far-field script output (as sec42_studies)."""
    sd = far_field_saddle(theta, OMEGA, Uinf, c0, lower=lower)
    return np.abs(np.ravel(FF)) * np.sqrt(sd.K0 / (2 * np.pi)) \
        * (np.cos(theta) ** 2 + sd.beta ** 2 * np.sin(theta) ** 2) ** -0.25


# --------------------------------------------------------------------------- boundary conditions
@pytest.fixture(scope="module", params=[(1.5, 1.0), (1.0, 5.0)], ids=["asym-w1", "sym-w5"])
def edge_field(request):
    """Total field just above (x2 = 0) and just below (x2 = -1e-9) the plate/wake
    line, for a vortex sheet at y2 = 0.5 over a hard plate."""
    Uinf_m, om = request.param
    x1 = np.array([-2, -1, -0.4, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.5, 1, 2.0])
    r = trailing_edge.trailing_edge_lined(
        Up=u_parabolic(0.3, 1, 1), Um=u_parabolic(0.3, Uinf_m, 1), cp=C5, cm=C5, omega=om, Z=np.inf,
        y2=0.5, x1=x1, x2=np.array([-1e-9, 0.0]), N_bl=600, N_C=500, N_CL=250, prec_N=2500)
    T = r.I_total
    # scale: max |p| upstream (finite even if the downstream field is not)
    return x1, T, np.max(np.abs(T.p[x1 < 0]))


class TestTrailingEdgeBoundaryConditions:
    TOL = 1e-3     # relative to max |p| on the upstream line

    def test_hard_plate_no_normal_velocity_above(self, edge_field):
        x1, T, s = edge_field
        up = x1 < 0
        assert np.max(np.abs(T.v[up, 1])) / s < self.TOL

    def test_hard_plate_no_normal_velocity_below(self, edge_field):
        x1, T, s = edge_field
        up = x1 < 0
        assert np.max(np.abs(T.v[up, 0])) / s < self.TOL

    def test_downstream_field_finite(self, edge_field):
        x1, T, _ = edge_field
        dn = x1 > 0
        assert np.all(np.isfinite(T.p[dn])) and np.all(np.isfinite(T.v[dn]))

    def test_wake_pressure_continuous(self, edge_field):
        x1, T, s = edge_field
        dn = x1 > 0
        assert np.max(np.abs(T.p[dn, 1] - T.p[dn, 0])) / s < self.TOL

    def test_wake_normal_velocity_continuous(self, edge_field):
        x1, T, s = edge_field
        dn = x1 > 0
        assert np.max(np.abs(T.v[dn, 1] - T.v[dn, 0])) / s < self.TOL

    def test_kutta_pressure_jump_vanishes_at_edge(self, edge_field):
        """With the Kutta condition the pressure jump across the plate behaves
        like |x1|^(1/2) as x1 -> 0- (without it, like |x1|^(-1/2))."""
        x1, T, s = edge_field
        sel = (x1 < 0) & (x1 >= -0.4)
        dp = np.abs(T.p[sel, 1] - T.p[sel, 0])
        slope = np.polyfit(np.log(-x1[sel]), np.log(dp), 1)[0]
        assert 0.3 < slope < 0.7, f"log-log slope of |dp| near the edge: {slope:.2f}"


# --------------------------------------------------------------------------- uniform flow
UNIFORM = u_parabolic(1.0, 1, 1)     # U = 1 everywhere: no boundary-layer shear


@pytest.fixture(scope="module")
def uniform_far_field():
    th1 = np.radians(np.linspace(5, 175, 18))
    th2 = 2 * np.pi - th1[::-1]
    return th1, th2, trailing_edge.far_field_trailing_edge_lined(
        Up=UNIFORM, Um=UNIFORM, cp=C5, cm=C5, omega=OMEGA, Z=np.inf, theta1=th1, theta2=th2,
        pressure_upstream=1.0, kappa1=KAPPA1, N_bl=600, prec_N=2500)


class TestUniformFlowTrailingEdge:
    """With U = U_inf on both sides the hard-plate wake kernel is
    K = Dw/(D1 D2) = -2C/gamma exactly (C = i(omega - U k1)), so K- is 1/gamma-
    and K+ is C/gamma+, each up to a constant."""

    def test_kernel_factors_match_gamma_factors(self, uniform_far_field):
        K = uniform_far_field[2].K
        _, gp, gm = gamma_ff(0.0, OMEGA, 5.0, 1.0)
        below = np.array([-3 - 3j, -0.5 - 4j, 2 - 3.5j, 5 - 5j])     # below every contour
        above = np.array([-2 + 4j, 1 + 5j, 4 + 4.5j])
        rm = K.m(below) * gm(below)
        rp = K.p(above) * gp(above) / (1j * (OMEGA - above))
        assert np.allclose(rm, rm[0], rtol=1e-3), rm / rm[0]
        assert np.allclose(rp, rp[0], rtol=1e-3), rp / rp[0]

    def test_far_field_matches_closed_form_directivity(self, uniform_far_field):
        """|FF| proportional to |sin Theta| / |gamma-(k_s) (k_s - kappa1)| above the plate
        (sin(theta/2)-like cardioid at low Mach)."""
        th1, _, r = uniform_far_field
        sd = far_field_saddle(th1, OMEGA, 1.0, 5.0)
        gm = gamma_ff(0.0, OMEGA, 5.0, 1.0)[2]
        closed = np.abs(np.sin(sd.Theta)) / np.abs(gm(sd.saddle + 0j) * (sd.saddle - KAPPA1))
        ratio = np.abs(r.FF.p[:th1.size]) / closed
        assert np.allclose(ratio, ratio[0], rtol=1e-3), ratio / ratio[0]

    def test_far_field_symmetric_about_plate(self, uniform_far_field):
        th1, _, r = uniform_far_field
        upper = np.abs(r.FF.p[:th1.size])
        lower = np.abs(r.FF.p[th1.size:])[::-1]     # theta -> 2 pi - theta
        assert np.allclose(upper, lower, rtol=1e-3)


# --------------------------------------------------------------------------- asymmetric far field
UP, UM2 = u_parabolic(0.3, 1, 1), u_parabolic(0.3, 2.0, 1)


def critical_angle(Uinf_m, Uinf_p=1.0, c0=5.0):
    """Upper observer angle at which the saddle crosses the lower stream's
    branch point omega/(c0 + U_inf^-)."""
    M, k0 = Uinf_p / c0, OMEGA / c0
    beta, K0 = np.sqrt(1 - M ** 2), k0 / (1 - M ** 2)
    Th = np.arccos((OMEGA / (c0 + Uinf_m) + M * k0 / (1 - M ** 2)) / K0)
    return np.arctan(np.tan(Th) / beta)


@pytest.fixture(scope="module")
def asymmetric_far_field():
    thc = critical_angle(2.0)
    h = np.radians([1e-3, 3e-3, 1e-2])
    th1 = np.sort(np.concatenate([np.radians(np.linspace(5, 175, 18)), [thc], thc - h, thc + h]))
    th2 = np.radians(np.linspace(185, 355, 18))
    r = trailing_edge.far_field_trailing_edge_lined(
        Up=UP, Um=UM2, cp=C5, cm=C5, omega=OMEGA, Z=np.inf, theta1=th1, theta2=th2,
        pressure_upstream=1.0, kappa1=KAPPA1, N_bl=800, prec_N=3000)
    return th1, th2, thc, h, r


class TestFarFieldConsistency:
    def test_upstream_and_downstream_representations_agree(self, asymmetric_far_field):
        """The '+' (x1 < 0) and '-' (x1 > 0) forms of the scattered field are equal
        when K = K+ K- holds; any factorisation error shows up as a mismatch."""
        _, _, _, _, r = asymmetric_far_field
        for p, m in ((r.FF1p.p, r.FF1m.p), (r.FF2p.p, r.FF2m.p)):
            assert np.allclose(p, m, rtol=1e-4)

    def test_critical_angle_square_root_cusp(self, asymmetric_far_field):
        """A(theta_c) - A(theta_c - h) ~ sqrt(h) on the evanescent side, at the angle
        given by the branch point of the faster lower stream."""
        th1, _, thc, h, r = asymmetric_far_field
        A = physical_amplitude(r.FF.p[:th1.size], th1, 1.0, 5.0)
        Ac = A[np.argmin(np.abs(th1 - thc))]
        below = np.array([A[np.argmin(np.abs(th1 - (thc - hh)))] for hh in h])
        s = np.abs(Ac - below) / np.sqrt(h)
        assert np.allclose(s, s[0], rtol=0.1), s


@pytest.mark.parametrize("theta, lower", [(np.radians(100.0), False), (np.radians(240.0), True)],
                         ids=["above-100deg", "below-240deg"])
def test_far_field_matches_near_field_asymmetric(theta, lower):
    """Saddle-point far field against direct near-field inversion at r = 30 for
    an asymmetric flow (U_inf^- = 1.5), above and below the plate."""
    from shear_flow_scattering.scripts.common import incident_sheet, wall_operators
    Um, y2, r = u_parabolic(0.3, 1.5, 1), 0.5, 30.0
    x1, x2 = r * np.cos(theta), r * np.sin(theta)
    nf = trailing_edge.trailing_edge_lined(Up=UP, Um=Um, cp=C5, cm=C5, omega=OMEGA, Z=np.inf, y2=y2,
                                           x1=np.array([x1]), x2=np.array([x2]),
                                           N_bl=500, N_C=400, N_CL=150, prec_N=1500)
    near = np.sqrt(r) * abs(nf.I_scattered.p[0, 0])
    _, wall = incident_sheet(None, None, y2, OMEGA / UP.f(y2), 0.0, OMEGA, UP, C5, 1.0, 500,
                             wall_operators(OMEGA, UP).hard)
    th1, th2 = ([theta], [np.pi + 0.3]) if not lower else ([0.5], [theta])
    ff = trailing_edge.far_field_trailing_edge_lined(Up=UP, Um=Um, cp=C5, cm=C5, omega=OMEGA, N_bl=500,
                                                     prec_N=1500, theta1=th1, theta2=th2,
                                                     pressure_upstream=wall.p, kappa1=OMEGA / UP.f(y2))
    FF = ff.FF2p.p if lower else ff.FF1p.p
    far = physical_amplitude(FF, np.array([theta]), Um.f(1.0) if lower else 1.0, 5.0, lower=lower)[0]
    assert np.isclose(near, far, rtol=0.03), (near, far)


# --------------------------------------------------------------------------- wake modes
class TestWakeModes:
    """Zeros of D_w found by the sampled argument principle in rectangles that
    avoid the branch cuts (``dw_zeros``)."""

    @staticmethod
    def _zeros(Uinf_m, omega, height=2.0):
        from shear_flow_scattering.scripts.trailing_edge import dw_zeros, plate_operators
        Up, Um = u_parabolic(0.3, 1, 1), u_parabolic(0.3, Uinf_m, 1)
        return dw_zeros(plate_operators(omega, np.inf, Up, Um), omega, 0.0, Up, C5, 1.0, Um, C5, 1.0, height)

    def test_thesis_mode(self):
        """Thesis Fig. 6.2.2 / paper Fig. 3 flow: kappa_w = 0.8720 +- 0.2618i, as quoted."""
        z = self._zeros(1.5, 1.0)
        z = z[np.argsort(z.imag)]                     # conjugate pair: order by imaginary part
        assert z.size == 2
        assert np.allclose(z, [0.87200 - 0.26184j, 0.87200 + 0.26184j], atol=1e-4)

    @pytest.mark.parametrize("Uinf_m, omega, expected", [
        (1.0, 1.0, [1.1833 - 0.2694j, 1.1833 + 0.2694j]),
        (2.0, 1.0, [0.6675 - 0.2559j, 0.6675 + 0.2559j]),
        (1.0, 5.0, []),
    ], ids=["sym-w1", "U2-w1", "sym-w5"])
    def test_other_flows(self, Uinf_m, omega, expected):
        """Including the cases where the original single Newton start returned NaN
        or landed on the end of the critical-layer cut."""
        z = self._zeros(Uinf_m, omega)
        e = np.array(expected, dtype=complex)
        z, e = z[np.argsort(z.imag)], e[np.argsort(e.imag)]
        assert z.size == e.size and np.allclose(z, e, atol=1e-3)

    def test_no_spurious_zero_near_cut_end(self):
        """Uniformly stepped semicircles produce a zero near omega/U0 = 3.33 that
        moves with resolution; the graded paths must not."""
        z = self._zeros(1.5, 1.0)
        assert not np.any(np.abs(z - 1.0 / 0.3) < 0.3)
