"""Trailing-edge scripts (``scripts.trailing_edge``): single-sheet and
distributed forcing, near and far field, wall pressure, and the check that the
saddle-point far field agrees with direct near-field inversion."""
import numpy as np
import pytest

from shear_flow_scattering._utils import Struct
from shear_flow_scattering.scripts import trailing_edge

from _helpers import assert_finite

X1 = np.linspace(-2, 2, 9)
SMALL = dict(N_bl=150, N_C=40)


class TestTrailingEdgeNearField:
    @pytest.mark.slow
    @pytest.mark.parametrize("Z, incident", [(2 - 2j, "regularised"), (np.inf, "regularised"),
                                             (np.inf, "original")], ids=["lined", "hard", "hard-original"])
    def test_runs_finite(self, Z, incident):
        x2 = np.linspace(-1.5, 1.5, 9)
        r = trailing_edge.trailing_edge_lined(Z=Z, x1=X1, x2=x2, N_CL=60, prec_N=300,
                                              omega=1.0, incident=incident, **SMALL)
        assert_finite(r.I_total, (X1.size, x2.size))


class TestTrailingEdgeFarField:
    @pytest.mark.slow
    def test_scripts_run_finite(self):
        r = trailing_edge.far_field_trailing_edge_lined(N_theta_p=8, N_theta_m=8, N_bl=150, prec_N=300)
        assert all(np.all(np.isfinite(r.FF[f])) for f in ("phi", "p", "v"))
        r = trailing_edge.wall_pressure_far_field_trailing_edge_lined(N_theta_p=8, N_theta_m=8,
                                                                     N_bl=150, prec_N=300)
        assert r.FF.p.shape == (16,)
        ff = trailing_edge.wall_pressure_point_source_far_field(
            5.0, trailing_edge._study_profiles()[0], np.linspace(0.1, 0.9, 3), N_bl=150, prec_N=300)
        assert ff.shape == (3,) and np.all(np.isfinite(ff))

    @pytest.mark.slow
    def test_matches_near_field_inversion(self):
        """Saddle-point far field (with the physical normalisation) against direct
        near-field inversion at large r, same forcing (regularised vortex sheet).
        Guards the whole far-field pipeline: factorisation (winding, contour),
        saddle evaluation and amplitude scaling."""
        from shear_flow_scattering.profiles import c_const, u_parabolic
        from shear_flow_scattering.scripts.common import incident_sheet, wall_operators
        from shear_flow_scattering.steepest_descent import far_field_saddle
        U, c, om, y2 = u_parabolic(0.3, 1, 1), c_const(5.0), 5.0, 0.5
        th, r = np.radians(120), 30.0
        x1, x2 = r * np.cos(th), r * np.sin(th)
        nf = trailing_edge.trailing_edge_lined(Up=U, Um=U, cp=c, cm=c, omega=om, Z=np.inf, y2=y2,
                                               x1=np.array([x1]), x2=np.array([-0.5, x2]),
                                               N_bl=500, N_C=400, N_CL=150, prec_N=1500)
        near = np.sqrt(r) * abs(nf.I_scattered.p[0, 1])
        _, wall = incident_sheet(None, None, y2, om / U.f(y2), 0.0, om, U, c, 1.0, 500,
                                 wall_operators(om, U).hard)
        ff = trailing_edge.far_field_trailing_edge_lined(Up=U, Um=U, cp=c, cm=c, omega=om, N_bl=500,
                                                         prec_N=1500, theta1=[th], theta2=[np.pi + 0.3],
                                                         pressure_upstream=wall.p, kappa1=om / U.f(y2))
        sd = far_field_saddle([th], om, 1.0, 5.0)
        far = abs(np.ravel(ff.FF1p.p)[0]) * np.sqrt(sd.K0 / (2 * np.pi)) \
            * (np.cos(th) ** 2 + sd.beta ** 2 * np.sin(th) ** 2) ** -0.25
        assert np.isclose(near, far, rtol=0.03)


@pytest.fixture(scope="module")
def distributed():
    """One low-resolution distributed trailing-edge solve, shared by the tests below."""
    x2 = np.linspace(-1.2, 1.2, 9)
    r = trailing_edge.distributed_trailing_edge(x1=X1, x2=x2, N_distribution=6, N_CL=60, prec_N=300, **SMALL)
    return r, x2


@pytest.mark.slow
class TestDistributedTrailingEdge:
    def test_runs_and_reruns_with_new_vorticity(self, distributed):
        r, x2 = distributed
        assert_finite(r.I_total, (X1.size, x2.size))
        r2 = trailing_edge.rerun_with_vorticity(r.state, r.y2 * (1 - r.y2))
        assert_finite(r2.I_total, (X1.size, x2.size))
        assert not np.allclose(r.I_total.p, r2.I_total.p)

    def test_weighting_toggle(self, distributed):
        r, x2 = distributed
        oc = (1 - r.y2) ** 2
        th = trailing_edge.rerun_with_vorticity(r.state, oc)  # default: as MATLAB
        assert th.weighting == "thesis"
        jac = trailing_edge.rerun_with_vorticity(r.state, oc, weighting="jacobian")
        assert_finite(jac.I_total, (X1.size, x2.size))
        assert not np.allclose(th.I_total.p, jac.I_total.p)
        with pytest.raises(ValueError):
            trailing_edge.rerun_with_vorticity(r.state, oc, weighting="nonsense")


class TestEdgeResolutionCheck:
    def test_quiet_when_outer_nodes_small_warns_when_dominant(self, caplog):
        st = Struct(y2=np.array([0.1, 0.5, 0.9]))
        assert trailing_edge._check_edge_resolution(st, np.array([1.0, 1.0, 0.01]), "jacobian") < 0.1
        assert "under-resolved" not in caplog.text
        assert trailing_edge._check_edge_resolution(st, np.array([1.0, 1.0, 5.0]), "jacobian") > 0.1
        assert "under-resolved" in caplog.text


class TestWallPressure:
    @pytest.mark.slow
    def test_scaled_point_source_script(self):
        r = trailing_edge.wall_pressure_point_source_scaled(y2=[0.2, 0.4], omega=[1.0, 10.0], N_bl=200)
        assert np.all(np.isfinite(r.Full_Pressure_base))
