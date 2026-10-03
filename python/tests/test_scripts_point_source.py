"""Point-source driver scripts (``scripts.point_source``), end to end at low
resolution: they run and return finite fields of the right shape."""
import numpy as np
import pytest

from shear_flow_scattering import FIELDS
from shear_flow_scattering.scripts import point_source

from _helpers import assert_finite

X1 = np.linspace(-2, 2, 9)
X2 = np.linspace(0, 1.5, 7)
SMALL = dict(N_bl=150, N_C=40)


class TestPointSourceNearField:
    @pytest.mark.slow
    @pytest.mark.parametrize("bc", ["hard", "pressure", "impedance"])
    def test_generic_boundary_conditions(self, bc):
        r = point_source.point_source_generic(bc=bc, x1=X1, x2=X2, N_CL=60, **SMALL)
        assert_finite(r.I_tot, (X1.size, X2.size))

    def test_hard(self):
        r = point_source.point_source_hard(x1=X1, x2=X2, N_CL=60, **SMALL)
        assert_finite(r.I_tot, (X1.size, X2.size))

    @pytest.mark.slow
    def test_wake(self):
        x2 = np.linspace(-1.5, 1.5, 9)
        r = point_source.point_source_wake(x1=X1, x2=x2, N_CL=60, **SMALL)
        assert_finite(r.I_tot, (X1.size, x2.size))
        assert np.all(np.isfinite(r.dispersion_zeros))


class TestPointSourceFarField:
    def test_generic(self):
        r = point_source.far_field_generic_point_source(N_theta=15, N_bl=150)
        assert_finite(r.FF, (15,))

    def test_wake(self):
        r = point_source.far_field_wake_point_source(N_theta_p=8, N_theta_m=8, N_bl=150)
        # theta = 0, pi, 2pi put the saddle on the branch point where D = 0 for this
        # (uniform-flow) default: 0 * inf = NaN there, exactly as in MATLAB.
        interior = np.abs(np.sin(r.theta)) > 1e-9
        assert_finite({f: r.FF[f][interior] for f in FIELDS})


class TestPressureZeroTracker:
    def test_zeros_finite_and_in_upper_half_plane(self):
        c0, z = point_source.pressure_zero_tracker(c0=[6, 7], N_bl=200)
        assert np.all(np.isfinite(z)) and np.all(z.imag > 0)
