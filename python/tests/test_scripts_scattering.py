"""Surface-junction scattering scripts (``scripts.scattering``): hard->soft,
hard->soft->hard and impedance junctions, near and far field, end to end at low
resolution."""
import numpy as np
import pytest

from shear_flow_scattering.scripts import scattering

from _helpers import assert_finite

X1 = np.linspace(-2, 2, 9)
X2 = np.linspace(0, 1.5, 7)
SMALL = dict(N_bl=150, N_C=40)


@pytest.mark.slow
class TestHardSoft:
    @pytest.mark.parametrize("incident", ["regularised", "original"])
    def test_near_and_far_field(self, incident):
        r = scattering.hard_soft_continuous(x1=X1, x2=X2, prec_N=300, incident=incident, **SMALL)
        assert_finite(r.I_full, (X1.size, X2.size))
        r = scattering.far_field_hard_soft(N_theta=15, N_bl=150, prec_N=300, incident=incident)
        assert_finite(r.FF, (15,))


@pytest.mark.slow
class TestHardSoftHard:
    def test_near_and_far_field(self):
        x1 = np.linspace(-7, 3, 11)
        r = scattering.hard_soft_hard(x1=x1, x2=X2, N_CL=60, prec_N=300, **SMALL)
        assert_finite(r.I_full, (x1.size, X2.size))
        r = scattering.far_field_hard_soft_hard(N_theta=15, N_bl=150, prec_N=300)
        assert_finite(r.FF, (15,))


@pytest.mark.slow
class TestImpedanceJunction:
    def test_near_and_far_field(self):
        # x2 grid avoids y2 = 0.5 exactly: there phi, dphi, v of the incident field are
        # undefined (NaN by design; p is the finite limit)
        r = scattering.imp_imp_continuous(x1=X1, x2=np.linspace(0, 1.5, 8), N_CL=60, prec_N=300, **SMALL)
        assert_finite(r.I_full, (X1.size, 8))
        r = scattering.far_field_imp_imp(R=2.0, N_theta=15, N_bl=150, prec_N=300)
        assert_finite(r.FF, (15,))
