"""Distributed-vorticity weighting (``distributed``): thesis vs Jacobian."""
import numpy as np
import pytest

from shear_flow_scattering import profiles, quadrature as q
from shear_flow_scattering.distributed import distribution_weight, kappa1_quad

PARABOLA = profiles.u_parabolic(0.3, 1, 1)
C5 = profiles.c_const(5.0)


class TestDistributionWeight:
    def test_jacobian_differs_from_thesis_by_jacobian_squared(self):
        """'jacobian' differs from 'thesis' by exactly (dy2/dkappa1)^2."""
        om, y = 2.0, np.array([0.0, 0.3, 0.8])
        dy_dk = -PARABOLA.f(y) ** 2 / (om * PARABOLA.df(y))
        assert np.allclose(distribution_weight(om, PARABOLA, C5, y, "jacobian"),
                           distribution_weight(om, PARABOLA, C5, y, "thesis") * dy_dk ** 2)

    def test_rejects_unknown_weighting(self):
        with pytest.raises(ValueError):
            distribution_weight(1.0, PARABOLA, C5, np.array([0.3]), "other")

    def test_jacobian_weighting_equals_y2_integral(self):
        """With 'jacobian', sum_j g w_j omega_c = int_0^delta c^2 omega_c / U' dy2
        (the literal source integral, up to the sign convention of dkappa1)."""
        om = 1.0
        k, w, y = kappa1_quad(om, PARABOLA, 1.0, 400)
        oc = lambda y: (1 - y) ** 2
        lhs = np.sum(distribution_weight(om, PARABOLA, C5, y, "jacobian") * oc(y) * w)
        yy, wy = q.lgwt(400, 0, 1)
        rhs = np.sum(25 * oc(yy) / PARABOLA.df(yy) * wy)
        assert np.isclose(abs(lhs), abs(rhs), rtol=1e-4)
