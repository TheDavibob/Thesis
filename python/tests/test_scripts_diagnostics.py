"""Diagnostic scripts (``scripts.diagnostics``)."""
import numpy as np

from shear_flow_scattering.scripts import diagnostics


class TestScrappyGreens:
    def test_finite_on_deformed_path(self):
        # I=0 (the MATLAB default) integrates straight through the critical point
        # U(x2) = omega/kappa1 and is singular there -- deform the path instead.
        g = diagnostics.scrappy_greens(x2=np.linspace(0.05, 0.95, 5), I=1, N_bl=200)
        assert np.all(np.isfinite(g.W))


class TestScrappyWronskian:
    def test_finite(self):
        w = diagnostics.scrappy_wronskian(shift=0.01, y2=[0.2, 0.5], N_bl=200)
        assert np.all(np.isfinite(w.W))
