"""Helpers shared by the test modules."""
import numpy as np

from shear_flow_scattering import FIELDS


def assert_finite(S, shape=None):
    """Every field of ``S`` (phi, dphi, p, v) is finite, and has ``shape`` if given."""
    for f in FIELDS:
        a = np.asarray(S[f])
        assert np.all(np.isfinite(a)), f
        if shape is not None:
            assert a.shape == shape, (f, a.shape)


def hard_wall_operators(omega, U):
    """(l1, l2, C0) for the hard-wall condition v = 0 at x2 = 0."""
    C0 = lambda k: 1j * (omega - U.f(0) * k)
    dC0 = lambda k: -1j * U.df(0) * k
    return (lambda k: -C0(k) ** 2), (lambda k: -3 * C0(k) * dC0(k)), C0
