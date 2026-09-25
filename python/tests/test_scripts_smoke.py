"""Run every driver script end-to-end at low resolution.

These check that the scripts execute and return finite fields of the right
shape -- *not* that they reproduce the thesis figures (that needs the
original resolutions, and MATLAB reference output to compare against).
"""
import numpy as np
import pytest

from shear_flow_scattering import FIELDS
from shear_flow_scattering.scripts import diagnostics, point_source, scattering, trailing_edge

X1 = np.linspace(-2, 2, 9)
X2 = np.linspace(0, 1.5, 7)
SMALL = dict(N_bl=150, N_C=40)


def _finite(S, shape=None):
    for f in FIELDS:
        a = np.asarray(S[f])
        assert np.all(np.isfinite(a)), f
        if shape is not None:
            assert a.shape == shape, (f, a.shape)


def test_point_source_generic_and_hard():
    for bc in ("hard", "pressure", "impedance"):
        r = point_source.point_source_generic(bc=bc, x1=X1, x2=X2, N_CL=60, **SMALL)
        _finite(r.I_tot, (X1.size, X2.size))
    r = point_source.point_source_hard(x1=X1, x2=X2, N_CL=60, **SMALL)
    _finite(r.I_tot, (X1.size, X2.size))


def test_point_source_wake():
    x2 = np.linspace(-1.5, 1.5, 9)
    r = point_source.point_source_wake(x1=X1, x2=x2, N_CL=60, **SMALL)
    _finite(r.I_tot, (X1.size, x2.size))
    assert np.all(np.isfinite(r.dispersion_zeros))


def test_far_field_point_sources():
    r = point_source.far_field_generic_point_source(N_theta=15, N_bl=150)
    _finite(r.FF, (15,))
    r = point_source.far_field_wake_point_source(N_theta_p=8, N_theta_m=8, N_bl=150)
    # theta = 0, pi, 2pi put the saddle on the branch point where D = 0 for this
    # (uniform-flow) default: 0 * inf = NaN there, exactly as in MATLAB.
    interior = np.abs(np.sin(r.theta)) > 1e-9
    _finite({f: r.FF[f][interior] for f in FIELDS})


def test_pressure_zero_tracker():
    c0, z = point_source.pressure_zero_tracker(c0=[6, 7], N_bl=200)
    assert np.all(np.isfinite(z)) and np.all(z.imag > 0)


def test_hard_soft():
    r = scattering.hard_soft_continuous(x1=X1, x2=X2, prec_N=300, **SMALL)
    _finite(r.I_full, (X1.size, X2.size))
    r = scattering.far_field_hard_soft(N_theta=15, N_bl=150, prec_N=300)
    _finite(r.FF, (15,))


def test_hard_soft_hard():
    x1 = np.linspace(-7, 3, 11)
    r = scattering.hard_soft_hard(x1=x1, x2=X2, N_CL=60, prec_N=300, **SMALL)
    _finite(r.I_full, (x1.size, X2.size))
    r = scattering.far_field_hard_soft_hard(N_theta=15, N_bl=150, prec_N=300)
    _finite(r.FF, (15,))


def test_impedance_junction():
    r = scattering.imp_imp_continuous(x1=X1, x2=X2, N_CL=60, prec_N=300, **SMALL)
    _finite(r.I_full, (X1.size, X2.size))
    r = scattering.far_field_imp_imp(R=2.0, N_theta=15, N_bl=150, prec_N=300)
    _finite(r.FF, (15,))


def test_trailing_edge_lined():
    x2 = np.linspace(-1.5, 1.5, 9)
    for Z in (2 - 2j, np.inf):
        r = trailing_edge.trailing_edge_lined(Z=Z, x1=X1, x2=x2, N_CL=60, prec_N=300,
                                              omega=1.0, **SMALL)
        _finite(r.I_total, (X1.size, x2.size))


def test_far_field_trailing_edge():
    r = trailing_edge.far_field_trailing_edge_lined(N_theta_p=8, N_theta_m=8, N_bl=150, prec_N=300)
    assert all(np.all(np.isfinite(r.FF[f])) for f in ("phi", "p", "v"))
    r = trailing_edge.wall_pressure_far_field_trailing_edge_lined(N_theta_p=8, N_theta_m=8,
                                                                 N_bl=150, prec_N=300)
    assert r.FF.p.shape == (16,)
    ff = trailing_edge.wall_pressure_point_source_far_field(
        5.0, trailing_edge._study_profiles()[0], np.linspace(0.1, 0.9, 3), N_bl=150, prec_N=300)
    assert ff.shape == (3,) and np.all(np.isfinite(ff))


def test_distributed_trailing_edge_and_rerun():
    x2 = np.linspace(-1.2, 1.2, 9)
    r = trailing_edge.distributed_trailing_edge(x1=X1, x2=x2, N_distribution=6, N_CL=60,
                                                prec_N=300, **SMALL)
    _finite(r.I_total, (X1.size, x2.size))
    r2 = trailing_edge.rerun_with_vorticity(r.state, r.y2 * (1 - r.y2))
    _finite(r2.I_total, (X1.size, x2.size))
    assert not np.allclose(r.I_total.p, r2.I_total.p)


def test_wall_pressure_scaled():
    r = trailing_edge.wall_pressure_point_source_scaled(y2=[0.2, 0.4], omega=[1.0, 10.0], N_bl=200)
    assert np.all(np.isfinite(r.Full_Pressure_base))


def test_diagnostics():
    # I=0 (the MATLAB default) integrates straight through the critical point
    # U(x2) = omega/kappa1 and is singular there -- deform the path instead.
    g = diagnostics.scrappy_greens(x2=np.linspace(0.05, 0.95, 5), I=1, N_bl=200)
    assert np.all(np.isfinite(g.W))
    w = diagnostics.scrappy_wronskian(shift=0.01, y2=[0.2, 0.5], N_bl=200)
    assert np.all(np.isfinite(w.W))
