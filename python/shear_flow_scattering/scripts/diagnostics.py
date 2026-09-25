"""Exploratory checks: Wronskians of the boundary-layer solutions.

Ports of ``Distributed_Source/Scrappy_Greens.m`` and ``Scrappy_Wronskian.m``.
(``NumericalAccuracy.m`` is :func:`..pridmore_brown.numerical_accuracy`;
``SmoothParabolaPlot.m`` is :func:`..profiles.smooth_parabola_plot`.)

For the decaying solution ``phi_d`` and the wall solution ``phi_l`` the
Wronskian ``W = phi_d' phi_l - phi_d phi_l'`` should vary smoothly with
height; these scripts tabulate it, which is how the normalisation of the
Green's function was checked.
"""
from __future__ import annotations

import numpy as np

from .._utils import Struct, row
from ..contours import path_to
from ..pridmore_brown import comp_dec_aug18a_skip, comp_generic_aug18a_skip
from ..profiles import c_const, u_parabolic
from .common import wall_operators


def _setup(U0, c0, omega):
    U, c = u_parabolic(U0, 1, 1), c_const(c0)
    return U, c, wall_operators(omega, U).hard


def _wronskian(x2, kappa1, shift_y2, I, U, c, L, delta, N_bl, omega, k3):
    x2 = row(x2)
    kappa1 = np.broadcast_to(kappa1, x2.shape).astype(complex)
    dec = comp_dec_aug18a_skip(kappa1, k3, omega, U, c, delta, N_bl, path_to(delta, x2 + shift_y2, I))
    inner = comp_generic_aug18a_skip(kappa1, k3, omega, U, c, L.l1, L.l2, N_bl,
                                     path_to(0.0, x2 - shift_y2, I))
    return Struct(x2=x2, A=dec.phi, B=dec.dphi, C=inner.phi, D=inner.dphi,
                  W=dec.dphi * inner.phi - dec.phi * inner.dphi)


def scrappy_greens(x2=None, y2=0.5, I=0, U0=0.3, c0=5.0, delta=1.0, N_bl=2000, omega=1.0,
                   k3=0.0, shift_y2=0.0):
    """``Scrappy_Greens.m`` -- ``phi_d``, ``phi_l`` and their Wronskian across
    the layer at fixed ``kappa1 = omega / U(y2)`` (hard wall).

    Note: for ``I != 0`` *and* ``shift_y2 != 0`` the MATLAB semicircle did
    not start exactly at ``delta``; here it does.  Identical for the default
    ``shift_y2 = 0``.
    """
    U, c, L = _setup(U0, c0, omega)
    x2 = np.arange(0, 1 + 1e-9, 0.01) if x2 is None else x2
    return _wronskian(x2, omega / U.f(y2), shift_y2, I, U, c, L, delta, N_bl, omega, k3)


def scrappy_wronskian(shift, y2=None, I=0, U0=0.3, c0=5.0, delta=1.0, N_bl=2000, omega=1.0,
                      k3=0.0, shift_y2=0.001):
    """``Scrappy_Wronskian.m`` -- the Wronskian *at the source* for each
    source height, with ``kappa1 = omega/U(y2) + i shift``.

    .. warning:: **Missing inputs in the original.**  The script used ``I``,
       ``delta``, ``U``, ``c``, ``omega``, ``k3``, ``N_bl``, ``l1``, ``l2``
       and ``shift`` from the workspace (presumably left by
       ``Scrappy_Greens``).  Those defaults are used here; ``shift`` (the
       imaginary offset of ``kappa1``) was never set anywhere and is required.
    """
    U, c, L = _setup(U0, c0, omega)
    y2 = np.arange(0.01, 0.99 + 1e-9, 0.01) if y2 is None else row(y2)
    kappa1 = omega / U.f(y2) + 1j * shift
    return _wronskian(y2, kappa1, shift_y2, I, U, c, L, delta, N_bl, omega, k3)
