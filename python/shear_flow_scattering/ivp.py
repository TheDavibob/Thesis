"""Fixed-step fourth-order Runge-Kutta, vectorised over a parameter ``k``.

Ports of ``FixedStep/Fixed_Step_IVP*.m``.  Solves ``dphi/dx = f(x, phi; k)``
from ``x = 0`` to ``x_max`` for *many* parameter values at once.

Shape conventions (kept from MATLAB):

* ``phi0(k)`` returns an array of shape ``(phi_N, k_N)``;
* ``f(x, phi, k)`` takes/returns ``(phi_N, k_N)`` arrays;
* the full solution has shape ``(x_N, phi_N, k_N)``.

The ``_skip`` variants keep only the end point (much less memory).  The
``_gridx`` variants let every ``k`` have its own end point ``x_max[k]``
(each column uses its own step ``H[k]``; ``x`` passed to ``f`` is then an
array of length ``k_N``).
"""
from __future__ import annotations

from typing import Callable

import numpy as np


def _rk4_step(f, x, phi, h, k):
    c1 = h * f(x, phi, k)
    c2 = h * f(x + h / 2, phi + c1 / 2, k)
    c3 = h * f(x + h / 2, phi + c2 / 2, k)
    c4 = h * f(x + h, phi + c3, k)
    return phi + (c1 + 2 * c2 + 2 * c3 + c4) / 6


def fixed_step_ivp(f: Callable, phi0: Callable, x_max, x_N: int, k):
    """``Fixed_Step_IVP.m``.  Returns ``(phi, x)`` with ``phi.shape == (x_N, phi_N, k_N)``."""
    x = np.linspace(0, x_max, x_N)
    h = x[1] - x[0]
    k = np.ravel(k)
    IC = np.asarray(phi0(k), dtype=complex)
    phi = np.zeros((x_N,) + IC.shape, dtype=complex)
    phi[0] = IC
    for j in range(1, x_N):
        phi[j] = _rk4_step(f, x[j - 1], phi[j - 1], h, k)
    return phi, x


def fixed_step_ivp_skip(f: Callable, phi0: Callable, x_max, x_N: int, k):
    """``Fixed_Step_IVP_skip.m``: as :func:`fixed_step_ivp`, returns only ``phi(x_max)``."""
    x = np.linspace(0, x_max, x_N)
    h = x[1] - x[0]
    k = np.ravel(k)
    phi = np.asarray(phi0(k), dtype=complex)
    for j in range(1, x_N):
        phi = _rk4_step(f, x[j - 1], phi, h, k)
    return phi


def _grid(x_max, x_N, k):
    K = np.ravel(k)
    X_max = np.ravel(x_max)
    if X_max.size != K.size:
        raise ValueError("Discrepancy between number of parameter values and end point values")
    X = X_max[:, None] * np.linspace(0, 1, x_N)[None, :]  # (k_N, x_N)
    return K, X, X[:, 1] - X[:, 0]


def fixed_step_ivp_gridx(f: Callable, phi0: Callable, x_max, x_N: int, k):
    """``Fixed_Step_IVP_GridX.m``: per-``k`` end points.

    Returns ``(phi, X, phi_end)``; ``X`` has shape ``(k_N, x_N)``,
    ``phi`` ``(x_N, phi_N, k_N)`` and ``phi_end = phi[-1]``.
    """
    K, X, H = _grid(x_max, x_N, k)
    IC = np.asarray(phi0(K), dtype=complex)
    phi = np.zeros((x_N,) + IC.shape, dtype=complex)
    phi[0] = IC
    for j in range(1, x_N):
        phi[j] = _rk4_step(f, X[:, j - 1], phi[j - 1], H, K)
    return phi, X, phi[-1].copy()


def fixed_step_ivp_gridx_skip(f: Callable, phi0: Callable, x_max, x_N: int, k):
    """``Fixed_Step_IVP_GridX_skip.m``: per-``k`` end points, end values only."""
    K, X, H = _grid(x_max, x_N, k)
    phi = np.asarray(phi0(K), dtype=complex)
    for j in range(1, x_N):
        phi = _rk4_step(f, X[:, j - 1], phi, H, K)
    return phi


# ---------------------------------------------------------------------------
# FixedStep/Test -- damped oscillator y'' + 2k y' + y = 0, y(0)=1, y'(0)=0
# ---------------------------------------------------------------------------
def _oscillator(x, phi, k):
    return np.stack([phi[1], -phi[0] - 2 * k * phi[1]])


def _oscillator_ic(k):
    return np.stack([np.ones(k.size), np.zeros(k.size)])


def fixed_step_ivp_test(x_max, x_N, k):
    """``Fixed_Step_IVP_test.m``."""
    return fixed_step_ivp(_oscillator, _oscillator_ic, x_max, x_N, k)


def fixed_step_ivp_test_gridx(x_max, x_N, k):
    """``Fixed_Step_IVP_test_GridX.m``."""
    return fixed_step_ivp_gridx(_oscillator, _oscillator_ic, x_max, x_N, k)


def fixed_step_ivp_test_gridx_skip(x_max, x_N, k):
    """``Fixed_Step_IVP_test_GridX_skip.m``."""
    return fixed_step_ivp_gridx_skip(_oscillator, _oscillator_ic, x_max, x_N, k)
