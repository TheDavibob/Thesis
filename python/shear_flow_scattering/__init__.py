"""NumPy port of the MATLAB code for the thesis *Acoustic Scattering in Shear Flow*.

Module map (MATLAB directory -> Python module):

=====================================  =========================================
``CommonTools`` (logs, powers)         :mod:`.complex_funcs`
``CommonTools/GaussQuad``, ``lgwt``    :mod:`.quadrature`
``CommonTools/Integration``,
``GenericInversion/Contours_Weights``  :mod:`.contours`
``CommonTools/Two_Solutions_Complex``  :mod:`.frobenius` (partly reconstructed)
``CommonTools/PrincipalValue``,
``ZP finder``                          :mod:`.zeros`
``FixedStep/Flow_Structures``          :mod:`.profiles`
``FixedStep/Fixed_Step_IVP*``          :mod:`.ivp`
``FixedStep/PridmoreBrown``            :mod:`.pridmore_brown`,
                                       :mod:`.steepest_descent`
``Factorisation/Fact_*``               :mod:`.factorisation`
``Factorisation/*Kernel*``             :mod:`.kernels`
``GenericInversion`` (functions)       :mod:`.inversion`, :mod:`.distributed`
``GenericInversion`` (scripts)         :mod:`.scripts`
=====================================  =========================================

See ``README.md`` for missing files, bugs found and fixed, and differences
from MATLAB behaviour.
"""
from ._utils import FIELDS, Path, Struct
from . import (complex_funcs, contours, distributed, factorisation, frobenius, inversion, ivp,
               kernels, pridmore_brown, profiles, quadrature, steepest_descent, zeros)

__all__ = ["FIELDS", "Path", "Struct", "complex_funcs", "contours", "distributed",
           "factorisation", "frobenius", "inversion", "ivp", "kernels", "pridmore_brown",
           "profiles", "quadrature", "steepest_descent", "zeros"]
