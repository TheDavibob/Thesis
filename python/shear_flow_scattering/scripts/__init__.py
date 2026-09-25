"""The MATLAB driver *scripts*, rewritten as functions with keyword defaults.

* :mod:`.point_source`  -- ``GenericInversion/PointSource_Cts`` (no scattering)
* :mod:`.scattering`    -- hard/soft, hard/soft/hard, impedance junctions
* :mod:`.trailing_edge` -- lined trailing edge, distributed vorticity,
  wall-pressure studies
* :mod:`.diagnostics`   -- the "Scrappy" Wronskian checks

Most are expensive at their original default resolutions (minutes to hours);
pass coarser grids / fewer quadrature points to experiment.
"""
