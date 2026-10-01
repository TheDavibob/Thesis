# shear-flow-scattering

A NumPy port of the MATLAB code in `../code` for the thesis *Acoustic
Scattering in Shear Flow* (see `../Thesis_Final.pdf`).

The physics: a disturbance convecting in a sheared boundary layer `U(x2)`
scatters off a change in wall condition (hard → pressure release, impedance
junctions, a trailing edge into a wake). The code solves the Pridmore-Brown
equation across the boundary layer with a fixed-step RK4, bending the
integration path into the complex `x2`-plane to go round the critical layer
`U(x2) = omega/k1`. It factorises the Wiener-Hopf kernel numerically, then
inverts the Fourier transform along hyperbolic "steepest-descent-like"
contours, around the critical-layer branch cut, and through residues at
instability modes.

Every MATLAB file has a Python counterpart (table below), except where
dependencies are missing from the repository (listed under
[What's missing](#whats-missing-from-the-matlab-repository)).

## Quick start

```bash
cd python
uv sync                      # numpy + pytest; add `--extra plot` for matplotlib
uv run pytest                # ~2 min: known-answer tests + low-res runs of every script
```

```python
import numpy as np
from shear_flow_scattering.scripts import scattering

# FarField_HardSoft_Continuous_Nov2018.m, default parameters
r = scattering.far_field_hard_soft()
np.abs(r.FF.p)               # far-field pressure directivity over r.theta
```

`examples/fig_6_2_2.py` regenerates thesis Fig. 6.2.2, the hard trailing
edge. Its docstring explains how the unstated parameters were pinned down
from the wake mode the thesis quotes.

The only runtime dependency is NumPy (SciPy isn't needed). Matplotlib is
optional and only used by the few functions that took a `plot` flag.
Progress messages that were `disp`/`waitbar` in MATLAB now go through the
`logging` module (`logging.basicConfig(level=logging.INFO)` to see them).

## Layout

| MATLAB | Python |
|---|---|
| `lgwt/lgwt.m`, `CommonTools/GaussQuad/*` | `quadrature.py` (`lgwt`, `gauss_legendre/hermite/laguerre`, plus `adaptive_integral` replacing MATLAB's `integral`) |
| `CommonTools/LogA, SqrtA, PowerA, PowerB, Log_CurvedBranch, LogCts, LogCtsF, NaN2Zero` | `complex_funcs.py` |
| `CommonTools/heaviside.m` | `_utils.heaviside` |
| `CommonTools/Integration/*Contour*.m`, `Truncate.m` | `contours.py` |
| `CommonTools/Integration/AddFactB.m` | `factorisation.add_fact_b` (stub, dependencies missing) |
| `CommonTools/PrincipalValue/PVTrapz.m`, `Poles.m` | `zeros.pv_trapz`, `zeros.poles` |
| `CommonTools/Two_Solutions_Complex.m`, `param_setup.m` | `frobenius.py` (helpers reconstructed) |
| `Factorisation/Fact_*.m`, `Prec_Setup.m`, `Contour_Setup.m` | `factorisation.py`, `contours.default_factorisation_contour` |
| `Factorisation/{Generic_BC,Generic_Wake,NoSlip_HardSoft}/*.m` | `kernels.py` |
| `FixedStep/Fixed_Step_IVP*.m`, `FixedStep/Test/*` | `ivp.py` |
| `FixedStep/Flow_Structures/*.m` | `profiles.py` (`Struct_U_Parabolic` → `u_parabolic` etc.) |
| `FixedStep/Common_Functions/Gamma_FF.m` | `pridmore_brown.gamma_ff` |
| `FixedStep/Common_Functions/CL_Ellipse_Contour.m` | `contours.cl_ellipse_contour` |
| `FixedStep/PridmoreBrown/*.m`, `Complex_Contour/*.m` | `pridmore_brown.py` (`Dec_Aug18A` → `dec_aug18a`, `Comp_Dec_Aug18A_skip` → `comp_dec_aug18a_skip`, ...) |
| `FixedStep/PridmoreBrown/Steepest_Descent_Contour*.m`, `Sep_18_Inversion/*` | `steepest_descent.py` |
| `FixedStep/PridmoreBrown/NumericalAccuracy.m` | `pridmore_brown.numerical_accuracy` |
| `GenericInversion/Fourier_Nov2018.m`, `Integral_Nov2018.m` | `inversion.fourier_nov2018`, `inversion.integral_nov2018` |
| `GenericInversion/Contours_Weights/*.m` | `contours.py` |
| `GenericInversion/PointSource_Cts/*` (functions) | `inversion.py` |
| `GenericInversion/Scattering/Y_*.m`, `Dirty_Dispersion.m` | `inversion.y_continuous_scattering` (the four variants merged), `inversion.dirty_dispersion` |
| `GenericInversion/Scattering/Distributed_Source/*` (functions) | `distributed.py` |
| `PointSource_*_Nov18.m`, `PointSource_Wake_Feb19.m`, `PressureZeroTracker.m`, `FarField/*` (scripts) | `scripts/point_source.py` |
| `HardSoft*`, `HardSoftHard*`, `ImpImp*`, and their `FarField_*` versions | `scripts/scattering.py` |
| `TrailingEdge_Lined.m`, `FarField_TrailingEdge_Lined.m`, `WallPressure_FarField_TrailingEdge_Lined.m`, `Distributed_TrailingEdge*.m`, `WallPressure_PointSource_*.m` | `scripts/trailing_edge.py` |
| `Scrappy_Greens.m`, `Scrappy_Wronskian.m` | `scripts/diagnostics.py` |
| `ZP finder/*` | `zeros.py` |
| `SmoothParabolaPlot.m` | `profiles.smooth_parabola_plot` |
| `startup.m`, `contents.m` | not needed (package imports) |

### Conventions carried over from MATLAB

* **Structs.** `U.f`, `dec.phi`, `C.k`/`C.wk` and the like are kept as attribute
  access on `_utils.Struct` (a `SimpleNamespace` that also allows
  `s["phi"]`, like MATLAB's `s.(name)`). Every solver returns `phi`, `dphi`,
  `p`, `v`.
* **`lgwt` ordering.** `quadrature.lgwt` returns nodes **descending**, like
  MATLAB's. The original code reverses them in many places, and the port keeps
  those reversals, so this ordering must not change.
* **Array shapes.** RK4 solutions are `(x_N, phi_N, k_N)`, `Y` fields are
  `(len(x2), len(k))`, and assembled fields are `(len(x1), len(x2))`, all as
  in MATLAB.
* **Scripts became functions.** Every hard-coded parameter at the top of a
  MATLAB script is now a keyword argument with the same default. The results
  come back in a `Struct` instead of being left in the workspace.

### Deliberate structural changes

* **Vectorisation.** The `Y_*` functions looped over `x2`, calling the
  boundary-layer solver once per row. Now every `(x2, k1)` pair runs in one
  RK4 sweep, each with its own integration path. The arithmetic is the same
  (the tests check it against the row-by-row version). `max_batch` limits
  memory use. `Integral_Nov2018`, `Upstream_Distributed_Integrated` and
  `Distributed_WallPressure` are vectorised as well.
* **Path helpers.** The many hand-written `z.f`/`z.df` semicircles
  (`zdy`, `zly`, `Contour_Semicircle`, ...) are all instances of
  `contours.semicircle_path(a, b, half_plane)`. `contours.path_to` wraps the
  repeated `if HP == 0 straight else semicircle` block.
* **Saddle points.** The far-field scripts built an `N_theta × N_theta`
  meshgrid only to take its diagonal (`diag(k1.Theta)`).
  `steepest_descent.far_field_saddle` computes the diagonal directly.
* **Duplicated scripts merged.** `WallPressure_FarField_TrailingEdge_Lined`
  is a thin wrapper round `far_field_trailing_edge_lined`, the four
  `Y_Continuous_Scattering_*` files are one function with `HP`/`truncate`
  flags, and `Distributed_TrailingEdge_VarOmegaC.m`, which relied on the
  workspace, is now `rerun_with_vorticity(state, omega_c)`.
* **Analytic continuation** (`Fact_Mult_Cont*`) evaluates only the branch it
  needs at each point. MATLAB evaluated both and blended them with 0/1 masks,
  which let a `NaN` from the unused branch leak through (`NaN*0 = NaN`), and
  which also called the expensive kernel twice.

## What's missing from the MATLAB repository

Code that is called but isn't in the repo:

| Missing | Called from | What the port does |
|---|---|---|
| `FullIntervalIntA.m`, `CauchyKernel.m` | `CommonTools/Integration/AddFactB.m` | `add_fact_b` raises `NotImplementedError`. The routine is legacy (2017); `fact_additive` superseded it. |
| `ContourSweep.m`, `BranchCuts.m` | `ZP finder/ZerosPolesBranchCuts.m` | `zeros_poles_branch_cuts` raises `NotImplementedError`. (The original also has `numel(B,1)`, line 24, which isn't valid for the purpose.) `branch_cut_integral` is ported. |
| `TaylorCoeff_Complex.m`, `Two_Series.m`, `Taylor_Series.m` | `CommonTools/Two_Solutions_Complex.m` | **Reconstructed** in `frobenius.py` (standard Frobenius method, with the log term for repeated or integer-separated roots). Tested against Bessel `J0` and ODE residuals. The meaning of the returned type flag `I` and the free constant `d_m` are guesses. Nothing else calls this function. |
| `PlotRI.m` | `WallPressure_PointSource_Scaled.m` | A small replacement (`trailing_edge._plot_ri`, real parts solid and imaginary parts dashed). |

Inputs that the scripts took from the MATLAB workspace (never defined in the
file itself):

| Variable | Script | Port |
|---|---|---|
| `R` | `FarField_ImpImp_Continuous_Feb2019.m:16` | required argument of `far_field_imp_imp` |
| `U0` (line commented out) | `FarField_TrailingEdge_Lined.m:13` | defaults to 0.3, the value in the commented line |
| `shift`, `I`, `delta`, `U`, `c`, `omega`, `k3`, `N_bl`, `l1`, `l2` | `Scrappy_Wronskian.m` | `Scrappy_Greens` defaults; `shift` is required |
| `Omega`, `U`, `y2`, `c0`, `delta` | `WallPressure_PointSource_FarField.m` (set by its `_routine`) | explicit arguments |
| everything | `Distributed_TrailingEdge_VarOmegaC.m` | `rerun_with_vorticity(state, ...)` |
| `Velocity_Upstream` (when `Z = 0`) | `FarField_TrailingEdge_Lined.m` | raises `NotImplementedError` for `Z = 0` |

Unfinished or placeholder parts of the original, ported as they are:

* `Sep18_Integrate_SD.m` is headed "UNFINISHED: A HECK OF A LOT OF AWKWARD WORK
  IS NEEDED". Its Wronskian normalisation is marked "not yet correct", and
  the hard-wall solution above the boundary layer is "FUTURE WORK".
* `Steepest_Descent_Contour_Anon.m` is described as "Not finished, essentially a
  check the mathematics is correct".
* In the point-source scripts the scaling constant `A = 1` is "to be fixed"
  and never applied.
* The residue block for dispersion zeros in `HardSoft_Continuous_Nov2018.m` is
  commented out in the original and isn't ported.
* `FFKernel_Impedance_Historic.m` says "THIS IS WRONG WRONG WRONG DO NOT USE".
  It is ported with a warning.
* The `.mat` files (`WallPressure_HardWall.mat`, `Wall_FarField_HardWall.mat`)
  are saved results that no code reads, and weren't converted. `.fig`, `.png`
  and `.tex` outputs are ignored. `CommonTools/LogCts.m~` and `LogCtsF.m~` are
  editor backups: they only differ in comments, one of which reads "WORTH
  NOTING I'M REALLY CONFUSED".

## Bugs found in the MATLAB, and fixed in the port

Each fix is also flagged in the relevant docstring.

| Where (MATLAB) | Problem | Fix |
|---|---|---|
| `FFKernel_HardSoft_NoSlip.m:10` | *Porting* trap rather than a MATLAB bug: `log((sigma+omega)/(omega-sigma))` is negative whenever `U'(0) > omega`. MATLAB returns a complex log; `np.log` returns NaN, which made every no-slip hard→soft factorisation NaN. | `np.log(complex(...))` (same in the historic impedance kernel). A regression test covers it. |
| `ZeroFinder_Recursive.m:178` | The lower-half-plane check `imag(z) < prec_real` is true for **every** LHP zero, so they were all rejected and the routine errored after 10 nudges. | Rejects `imag(z) > -prec_real`. |
| `ZeroFinder_Recursive.m:85` (and similar) | After deflation `F.df` was rebuilt using the *already deflated* `F.f`, so the derivative was wrong. Newton still often converged, just slowly. | The correct deflated derivative is used. |
| `Zeros_Box.m:23` | `N = varargin(1)` gives a cell array, not a number. | Fixed (`N` is a keyword argument). |
| `HardSoftHard_Continuous_Feb2019.m:287, 319` | `I_s1 + I_s2_1 + I_s2_1`: drops the downstream forcing term `I_s2_2` and counts `I_s2_1` twice (both the acoustic and the critical-layer sums). | `I_s1 + I_s2_1 + I_s2_2`. |
| `Dec_Aug18A_GridX.m:14` | `phi` was reversed along the grid but `dec.x` wasn't, so `phi`, `p` and `v` were paired with the wrong `x2`. | `dec.x` is reversed too (a test checks it). |
| `Oct_18_FF_SD.m:28` | The path `z = delta(1 - t)` had `dz = -delta*t`; it should be `-delta`. | Fixed. |
| `EllipseContour.m:13` | `dCp = 1 + dg` is missing the factor `1i`. | Fixed. |
| `TanhContourB.m:13` | The left-half derivative `(bm + a) sech^2/tm` should be `-(bm - a) sech^2/tm`. | Fixed. |
| `Y_GreensWake.m:44` | The loop indices `jmin + 1 + (1:N_ybl)` are only correct if `x2` contains `-deltam` exactly (the default grid happened to). | Masks. |
| `Y_GreensWake.m:113` | `x2(-x2 >= deltap)` where `deltam` was meant; a size mismatch whenever `deltap ≠ deltam`. | Fixed. |
| `Upstream_Distributed_B.m:45, 49` | `numel(x2_y2 == 1)` should be `numel(x2_y2) == 1`, and the `x2 == y2` branch reads `dec_y2.v`, which `_A` never stores (it crashes if a grid point hits a source height). | `_A` now stores `v`. |
| `Integrand_Distributed.m:4` | Uses `y2` without it being an argument. | `y2` added as an argument. |
| `Distributed_TrailingEdge.m:121` | `min(omega/Up.f(deltap), omega/Up.f(deltap))` uses `Up` twice; `Um` was meant for the second (as in the lined version). | Fixed. |
| scripts, e.g. `HardSoft_Continuous_Nov2018.m:24-25`; `FFKernel_HardSoft_Slip.m` | The branch points use `sqrt(k0^2 + k3^2(1-M^2))`; the zeros of `gamma` need `-k3^2`. Only matters for `k3 ≠ 0`, and every script has `k3 = 0`. | Fixed in `pridmore_brown.branch_points`. Left in `ff_kernel_hardsoft_slip`, where it only enters when `W ≠ 0`. |
| `NewtonRaphson.m` | Has no iteration cap, so it can loop forever. | `max_iter=1000`, with a warning. |

## Suspicious, but left as the original has it

* **Outer-region velocity.** Outside the boundary layer, `Y.v = -C^2 phi'`
  drops the `-3 C C' phi` term, i.e. it assumes `U'(delta) = 0`. That holds
  for the parabolic profiles but **not** for `Struct_U_Linear`, which
  `TrailingEdge_Lined` and `PointSource_Hard_Nov18` use by default. `v` is
  then discontinuous at `x2 = delta` by `3 C U'(delta) k1 phi` (in magnitude).
* **`WallPressure.m`.** `Q0` contains `(omega - U(y2) kappa1)^4` with
  `kappa1 = omega/U(y2)`, so it is identically zero. The function is
  superseded by `Distributed_WallPressure`.
* **Wake far field.** `FarField_WakePointSource_Feb19.m` evaluates the lower
  half (`theta ∈ [pi, 2 pi]`) with the upper flow (`Up`, `cp`, `deltap`). It
  makes no difference for the symmetric default.
* **Plate BC sign in `Distributed_TrailingEdge.m`.** It uses
  `L2 = Zed V2 - C0 P` with `DZ2 = -D(L2)`. For a hard plate (the only case
  it supports) that equals `TrailingEdge_Lined`/`Wake_Kernel`, but not for
  finite `Z`. The port uses the `Wake_Kernel` convention throughout.
* **`Z = 0` plates.** `A = I0/(1i*Zed*Km)` divides by `Zed = 0`, so the
  trailing-edge scripts can't handle a pressure-release plate. The port raises
  an error.
* **Branch points at the far-field endpoints.** At `theta = 0, pi` the saddle
  sits on the acoustic branch point, and where the dispersion relation
  vanishes there the result is `0 × ∞ = NaN` in MATLAB too. With the default
  (uniform-flow) wake this happens at the endpoints.
* **`Scrappy_Greens` with `I = 0`** (its default) integrates straight through
  the critical point `U(x2) = omega/kappa1` and is singular there.

* **Single-sheet incident field (all `Y_Greens*` at `k1 = kappa1`).** These
  evaluate `phi` *at* the sheet `y2`, which is the critical level for
  `kappa1 = omega/U(y2)`, where `phi ~ C^-3` is infinite. The result is a
  resolution-dependent number that grows like `N_bl^3` and is absorbed into
  `Q0 = 1e-10`. The paths above and below the sheet have lengths `delta - y2`
  and `y2`, so the incident pressure jumps across the sheet by about
  `((delta - y2)/y2)^3`: a factor of about 700 at `y2 = 0.1`, the default of
  `HardSoft_Continuous`, and 3.4 at `y2 = 0.4` (`HardSoftHard`). It is only
  consistent at `y2 = delta/2`, as in Fig. 6.2.2. The scattered field is
  consistent with the wall value, so only the incident part above the sheet is
  affected. The regular form uses one-sided limits `p = lim(-C^3 phi)`, as
  `Distributed_WallPressure` does with `shift_y2`. Not yet changed in the port.
* **Critical-layer branch for `D(kappa1)`.** At `kappa1 = omega/U(y2)`, the
  dispersion function needs `phi_d` continued through the source's own
  critical level, and the two ways round give complex conjugates. Causality
  (`Im omega < 0` at real `kappa1`) selects the upper-half-plane path
  (`HP = +1`); `D(kappa1 + i eta)` converges to it. `Distributed_TrailingEdge`
  uses `+1`, but `WallPressure_PointSource_Scaled` (Fig. 6.4.2) and
  `FarField_TrailingEdge_Lined` use `-1`. The moduli agree, so `|p_w|` is
  unaffected; the phases are conjugated.

## Distributed sources: the `weighting` toggle

`distributed_trailing_edge`, `rerun_with_vorticity` and
`upstream_distributed_integrated` take `weighting="thesis"` (the default:
eq. 6.3.10 and the MATLAB) or `"jacobian"` (the literal `y2 -> kappa1` change
of variables). The two differ by `(dy2/dkappa1)^2`: (6.3.10) applies the
Jacobian upside down, which weights near-wall vorticity up by about
`(omega U'/U^2)^2`, roughly 240 at the wall for the default profile. Both
weightings converge in `N_distribution`. The `"jacobian"` weight grows like
`1/U'^2` towards the edge, but the unit-source wall pressure `I0` vanishes
there like `(delta - y2)^4`. `examples/fig_weighting_comparison.py` plots
the two side by side.

## Other differences from MATLAB behaviour

* **Points exactly on a branch cut.** NumPy honours signed zeros
  (`sqrt(-4-0j) = -2j`) where MATLAB gives `+2i`. Values on a cut can differ;
  values off the cuts agree.
* **`integral`.** Replaced by adaptive Gauss-Legendre with MATLAB's default
  tolerances (AbsTol 1e-10, RelTol 1e-6). Results agree to tolerance, not bit
  for bit.
* **`interp1(..., 'pchip')`.** Replaced by a Fritsch-Carlson PCHIP
  (`_utils.pchip`), the same algorithm MATLAB uses.
* **`Zeros_Box`** returns `0` rather than `[]` for half-planes the box doesn't
  cover.

## Testing and what has been verified

`uv run pytest` runs 47 tests.

* **Known answers for the numerics** (`tests/test_core.py`):
  * Gauss rules and adaptive integration.
  * Branch-cut functions.
  * RK4 against the damped oscillator.
  * Profile derivatives against finite differences.
  * `gamma^2` against the dispersion relation.
  * **Uniform flow:** Pridmore-Brown reduces to `phi'' = gamma^2 phi`, so the
    decaying solution must be exactly `exp(-gamma x2)`, along straight and
    semicircular paths.
  * Path independence away from the critical layer.
  * The GridX bug fix.
  * Additive and multiplicative Wiener-Hopf splits of functions with known
    factors, including the analytic continuation.
  * **The actual hard→soft kernel:** `K+ K- = K` between separated contours,
    to about `1e-6` with slip and `1e-5` without (limited by `N_bl`).
  * Zero finding, Frobenius vs `J0`, vectorised vs row-by-row `Y`, and
    continuity of the Green's function at the source.
* **Every script runs end to end** at low resolution and gives finite fields
  of the right shape (`tests/test_scripts_smoke.py`).

**Not verified:** there are no MATLAB reference outputs in the repository, so
nothing has been compared number-for-number against the original, and the
thesis figures haven't been regenerated. The original default resolutions
are expensive (`N_bl` up to 10 000, `prec.N` up to 10 000, 400 × 400 grids);
expect minutes to hours per full-resolution script. To validate the port,
run the same parameters through MATLAB and compare, starting with
`far_field_hard_soft` and `hard_soft_continuous`, which exercise most of the
pipeline.
