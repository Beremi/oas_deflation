# TS-N65 Nonlinear Solver Experiment Handoff

Generated: 2026-05-22

Repository: `/home/beremi/repos/oas_deflation`

Branch at handoff: `nonlinear_solver_testing`

Previous checkpoint commit before CP6 queue: `317c23138dc7eab098816e1c16ac9619e00089b2`

## Purpose

This handout summarizes the TS-N65 nonlinear solver experiments performed so
far, including what was tested, why it was tested, what happened, and why each
non-promoted branch failed. It is intended for an expert reviewer in nonlinear
finite element/discrete element solution strategies, path-dependent softening,
damage/plasticity tangents, and continuation methods.

The short conclusion is:

- The strict TS-N65 local baseline is reproducible and remains the best complete
  reference run.
- No tested nonlinear globalization, indirect control, arc-length, or tangent
  branch has produced a valid speedup of the strict baseline.
- The most useful positive finding is diagnostic: the CSL active-damage tangent
  defect is real and can be fixed locally, but the locally better tangent does
  not by itself close the full TS-N65 phase benchmark.
- The dominant cost remains late-step load-control convergence in steps 6-8.
  Some methods reduce step 6, but then worsen or fail step 7.

## Baseline Definition

The current source of truth is:

```text
docs/tsn65-baseline-replication.md
```

Strict local baseline settings:

```text
time_step = min_time_step = max_time_step = 1.25e-3
total_time = 1e-2
stiff_matrix_type tangent
stiffness_matrix_iter_update 10
stiffness_matrix_step_update 1
solver_type DeflatedFGMRES
full DFGMRES/HYPRE block active
limit_tolerance 0
OMP_NUM_THREADS=16
```

Target result:

```text
steps: 8/8
iteration sequence: 6,6,10,13,17,183,187,163
total nonlinear rows: 585
no NaNs
no fallback acceptance
no nonlinear convergence failure
```

The full DFGMRES/HYPRE block is part of the baseline. A reduced solver deck does
not reproduce the baseline and can fail with a different sequence and DFGMRES
true-tolerance warnings.

Baseline step table from local replication:

| step | time | rows | residual | displacement | energy |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.00125 | 6 | 7.705844e-06 | 6.743608e-04 | 7.134058e-05 |
| 2 | 0.0025 | 6 | 4.198577e-04 | 1.021696e-04 | 6.255293e-04 |
| 3 | 0.00375 | 10 | 7.308889e-04 | 1.346596e-04 | 9.259996e-04 |
| 4 | 0.005 | 13 | 6.064366e-04 | 1.451312e-04 | 8.371949e-04 |
| 5 | 0.00625 | 17 | 6.189476e-04 | 1.136484e-04 | 9.322705e-04 |
| 6 | 0.0075 | 183 | 6.428869e-04 | 8.777311e-05 | 9.986343e-04 |
| 7 | 0.00875 | 187 | 6.721379e-04 | 1.022238e-04 | 9.840756e-04 |
| 8 | 0.01 | 163 | 5.838771e-04 | 1.075313e-04 | 9.817164e-04 |

Key interpretation:

- Steps 1-5 are cheap: 52 total rows.
- Steps 6-8 dominate: 533 of 585 rows.
- Any speedup has to reduce the late-step tail without losing phase closure.

## Evidence Map

Primary compact reports and data:

| topic | path |
| --- | --- |
| Baseline replication recipe | `docs/tsn65-baseline-replication.md` |
| Roadmap/checkpoint status | `docs/nonlinear_solver_roadmap.md` |
| Original stabilization handout | `docs/nonlinear-solver-stabilization-handoff.md` |
| Baseline local report | `results/tsn65-full-baseline-replication-20260520-123252/report.md` |
| Zip/input inspection | `results/tsn65-zip-runs-20260519-190000/report.md` |
| Tangent attribution | `results/nonlinear-solver-tangent-attribution-20260519-113815/report.md` |
| Material tangent audit | `results/nonlinear-solver-material-tangent-audit-20260519-143422/report.md` |
| CSL analytical tangent experiment | `results/tsn65-csl-analytic-tangent-20260519-230004/report.md` |
| CP1 globalization sweep | `results/tsn65-globalization-sweep-20260520-cp1/report.md` |
| CP2 adaptive cutback | `results/tsn65-adaptive-cutback-20260520-cp2/report.md` |
| CP3 indirect control | `results/tsn65-indirect-control-20260520-cp3/report.md` |
| CP4 arc-length prototype | `results/arc-length-cp4-20260521-cp4/report.md` |
| CP5 TS-N65 arc-length | `results/tsn65-arc-length-cp5-20260521/report.md` |
| CP5 radius calibration | `results/tsn65-arc-length-cp5-radcal-20260521/report.md` |
| CP5b gauge arc-length | `results/tsn65-gauge-arc-cp5b-20260521-094117/report.md` |
| CP5c first-iteration screen | `results/tsn65-initial-guess-screen-20260521/report.md` |
| CP5c guarded predictor screen | `results/tsn65-initial-guess-predictor-20260521/report.md` |
| Latest baseline gate | `results/tsn65-cp5b-baseline-gate-20260521/report.md` |

## Experiment Ledger

### 1. TS-N65 zip/input inspection and local solver substitution

Why it was tested:

- The user supplied benchmark archives. We needed to confirm the TS-N65 deck
  was valid and runnable locally.

What was tested:

- Extracted `TS-N_65.zip` and `Dogbone.zip`.
- Attempted the original TS-N65 deck.
- Replaced unavailable `PardisoLDLT` with local `DeflatedFGMRES`/HYPRE setup.

Result:

- The original TS-N65 solver deck requests `PardisoLDLT`, which is unavailable
  in this build:

```text
Solver type PardisoLDLT is not implemented
```

- The benchmark deck itself loads successfully with approximately:

```text
DoF: 567923
constraint DoF: 19254
directly prescribed DoF: 13
elements: 514253
aux nodes: 18513108
```

Failure/useful finding:

- Original solver settings cannot be used locally as-is.
- The replacement baseline must be defined by the full `DeflatedFGMRES`/HYPRE
  deck, not by a shortened nonlinear solver block.

### 2. Strict baseline replication

Why it was tested:

- Before any solver changes, we needed a reproducible local target.

What was tested:

- Full TS-N65 deck.
- `DeflatedFGMRES` plus HYPRE settings from the saved baseline.
- 16 threads.
- `limit_tolerance 0`.
- Fixed 8 quarter steps.

Result:

- Successful replication on current checkpoint code:

```text
8/8 steps
6,6,10,13,17,183,187,163
585 rows
```

- Latest baseline gate after the CP5b arc-length code still matches exactly:

```text
G1-baseline: 8 steps, 585 rows, tail 6-8 = 533, duration 00:22:47.404
```

Failure/useful finding:

- This is the hard gate. Any proposed speedup must complete this same target or
  clearly state that it is tracing a different path.
- Runtime varies by machine, so row sequence is the primary criterion.

### 3. Reduced solver deck failure

Why it was tested:

- Early replication attempts used a shorter `solver.inp` block. We needed to
  know whether the strict baseline depended on omitted DFGMRES/HYPRE settings.

What was tested:

- Reduced nonlinear block with `DeflatedFGMRES`, but missing multiple DFGMRES
  and HYPRE controls.

Result:

- The run no longer matched baseline:

```text
early sequence changed to 8,9,12,14,18
step 6 emitted DFGMRES true-tolerance warnings
NaNs began in step 6
run stopped at max_iterations = 300
```

Failure/useful finding:

- The full linear/preconditioner block is necessary for reproduction.
- This does not prove the linear solver is the primary nonlinear bottleneck,
  but it proves that preconditioner configuration changes nonlinear behavior.

### 4. Tangent attribution at TS-N65 hard state

Why it was tested:

- Step 6 is the first very slow baseline step. We needed to determine whether
  the assembled tangent is locally consistent near that state.

What was tested:

- Global nonlinear tangent finite-difference check at step 6, iteration 10.
- Element-level top-contributor attribution.
- Matrix type `tangent`, then diagnostic matrix type `consistent`.

Key results:

Global Newton-direction check with production `tangent`:

| direction | relerr | cosine | `||Kp||` | `||fd||` |
| --- | ---: | ---: | ---: | ---: |
| Newton | 5.320964 | 0.707442 | 1.186577e7 | 1.983825e6 |
| Random | 0.0329979 | 0.999465 | 2.92288e9 | 2.91168e9 |

Top element rows:

```text
element_name = LDPMTetra
material_id = 0
material_name = CSL material
material_status_names = CSL mat. statusx12
```

Top-50 local summary:

| matrix | mean relerr | max relerr | mean cosine | max mismatch |
| --- | ---: | ---: | ---: | ---: |
| production `tangent` | 3.75268 | 7.00578 | -0.308147 | 358831 |
| numerical `consistent` diagnostic | 0.133558 | 0.466171 | 0.986366 | 10748.7 |

Failure/useful finding:

- The Newton direction is the problematic direction, not every direction.
- The top hard-state local mismatch is in `LDPMTetra` elements backed by CSL
  material/status objects.
- A numerical `consistent` diagnostic greatly improves the element-local
  tangent mismatch, but that alone does not prove a production solver speedup.

### 5. Material tangent audit

Why it was tested:

- The hard-state attribution pointed to CSL material statuses. We needed
  material-level one-point finite-difference checks.

What was tested:

- Local `OAS_material_tangent_audit`.
- CSL and LDPM material probes.
- Cases:
  - `elastic_loading`
  - `damage_growth`
  - `damage_growth_frozen`
  - `damaged_unloading`
  - `damaged_unloading_frozen`
- Branches:
  - `elastic`
  - `secant`
  - `tangent`
  - `archived_csl_damage_tangent`
  - `consistent`

Key CSL rows:

| case | evaluation | branch | relerr | cosine | interpretation |
| --- | --- | --- | ---: | ---: | --- |
| elastic_loading | actual | tangent | 2.08e-13 | 1.0 | closes |
| damage_growth | actual | secant | 3.4089 | -0.9963 | bad active-damage derivative |
| damage_growth | actual | tangent | 3.4089 | -0.9963 | same mismatch in legacy tangent |
| damage_growth | actual | archived_csl_damage_tangent | 2.38e-5 | 1.0 | local fix closes |
| damage_growth | actual | consistent | 3.01e-5 | 1.0 | numerical reference closes |
| damage_growth_frozen | frozen | tangent | 6.12e-13 | 1.0 | frozen path closes |
| damaged_unloading | actual | tangent | 7.77e-13 | 1.0 | unloading closes |

Key LDPM rows:

| case | branch | relerr | cosine | interpretation |
| --- | --- | ---: | ---: | --- |
| elastic_loading | secant | 186.796 | 0.99996 | diagnostic mismatch |
| elastic_loading | consistent | 2.32e-6 | 1.0 | numerical reference closes |
| damage_growth | secant | 8.76e-12 | 1.0 | closes in this local path |
| damaged_unloading | secant | 9.66e-12 | 1.0 | closes |

Failure/useful finding:

- The CSL active-damage tangent defect is real.
- The defect is branch-specific: elastic, frozen, and damaged unloading are
  not the failing CSL cases.
- The archived analytical branch and numerical consistent branch close the
  local CSL active damage-growth probe.
- The LDPM local elastic mismatch remains a diagnostic item, but the TS-N65
  hard-state top rows were attributed to CSL statuses, so CSL is the primary
  local tangent suspect for this benchmark.

### 6. Analytical CSL active-damage tangent experiment

Why it was tested:

- Since the material audit found a real CSL active-damage tangent defect, we
  tested whether an analytical derivative improves TS-N65 phase closure.

What was tested:

- Analytical CSL active-damage tangent:

```text
sigma = (1 - damage(eps)) * D_elastic * eps
D_tangent = (1 - damage) * D_elastic - (D_elastic * eps) outer d_damage/d_eps
```

- Temporarily wired into production `tangent`.
- Later archived as `stiff_matrix_type archived_csl_damage_tangent`.

Results:

1. Mid-step rebuild variant:

```text
stiff_matrix_type tangent
stiffness_matrix_iter_update 10
```

- Did not reach step 6.
- Became unstable in step 3 immediately after the iteration-10 matrix rebuild.
- Last useful row:

```text
step 3 iteration 15 residual 9.897274e-01 displacement 9.449099e-01 energy 1.552062e+00
```

2. Step-start matrix only:

```text
stiffness_matrix_iter_update -1
```

- Steps 1-5 converged with rows `8,9,12,14,18`.
- Step 6 failed at `max_iterations = 300`.
- Total duration before failure: `00:10:56.446`.

3. Step-start analytical tangent plus actual-state backtracking:

- Steps 1-5 converged with rows `8,9,12,14,18`.
- Step 6 converged in `54` rows with repeated `alpha=0.5`.
- Step 7 stalled and was manually stopped near:

```text
iteration 75 residual 5.843750e-03 displacement 3.003728e-04 energy 5.399725e-03
```

Failure/useful finding:

- The analytical CSL tangent fixes the local material derivative and strongly
  improves local element diagnostics.
- It is too aggressive or indefinite for the current strict load-control Newton
  path when rebuilt mid-step.
- With line search, it improves step 6 but still stalls in step 7.
- Therefore the current blocker is no longer simply "CSL tangent is locally
  wrong." The phase closure needs a stabilized nonlinear/path-following method
  or model-level regularization.

### 7. Production `consistent` matrix smoke

Why it was tested:

- The numerical `consistent` diagnostic reduced hard-state local mismatch, so
  we checked whether production `stiff_matrix_type consistent` is viable.

What was tested:

```text
stiff_matrix_type consistent
total_time 1.25e-3
```

Result:

- The model loaded and assembled.
- The first linear solve was impractical with current DFGMRES/HYPRE:

```text
DeflatedFGMRES warning: performed 500 iterations and reached true relative residual 0.843808, required true tolerance is 0.1
```

- The process was stopped after the first nonlinear row. It was already using
  about 27 GB RSS.

Failure/useful finding:

- A locally better matrix can be globally unsuitable for the current linear
  solver/preconditioner setup.
- Any future true consistent-tangent production attempt needs linear solver and
  preconditioner work, or tangent regularization/blending.

### 8. CP1: Existing globalization sweep on legacy tangent

Why it was tested:

- Before adding more algorithms, test implemented globalization controls against
  the strict fixed-step baseline.

What was tested:

| variant | result |
| --- | --- |
| `G1-baseline` | complete baseline reproduction |
| `G1-backtracking-frozen-actual` | partial/fails |
| `G1-backtracking-actual` | partial/manual stop |

The queued but not run variants were:

```text
G1-bisection-frozen-actual
G1-fixed-damping-05
G1-adaptive-damping
G1-trust-stepnorm
G1-line-search-adaptive-K
```

They were left queued because the first two backtracking variants already
showed the main failure pattern: step 6 improves but step 7 gets worse.

Results:

| variant | exit | steps | rows | tail 6-8 | min alpha | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| G1-baseline | 0 | 8 | 585 | 533 |  | neutral |
| G1-backtracking-frozen-actual | 124 | 8 | 495 before timeout | 443 before timeout | 0.0625 | failed |
| G1-backtracking-actual | manual stop | 7 | 161 partial | 109 partial | 0.125 | partial worse |

Detailed observations:

- `G1-backtracking-frozen-actual` reduced step 6 from `183` rows to `63`, but
  step 7 grew from `187` rows to `270` and took `00:50:36.666`; step 8 was
  still unconverged when the 90-minute runner timeout fired.
- `G1-backtracking-actual` also reduced step 6 to `63` rows but showed the same
  step-7 pathology. It was manually stopped at step 7 row `45` with:

```text
residual 6.728244e-03
displacement 2.237895e-03
energy 1.310600e-02
```

Failure/useful finding:

- Backtracking can improve step 6.
- The same strategy destabilizes or greatly slows step 7.
- This strongly suggests a nonmonotone transition or path/state issue rather
  than a simple need for smaller Newton increments everywhere.

### 9. CP2: Adaptive cutback/stagnation

Why it was tested:

- Avoid burning hundreds of iterations on stalled increments.
- Disable fallback acceptance and test whether smaller adaptive steps are
  enough to reach `total_time = 1e-2`.

Common adaptive settings:

```text
max_iterations 40
min_time_step 1.25e-5
max_time_step 1.25e-3
critical_step_decrease 0.5
step_decrease 0.5
step_increase 1.2
enlargeIt 6
shortenIt 15
limit_tolerance 0
nonlinear_line_search_cutback_on_fail 1
nonlinear_stagnation_cutback 1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25
```

Results:

| variant | exit | accepted | rows | furthest accepted time | min dt | cutbacks | fallback |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G2-stagnation-only | 1 | 10 | 183 | 0.007421875 | 1.25e-5 | 5 | 0 |
| G2-backtracking-actual | 1 | 8 | 168 | 0.0078125 | 1.25e-5 | 5 | 0 |

Last/failure characteristics:

- `G2-stagnation-only` failed after reaching `min_time_step`.
- `G2-backtracking-actual` reached further physical time but still failed after
  repeated line-search and stagnation cutbacks.
- No NaNs and no fallback acceptance occurred.

Failure/useful finding:

- Adaptive cutback is useful as a safety mechanism.
- It does not by itself solve the late TS-N65 load-control path.
- It reduces local row burn but can get trapped at the minimum time step before
  reaching the final target.

### 10. CP3: Existing indirect displacement control

Why it was tested:

- Before implementing arc-length, use the existing IDC support as a cheaper
  path-following/control experiment.

What was tested:

- Historical `src/benchmark/Indirect_Control`.
- TS-N65 coordinate-based displacement gauge using intended `v01 uz` measure:

```text
indirect_control 2
ic_xcoords 2.1 2.1
ic_ycoords 0 0
ic_zcoords 0 -0.4
ic_directions 2 2
ic_displ_weights -1 1
ic_force_weights 0 0
```

- Sign calibration with flipped weights.
- Raw final target and calibrated final target.

Important implementation note:

- OAS function indices are zero-based.
- The runner appends a local IDC function and uses the actual appended index.

Results:

| variant | exit | steps | rows | v01 result | verdict |
| --- | ---: | ---: | ---: | ---: | --- |
| IDC benchmark | 1 | 146 | 21910 | n/a | failed legacy benchmark |
| TS-N65 minus calibration | 0 | 1 | 25 | 1.575315e-6 | direction selected |
| TS-N65 plus calibration | 1 | 1 | 0 | n/a | wrong sign failed immediately |
| TS-N65 raw IDC strict | 0 | 8 | 193 | 1.233328e-5 | completed but undershot |
| TS-N65 calibrated final v01 | -15 | 7 | 218 | 1.144418e-5 | stalled/manual stop |

Strict target:

```text
v01 target = 1.519697e-5
raw IDC final v01 = 1.233328e-5
relative miss = 18.84%
```

Calibrated run stalled at step 7 iteration 69:

```text
residual 2.352053e-03
displacement 4.419994e-04
energy 4.186416e-03
```

Failure/useful finding:

- IDC can follow a lower-amplitude path cheaply, but that is not the strict
  physical target.
- Scaling the IDC target to recover the final `v01` reintroduces late-step
  nonlinear instability.
- The internal IDC coordinate gauge is not equivalent to the exported `LD.out`
  interpolation gauge.

### 11. CP4: Arc-length prototype on a small benchmark

Why it was tested:

- Implement a disabled-by-default arc-length prototype and verify that the
  infrastructure works before applying it to TS-N65.

What was implemented:

```text
nonlinear_control load|indirect|arc_length
arc_length_radius_initial
arc_length_radius_min
arc_length_radius_max
arc_length_psi
arc_length_shrink
arc_length_expand
arc_length_target_iterations
arc_length_max_iterations
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
```

V1 arc-length method:

- proportional external loading only;
- reference load `f_ext(load=1) - f_ext(load=0)`;
- Schur-complement corrector;
- optional line-search/snapshot rollback path.

Small benchmark result:

| variant | verdict | steps | rows | final lambda | notes |
| --- | --- | ---: | ---: | ---: | --- |
| legacy load default | passed | 3 | 6 | 0.3 | no keyword regression |
| explicit load control | passed | 3 | 6 | 0.3 | explicit control regression |
| arc-length prototype | passed | 3 | 6 | 0.3 | spherical arc path |
| arc-length line search | passed | 4 | 18 | 0.3 | alpha min 0.5 |

Failure/useful finding:

- Arc-length infrastructure works on a small proportional-load benchmark.
- Line-search/snapshot rejection is functional.
- Strict TS-N65 baseline gate after CP4 still matched exactly, so the disabled
  controls preserved legacy behavior.

### 12. CP5: TS-N65 arc-length application

Why it was tested:

- Apply the CP4 arc-length prototype to the actual target.

What was tested:

- Proportional-load arc-length on TS-N65.
- Finite-difference arc-length reference for prescribed-displacement driven
  TS-N65.
- Radius calibration.

Results:

| variant | result |
| --- | --- |
| CP5 proportional smoke | unsupported zero reference |
| finite-difference radius `1.25e-3` | failed in first arc step after 3 rows |
| finite-difference radius `1.35e-2` | partial/manual stop, 39 rows before first target was stably reached |

Important details:

- Direct proportional-load arc-length is not appropriate for TS-N65 because the
  deck is prescribed-displacement driven and has no nonzero reduced nodal
  reference load `f_ext(1)-f_ext(0)`.
- The finite-difference reference mode can derive a direction from prescribed
  BC/load parameter perturbation, but the spherical norm is badly scaled for
  this model.
- Radius `1.35e-2` calibrated the first predictor to `lambda = 0.00125`, but
  the corrector oscillated and was already slower than strict baseline step 1:

```text
arc run: 39 rows before first strict quarter-step was stably reached
baseline step 1: 6 rows
```

Failure/useful finding:

- Generic spherical arc-length is not competitive for this prescribed-BC
  TS-N65 baseline.
- A useful TS-N65 arc method likely needs a physically meaningful displacement
  constraint or better scaling, not a full-vector spherical norm.

### 13. CP5b: Gauge-constrained arc-length

Why it was tested:

- CP5 suggested that spherical arc-length is badly scaled. CP5b replaced the
  spherical constraint with a displacement gauge based on IDC.

What was implemented:

```text
arc_length_constraint gauge
arc_length_sign_strategy monotone_lambda
arc_length_auto_radius
arc_length_gauge_tolerance
```

Gauge mode:

- Arc-length load factor remains the equilibrium load parameter.
- Corrector constraint is a displacement-only IDC gauge.
- Force-weighted IDC gauges are rejected because the reaction-force gauge is not
  consistently linearized yet.

Tested variants:

| variant | exit | steps | rows | final v01 | target v01 | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| first-quarter uncalibrated | 0 | 1 | 25 | 2.167438e-6 | 2.613423e-6 | wrong gauge |
| first-quarter calibrated | 0 | 1 | 25 | 2.613423e-6 | 2.613423e-6 | target hit but slower |
| first-quarter calibrated line search | 0 | 1 | 25 | 2.613423e-6 | 2.613423e-6 | slower wall time |
| full calibrated | 130 | 7 | 235 | 1.14461e-5 | 1.519697e-5 | failed/manual stop |

Detailed results:

- Uncalibrated gauge mode hit the internal IDC target exactly but exported
  `LD.out` `v01` was `17.07%` low.
- Calibrated first-quarter mode matched the exported strict baseline first-step
  `v01 = 2.613423e-6`.
- Calibrated first-quarter mode took `25` rows and about `4:08`; strict
  baseline first step takes `6` rows and about `0:35`.
- Full calibrated mode reached reported step 7 and `235` rows before being
  stopped because step 7 drifted away from convergence. Last parsed row:

```text
step 7 iteration 89
residual 3.963647e-3
displacement 2.671126e-3
energy 4.556947e-3
```

Failure/useful finding:

- Gauge arc-length is mechanically usable enough to hit a prescribed
  displacement gauge.
- It is not a speed replacement for the strict TS-N65 baseline.
- The current gauge is not the same as the exported interpolation gauge used in
  `LD.out`, so calibration is needed and path comparability is imperfect.
- The full run followed a softer linear gauge path to the final target rather
  than the strict baseline loading history, so it is a continuation experiment,
  not a baseline replacement.

### 14. CP5c: Initial guess and first-iteration predictor screen

Why it was tested:

- The direct load-control solver starts each step from the last committed
  solution plus new prescribed BCs, with no free-DOF nonlinear predictor.
- The late TS-N65 cost could have been partly caused by starting each new load
  step too far from the equilibrium branch.

What was implemented/tested:

- Existing no-code first-correction direction screen:

```text
first_iteration_stiff_matrix_type elastic
first_iteration_stiff_matrix_type secant
first_iteration_stiff_matrix_type archived_csl_damage_tangent
```

- New disabled-by-default guarded load-step predictors:

```text
nonlinear_initial_guess off|last_step|two_step
nonlinear_initial_guess_alpha
nonlinear_initial_guess_start_step
nonlinear_initial_guess_max_norm_ratio
nonlinear_initial_guess_guard
nonlinear_initial_guess_guard_merit residual|energy|mixed
nonlinear_initial_guess_accept_ratio
nonlinear_initial_guess_frozen_eval
```

- `last_step` predicts the new step-start field from the last accepted reduced
  free-DOF increment:

```text
u_n^0 = u_{n-1} + alpha * (dt_n / dt_{n-1}) * delta_u_{n-1}
```

- The predictor runs only in direct load control, after prescribed Dirichlet
  BCs are updated and before the first step-start force evaluation.
- The guard evaluates residual merit at the no-predictor and predicted states
  and restores state when rejected.

No-code first-iteration results:

| variant | result |
| --- | --- |
| baseline | completed `8/8`, `585` rows |
| firstit elastic | worsened early rows and was stopped in step 6 |
| firstit secant | matched through step 6, then worsened in step 7 |
| firstit archived CSL | matched through step 5, then worsened in step 6 |

Guarded predictor results:

| variant | status | predictor accepted | useful finding |
| --- | ---: | ---: | --- |
| baseline | completed `8/8`, `585` rows | 0 | no-regression gate passed |
| alpha `0.25` | stopped in step 6 | 5/5 | tracked baseline, slightly worse by row 69 |
| alpha `0.50` | stopped in step 6 | 5/5 | lower first row, slightly worse by rows 30-40 |
| alpha `0.75` | stopped in step 6 | 5/5 | much lower first-row energy, slightly worse by rows 30-40 |
| alpha `1.00` | stopped in step 6 | 5/5 | strong early predictor effect, worse by row 69 |

Step-6 comparison highlights:

| variant | row 0 residual | row 0 energy | later comparison |
| --- | ---: | ---: | --- |
| baseline | `2.548287e-02` | `1.001264e-01` | row 69 residual/energy `3.831386e-03` / `7.130556e-03` |
| alpha `0.50` | `2.407001e-02` | `6.225239e-02` | row 39 slightly worse than baseline |
| alpha `0.75` | `2.372770e-02` | `2.848583e-02` | row 39 slightly worse than baseline |
| alpha `1.00` | `2.547188e-02` | `5.114121e-02` | row 69 slightly worse than baseline |

Failure/useful finding:

- Step-start predictors are accepted and can reduce the first residual/energy
  row, so the implementation is functional.
- The predictor benefit disappears quickly in the nonlinear iteration. By
  rows 30-40, or by row 69 for the longer screens, the variants match or
  slightly worsen the strict baseline tail.
- Therefore the expensive TS-N65 tail is not primarily caused by the lack of a
  free-DOF step-start displacement predictor. The bottleneck remains later
  path/tangent/history behavior inside the nonlinear step.

## Cross-Cutting Failure Interpretation

### A. The late-step tail is the real speed target

The first five strict baseline steps consume only 52 rows. Steps 6-8 consume
533 rows. Any method that changes early behavior but does not reduce the
late-step tail is not useful for the stated speed goal.

### B. Local tangent correctness is necessary but not sufficient

The CSL active-damage tangent mismatch is real:

```text
legacy tangent damage-growth relerr = 3.4089, cosine = -0.9963
archived analytical tangent relerr = 2.38e-5, cosine = 1.0
```

But strict TS-N65 phase closure did not improve:

- mid-step analytical tangent destabilized step 3;
- step-start analytical tangent failed in step 6;
- line-search analytical tangent closed step 6 but stalled in step 7.

Likely implication:

- The true softening tangent may be indefinite/aggressive for load control.
- A stabilized tangent, trust region, tangent limiter/blend, or model
  regularization may be needed.

### C. Line search is a partial improvement, not a solution yet

Line search/backtracking repeatedly produced the most interesting partial
speedup:

```text
step 6: 183 rows baseline -> about 63 rows with backtracking
```

But it then worsened step 7:

```text
step 7: baseline 187 rows -> 270 rows/timeout or manual stop
```

Likely implication:

- The solver can find better local increments in step 6.
- Step 7 may cross a material/history/localization transition where monotone
  residual decrease is too restrictive, the accepted path changes history, or
  the stiffness update strategy becomes harmful.

### D. Adaptive cutback is useful as safety, not speed

Adaptive cutback prevents unlimited row burn and fallback acceptance, but in
TS-N65 it hit `dtmin` before final time.

Likely implication:

- Smaller load steps alone are not enough, or `dtmin`/control path needs a more
  physical continuation variable.

### E. Existing IDC/gauge control does not match the exported physical gauge

The internal coordinate-based IDC gauge and exported `LD.out` `v01` gauge differ
significantly. In CP5b, the uncalibrated first-quarter run hit the internal
gauge but was `17.07%` low in exported `v01`.

Likely implication:

- A serious displacement-control experiment should constrain the exact same
  interpolated quantity used for reporting, or the comparison will remain
  ambiguous.

### F. Generic spherical arc-length is poorly scaled for TS-N65

The proportional-load arc-length path is inapplicable because TS-N65 is driven
by prescribed displacement/BC changes. The finite-difference reference path can
move the model but is too expensive and oscillatory.

Likely implication:

- If arc-length is pursued, it should be a problem-specific continuation method
  with a meaningful control measure and scaling, not the current full-field
  spherical norm.

### G. Production `consistent` tangent currently has a linear-solver barrier

The numerical `consistent` matrix improves local diagnostics but was not viable
with current DFGMRES/HYPRE settings.

Likely implication:

- Either tune/change the linear solver for consistent/indefinite matrices, or
  use a regularized/blended tangent rather than a raw consistent tangent.

## Useful Code Infrastructure Already Added

Diagnostics and parsing:

- `scripts/summarize_oas_run.py`
- `scripts/run_tsn65_globalization_sweep.py`
- `scripts/run_tsn65_indirect_control.py`
- `scripts/run_arc_length_checkpoint.py`
- `scripts/run_tsn65_arc_length_cp5.py`
- `scripts/run_tsn65_gauge_arc_cp5b.py`

Solver controls/infrastructure:

- nonlinear line search, damping, stagnation/cutback hooks;
- nonlinear material snapshot/rollback;
- material status clone/restore/hash support;
- guarded nonlinear initial-guess predictor infrastructure;
- tangent attribution diagnostics;
- `archived_csl_damage_tangent`;
- numerical `consistent` diagnostic branches;
- disabled-by-default arc-length prototype;
- gauge-constrained arc-length experimental mode.

Important safeguard:

- Latest strict baseline gate after these code changes still matches
  `6,6,10,13,17,183,187,163` and `585` rows.

## Questions For Expert Review

1. Is the strict `585`-row load-control baseline a physically meaningful path,
   or is it a numerically difficult but acceptable load-controlled traversal of
   an unstable/softening region?

2. Given that the locally correct CSL active-damage tangent destabilizes or
   stalls the strict phase run, should the production tangent be:

```text
legacy degraded stiffness,
full consistent tangent,
positive-definite/secant blend,
projected/limited tangent,
LM-regularized tangent,
or something material-model specific?
```

3. Is the step 6 to step 7 failure likely caused by:

```text
snapback/load-control instability,
damage localization,
history path pollution from trial evaluations,
linear preconditioner degradation,
or tangent indefiniteness?
```

4. Should the next continuation attempt control the exact exported `LD.out`
   interpolation gauge instead of the nearest-node IDC coordinate gauge?

5. Is it preferable to improve global nonlinear globalization first, or move to
   model-level stabilization such as local damage substepping, viscosity, or
   crack-band/nonlocal regularization?

6. If a raw consistent tangent is desired, what linear solver/preconditioner
   strategy should replace the current DFGMRES/HYPRE path?

## Suggested Next Experiments

These are not yet proven. They are proposed because they directly address the
failure modes above.

## CP6 Failure-Mode Isolation Implementation

Checkpoint implementation status:

- Added disabled-by-default material-state dumps for accepted nonlinear steps:
  `nonlinear_state_dump`, `nonlinear_state_dump_steps`,
  `nonlinear_state_dump_top_damage`, `nonlinear_state_dump_include_coordinates`,
  and `nonlinear_state_dump_directory`.
- Added `scripts/compare_tsn65_states.py` to compare two dumped step summaries
  and top-damage overlap.
- Added `stiff_matrix_type csl_stabilized_tangent` with bounded CSL
  active-damage correction controls: `csl_tangent_beta`,
  `csl_tangent_softening_limit`, `csl_tangent_active_only`, and
  `csl_tangent_log_stats`.
- Added true LM-style regularized Newton controls:
  `nonlinear_lm_regularization`, `nonlinear_lm_mu_initial`,
  `nonlinear_lm_mu_min`, `nonlinear_lm_mu_max`, `nonlinear_lm_mu_growth`,
  `nonlinear_lm_mu_shrink`, `nonlinear_lm_max_trials`,
  `nonlinear_lm_diag`, and `nonlinear_lm_accept`.
- Added CP6 variants to `scripts/run_tsn65_globalization_sweep.py` for
  state-equivalence gates, stabilized CSL tangent screens, and LM screens.

Implementation notes:

- The stabilized CSL tangent uses the legacy degraded stiffness as the base and
  subtracts a limited fraction of the active-damage rank-one correction:

```text
K = K0 - beta * scale * S
scale = min(1, gamma * ||K0|| / max(||S||, tiny))
```

- `beta = 0` recovers legacy degraded stiffness. `beta = 1` with no softening
  limit approaches the archived analytical active-damage tangent.
- LM v1 regularizes the currently assembled effective matrix:

```text
(Keff + mu * D) du = R
```

  with `D = abs(diag(Keff))` by default. `row_sum_diag` is also available.
  `elastic_diag` is parsed for compatibility and currently falls back to
  `abs_diag`.
- LM can be combined with the existing backtracking/bisection line search: the
  regularized direction is solved first, then the existing trial-evaluation and
  rollback machinery chooses an accepted alpha.

Initial smoke check:

- Build passed with `/tmp/oas_tsn65_full_baseline_build`.
- A small CSL benchmark accepted the new solver controls and wrote
  `state/step_001_accepted_summary.json`, `damage_hist.csv`,
  `top_damage.csv`, and `status_hashes.csv`.
- `scripts/compare_tsn65_states.py` successfully compared the smoke dump to
  itself with matching global material-state hash and top-damage overlap.

First TS-N65 CP6 gate:

- `D0-baseline-state-dump` completed the strict baseline unchanged:
  `6,6,10,13,17,183,187,163`, total `585` rows, no fallback acceptance,
  no warnings, and no NaNs.
- Step 6 accepted dump:
  `6171036` CSL statuses, `125947` active damage-growth statuses,
  `148183` damaged statuses, `damage_p99 = 0.48690121350286242`,
  `damage_increment_p99 = 0.1270366795664456`, global state hash
  `6622989350471731203`.
- Step 7 accepted dump:
  `6171036` CSL statuses, `135027` active damage-growth statuses,
  `184849` damaged statuses, `damage_p99 = 0.58104234307705516`,
  `damage_increment_p99 = 0.10831286900670489`, global state hash
  `14127883338335681505`.
- Compact report:
  `results/tsn65-cp6-state-dump-20260521/report.md`.

Completed TS-N65 CP6 queue:

```text
D1-backtracking-state-dump
D2-archived-csl-backtracking-state-dump
S1-cslbeta005-gamma005-stepstart
S2-cslbeta010-gamma005-stepstart
S3-cslbeta020-gamma005-stepstart
S4-cslbeta010-gamma010-stepstart
S5-cslbeta020-gamma010-stepstart
LM1-legacy-mu1e-4
LM2-legacy-mu1e-3
LM3-cslbeta005-gamma005-mu1e-4
LM4-cslbeta010-gamma005-mu1e-4
LM5-cslbeta010-gamma010-mu1e-4
```

Queue report:
`results/tsn65-cp6-queue-20260522-0652/report.md`.

Aggregate result table:

| variant | exit | steps | rows | tail 6-8 | duration | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| D1-backtracking-state-dump | 0 | 8 | 576 | 524 | 01:16:42.602 | not promoted |
| D2-archived-csl-backtracking-state-dump | -15 | 7 | 211 | 159 |  | failed/manual stop |
| S1-cslbeta005-gamma005-stepstart | 1 | 6 | 352 | 300 | 00:07:27.422 | failed |
| S2-cslbeta010-gamma005-stepstart | 1 | 6 | 352 | 300 | 00:07:27.738 | failed |
| S3-cslbeta020-gamma005-stepstart | 1 | 6 | 352 | 300 | 00:07:27.614 | failed |
| S4-cslbeta010-gamma010-stepstart | 1 | 6 | 352 | 300 | 00:07:27.689 | failed |
| S5-cslbeta020-gamma010-stepstart | 1 | 6 | 352 | 300 | 00:07:28.775 | failed |
| LM1-legacy-mu1e-4 | 1 | 6 | 72 | 25 | 00:19:57.716 | failed |
| LM2-legacy-mu1e-3 | 1 | 6 | 78 | 26 | 00:21:42.972 | failed |
| LM3-cslbeta005-gamma005-mu1e-4 | 1 | 6 | 72 | 25 | 00:20:14.966 | failed |
| LM4-cslbeta010-gamma005-mu1e-4 | 1 | 6 | 72 | 25 | 00:20:15.910 | failed |
| LM5-cslbeta010-gamma010-mu1e-4 | 1 | 6 | 72 | 25 | 00:20:17.283 | failed |
| LM6-legacy-relaxed102-mu1e-4 | 1 | 6 | 74 | 27 | 00:20:29.124 | failed |
| LM8-legacy-relaxed105-mu1e-4 | 1 | 6 | 74 | 27 | 00:20:28.532 | failed |
| LM9-legacy-relaxed110-mu1e-4 | 1 | 6 | 80 | 33 | 00:21:58.819 | failed |

CP6 queue conclusions:

- The apparent row-count improvement of `D1` is not a strict speedup. It
  completes `8/8` with `576` rows, but state comparison proves step 6 is a
  different material path: the global status hash differs, active CSL
  damage-growth statuses increase by about `17%`, p99 damage increment rises by
  about `23%`, and top-1000 damage overlap is only `264/1000`. Step 7 then
  worsens from `187` rows to `270`.
- `D2` is the same failure pattern with a more aggressive tangent branch:
  step 6 drops to `47` rows but the material state shifts similarly and step 7
  stalls; it was stopped manually.
- The step-start stabilized CSL tangent variants `S1`-`S5` all fail at step 6
  after `300` rows. Logs show the active-damage correction is inactive at the
  step-start rebuild (`active_damage_status_count = 0`), so this exact tangent
  screen does not apply the intended correction where it matters.
- The LM variants `LM1`-`LM5` all fail at step 6. They reach a near-residual
  solution quickly, but the energy error plateaus around `1.8e-3` and the next
  LM trial is rejected for all tested `mu` retries. The implementation is also
  expensive because each trial refactorizes and loses the baseline DFGMRES
  deflation basis (`basis_size=0` in setup logs).
- Relaxed-merit LM gates did not solve the plateau. A `2%` and `5%` allowed
  merit increase both still failed at step 6 iteration `27`; the `10%` gate
  pushed the failure to iteration `33`, but only by accepting controlled error
  growth from merit about `1.89` to `2.65`. The rejected row then exceeded the
  relaxed gate again. This is useful negative evidence: widening the acceptance
  criterion delays failure but does not create convergence.

CP6 verdict:

- Do not promote any queued variant.
- Do not continue broad beta/gamma/mu sweeps until the energy-merit plateau and
  repeated refactorization/deflation loss are addressed.
- The row-count speedups seen in line search are path-shifting, not strict-path
  acceleration.

### 1. Exact exported-gauge control

Implement a control variable that matches the reported `LD.out` `v01` gauge,
not the current internal IDC coordinate approximation.

Test:

- strict baseline first-quarter and final `v01` checkpoints;
- PWL target matching the strict baseline `v01` at all 8 quarter points;
- gauge control plus existing line search and stagnation cutback.

Pass condition:

- same final physical gauge and comparable load/displacement curve;
- fewer late-step rows than the strict baseline.

Implementation update:

- Added optional IDC coordinate interpolation with
  `ic_coordinate_interpolation 1`. When enabled, coordinate-based IDC points
  use the same owning-element/natural-coordinate interpolation pattern as the
  `DisplacementGauge` exporter, falling back to the historical nearest-node
  value only when an owning element is not found.
- Added exact exported-gauge variants to `scripts/run_tsn65_gauge_arc_cp5b.py`:
  `CP5c-exact-gauge-first-quarter`,
  `CP5c-exact-gauge-first-quarter-line-search`,
  `CP5c-exact-gauge-pwl-8step`, and
  `CP5c-exact-gauge-pwl-8step-line-search`.
- First-quarter smoke:
  `results/tsn65-exact-gauge-cp5c-20260522-123843/report.md`.
  The run hit the exported `LD.out` `v01` target exactly without the old
  calibration factor:

```text
target v01 = 2.613423e-6
exported v01 = 2.613423e-6
rows = 29
duration = 00:04:32.739
```

Interpretation:

- The exact interpolation control is now mechanically correct for the first
  quarter-step.
- It is still slower than the strict first baseline step (`6` rows), so it is
  not a speed candidate by itself.
- The eight-point PWL exact-gauge run without line search failed as a speed or
  path-control solution. It matched the first five exported `v01` checkpoints
  exactly, but step 6 developed the same energy-merit plateau and then drifted
  away. It was stopped after `75` step-6 rows:

```text
result root = results/tsn65-exact-gauge-cp5c-pwl-20260522-124701
rows before stop = 229
accepted checkpoint rows = 29,28,30,31,36
step 6 partial rows = 75
step 6 best near-miss = energy merit about 1.12 near row 44
last row = residual 1.678442e-03, displacement 3.711674e-04, energy 3.180927e-03
```

- This is valuable negative evidence: matching the exported displacement gauge
  does not remove the step-6 nonlinear pathology by itself. The line-search
  exact-gauge PWL variant also failed: it completed the first quarter in `29`
  rows, then step 2 immediately hit repeated
  `arc_length_line_search_failed` cutbacks before any accepted row. It was
  stopped after seven identical radius cutbacks.

```text
line-search root = results/tsn65-exact-gauge-cp5c-pwl-ls-20260522-131530
rows before stop = 29
accepted checkpoints = 1
cutbacks = 7
last radius = 7.8125e-6
failure reason = arc_length_line_search_failed at step 2
```

### 1b. Do not prioritize simple last-step predictors

The CP5c predictor screen already tested guarded last-step extrapolation with
alpha values `0.25`, `0.50`, `0.75`, and `1.0`. These predictors changed the
first row but did not reduce the late-step tail. Further predictor work should
only be considered if it uses a substantially different, physics-informed
state variable or exact control-gauge continuation, not a simple repeat of the
last free-DOF increment.

### 2. Tangent stabilization instead of raw analytical tangent

Test a controlled tangent blend/limiter:

```text
K = (1 - beta) * K_legacy + beta * K_archived_csl_damage_tangent
```

or an eigenvalue/diagonal/positive-definite regularization in active softening
zones.

Start with small `beta` values and fixed strict 8-quarter runs.

Pass condition:

- preserve `8/8` strict target;
- reduce rows in steps 6-8.

### 3. Trust-region or LM regularized Newton on the strict path

Backtracking helps step 6 but worsens step 7. A trust-region/LM method may be
better when the tangent direction is unreliable or indefinite.

Test:

- line-search off/on combinations;
- diagonal regularization based on `abs(diag(K))`;
- radius or `mu` logging;
- strict baseline gate after disabling controls.

Pass condition:

- step 6 improvement without step 7 degradation.

### 4. Adaptive matrix rebuild on small alpha/stagnation

The CP1 `G1-line-search-adaptive-K` variant was planned but not run in the
partial sweep. It remains a low-effort direct test.

Test:

```text
nonlinear_line_search backtracking
nonlinear_line_search_evaluation actual or frozen_then_actual
nonlinear_adaptive_matrix_update 1
nonlinear_rebuild_on_small_alpha 0.5
nonlinear_rebuild_on_merit_growth 1
```

Pass condition:

- line search keeps the step 6 benefit;
- adaptive rebuild prevents the step 7 deterioration.

### 5. Model-level stabilization diagnostics

Enter only after solver-level options remain negative, because these change the
physics/model.

Order:

1. material local substepping diagnostics;
2. damage increment limiter for diagnostics only;
3. viscosity/rate-regularization sensitivity;
4. crack-band/nonlocal/gradient regularization investigation.

Pass condition:

- report sensitivity;
- clearly state physics changed;
- show stable load/displacement behavior and improved convergence.

## Minimal Reproduction Commands

Build:

```bash
cmake -S . -B /tmp/oas_tsn65_full_baseline_build \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_VTK=OFF \
  -DUSE_HYPRE=ON
cmake --build /tmp/oas_tsn65_full_baseline_build --target OAS -j 16
```

Strict baseline gate:

```bash
scripts/run_tsn65_globalization_sweep.py \
  --variants G1-baseline \
  --root results/tsn65-cp5b-baseline-gate-20260521 \
  --jobs 16 \
  --timeout-seconds 7200
```

Summarize a solver log:

```bash
scripts/summarize_oas_run.py path/to/solver.out --baseline tsn65
```

Gauge arc-length CP5b screen:

```bash
scripts/run_tsn65_gauge_arc_cp5b.py \
  --root results/tsn65-gauge-arc-cp5b-20260521-094117 \
  --variants CP5b-gauge-first-quarter-calibrated \
  --timeout 1800 \
  --jobs 16
```

Material tangent audit:

```bash
cmake --build /tmp/oas_deflation_build_material_audit \
  --target OAS_material_tangent_audit OAS -j 16
/tmp/oas_deflation_build_material_audit/bin/OAS_material_tangent_audit \
  --output results/nonlinear-solver-material-tangent-audit-20260519-143422/material_tangent_audit.tsv
```

## Current Bottom Line

The experiments found useful technical facts, but not a speedup:

- The baseline is stable and reproducible.
- The CSL active-damage tangent defect is real, but the locally correct tangent
  is not a production fix by itself.
- Line search gives a real partial improvement in step 6 but transfers the
  difficulty to step 7; CP6 state dumps show this is a material-path shift, not
  a strict-path speedup.
- Last-step nonlinear initial guesses are accepted and improve the first row,
  but do not speed up the late strict TS-N65 tail.
- Step-start stabilized CSL tangent variants do not help because the active
  damage correction is inactive at the step-start matrix rebuild; all tested
  variants fail at step 6.
- LM/regularized Newton variants reach a near-residual solution quickly but
  fail the energy criterion around step 6 with `eneErr ~ 1.8e-3`, and current
  LM refactorization loses the DFGMRES deflation basis on every trial.
- Exact exported-gauge control now matches `LD.out` `v01` without calibration,
  but it is not a speedup. The first-quarter exact-gauge run took `29` rows.
  The eight-point PWL exact-gauge run matched checkpoints through step 5, then
  step 6 hit the same energy plateau and drifted away; adding line search made
  the method fail immediately at step 2 through repeated
  `arc_length_line_search_failed` cutbacks.
- Adaptive cutback, IDC, and current arc-length formulations are useful
  diagnostics, not speed replacements.
- The next expert-guided work should focus on the shared step-6 energy-merit
  plateau and the material/model response behind it. Solver-only branches have
  mostly become diagnostics; model-level regularization or material-integration
  diagnostics are now the likely next branch.
