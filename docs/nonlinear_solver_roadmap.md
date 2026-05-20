# TS-N65 Nonlinear Solver Roadmap

This document is the living checkpoint plan for the TS-N65 nonlinear solver
experiments. The baseline source of truth is
`docs/tsn65-baseline-replication.md`.

## Baseline Gate

Every checkpoint that changes solver behavior must be compared against the
strict TS-N65 full-deck baseline:

- 16 threads.
- Full DFGMRES/HYPRE block.
- `stiff_matrix_type tangent`.
- `min_time_step = max_time_step = time_step = 1.25e-3`.
- `total_time = 1e-2`.
- `limit_tolerance 0`.
- Target step sequence: `6, 6, 10, 13, 17, 183, 187, 163`.
- Target total nonlinear rows: `585`.

Runtime is machine-dependent. The iteration sequence and convergence status are
the primary replication criteria.

## Checkpoint Status

| checkpoint | status | purpose | required comparison |
| --- | --- | --- | --- |
| CP0 | complete | Parser, roadmap, report harness | Existing replicated baseline log |
| CP1 | pending | Existing globalization sweep on legacy `tangent` | Strict fixed-step baseline |
| CP2 | pending | Adaptive cutback/stagnation run | Strict baseline plus adaptive-step metrics |
| CP3 | pending | Existing indirect displacement control before arc-length | IDC benchmark and TS-N65 baseline target |
| CP4 | pending | Arc-length prototype on a small benchmark | Legacy load-control regression |
| CP5 | pending | TS-N65 arc-length application | Strict baseline and load-displacement curve |
| CP6 | pending | Model-level stabilization if needed | Sensitivity data, physics-change note |

## CP0: Experiment Harness

Add `scripts/summarize_oas_run.py` and validate it against the existing local
baseline replication:

```bash
scripts/summarize_oas_run.py \
  results/tsn65-full-baseline-replication-20260520-123252/runs/current-branch-full-baseline/solver.out \
  --baseline tsn65
```

Pass criteria:

- Parser reports step sequence `6,6,10,13,17,183,187,163`.
- Parser reports total nonlinear rows `585`.
- Parser reports `end_of_calculation=True`, no NaNs, no warnings, and no fallback
  acceptance.

## CP1: Existing Globalization Sweep On Legacy `tangent`

Run strict fixed-step variants with the full baseline deck. Do not introduce a
new algorithm in this checkpoint.

Variants:

- `G1-baseline`: controls off.
- `G1-backtracking-frozen-actual`: snapshot rollback, backtracking,
  `frozen_then_actual`, mixed merit, `min_alpha=0.03125`, `max_trials=6`,
  `cutback_on_fail=0`.
- `G1-backtracking-actual`: same as above, but actual trial evaluation.
- `G1-bisection-frozen-actual`: same as above, but bisection line search.
- `G1-fixed-damping-05`: fixed damping factor `0.5`.
- `G1-adaptive-damping`: adaptive damping, initial factor `0.5`.
- `G1-trust-stepnorm`: step-norm trust region, radius initial `1e-3`, min
  `1e-8`, max `1`, shrink `0.5`, expand `2`.
- `G1-line-search-adaptive-K`: backtracking frozen/actual plus adaptive matrix
  rebuild on small alpha and merit growth.

Promote only variants that complete `8/8` and reduce total rows or the step
6-8 tail compared with `585`.

## CP2: Adaptive Cutback/Stagnation

Use the best CP1 variant and enable adaptive stepping:

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

Pass criteria:

- Reaches `total_time = 1e-2`.
- No fallback accepted steps.
- No NaNs.
- Report minimum accepted `dt`, cutback count, accepted-step count, and final
  load/displacement gauges.

## CP3: Indirect Displacement Control Before Arc-Length

Run the existing `src/benchmark/Indirect_Control` benchmark first. Then test
TS-N65 with coordinate-based indirect control using the `v01 uz` gauge:

```text
indirect_control 2
ic_xcoords 2.1 2.1
ic_ycoords 0 0
ic_zcoords 0 -0.4
ic_directions 2 2
ic_displ_weights -1 1
ic_force_weights 0 0
ic_function 2
```

Also run the flipped sign. Select the sign whose first increment matches the
baseline load/displacement direction.

Pass criteria:

- IDC benchmark still passes.
- TS-N65 IDC reaches at least the same physical target as the strict 8-quarter
  baseline.
- Report `idc_time`, controlled displacement, load gauge, rows, cutbacks, and
  curve comparability.

## CP4: Arc-Length Prototype On A Small Benchmark

Add disabled-by-default controls:

```text
nonlinear_control load|indirect|arc_length
arc_length_radius_initial
arc_length_radius_min
arc_length_radius_max
arc_length_psi
arc_length_shrink 0.5
arc_length_expand 1.2
arc_length_target_iterations 8
arc_length_max_iterations 40
arc_length_constraint spherical
arc_length_sign_strategy previous_increment
```

V1 supports proportional external loading only:

```text
f_ref = f_ext(load=1) - f_ext(load=0)
```

Use a Schur-complement corrector:

- solve `K a = R`;
- solve `K b = f_ref`;
- compute `dlambda` from the linearized spherical arc-length constraint;
- trial-evaluate `(u + alpha du, lambda + alpha dlambda)` using the existing
  snapshot/rollback and line-search machinery.

Pass criteria:

- `nonlinear_control load` reproduces legacy behavior.
- Small benchmark runs with arc-length and writes `lambda`, radius, iterations,
  and controlled displacement.
- Rejected arc-length trials restore solver and material state.

## CP5: TS-N65 Arc-Length Application

Apply arc-length to TS-N65 only after CP4 passes.

Runs:

- strict baseline controls off;
- TS-N65 arc-length with legacy `tangent`;
- TS-N65 arc-length with best CP1 globalization;
- optional archived CSL tangent diagnostic variant, clearly marked as
  non-main-path.

Pass criteria:

- Reaches the same physical target as the strict baseline.
- No fallback acceptance.
- Report load factor, displacement gauges, radius changes, total rows, wall
  time, and comparison to `585`.

## CP6: Model-Level Stabilization

Enter only if CP5 still fails or needs tiny radii/alphas.

Order:

- material local substepping diagnostics;
- viscosity/rate-regularization sensitivity;
- nonlocal/crack-band regularization investigation.

These are modeling changes. Every report must state that physics changed and
include sensitivity data.

## Reporting Rules

- Store run folders under ignored `results/*/runs/`.
- Commit compact `report.md`, TSV summaries, scripts, and roadmap docs.
- Every report must include exact commit, executable path, solver deck, command,
  step table, total rows, wall time, warnings, NaNs, cutbacks, and verdict.
