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
| CP1 | partial | Existing globalization sweep on legacy `tangent` | Strict fixed-step baseline |
| CP2 | complete-negative | Adaptive cutback/stagnation run | Strict baseline plus adaptive-step metrics |
| CP3 | complete-negative | Existing indirect displacement control before arc-length | IDC benchmark and TS-N65 baseline target |
| CP4 | complete | Arc-length prototype on a small benchmark | Legacy load-control regression and TS-N65 strict baseline gate |
| CP5 | complete-negative | TS-N65 arc-length application | Strict baseline and load-displacement curve |
| CP5b | complete-negative | Gauge-constrained arc-length continuation | Strict baseline gate and TS-N65 gauge/row comparison |
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

Current checkpoint result:

- Report: `results/tsn65-adaptive-cutback-20260520-cp2/report.md`.
- Executed variants: `G2-stagnation-only`, `G2-backtracking-actual`.
- Both variants failed before `total_time = 1e-2` after reaching
  `min_time_step = 1.25e-5`.
- No NaNs and no fallback acceptance were observed.
- `G2-stagnation-only` reached furthest accepted time `7.421875e-3`.
- `G2-backtracking-actual` reached furthest accepted time `7.8125e-3`, with
  minimum line-search alpha `0.125`.
- Verdict: adaptive cutback and line search reduce the local row burn, but are
  not sufficient to pass TS-N65 load control. Proceed to CP3 indirect
  displacement control before arc-length.

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
ic_function 4
```

Also run the flipped sign. Select the sign whose first increment matches the
baseline load/displacement direction.

Note: OAS function indices are zero-based. In the TS-N65 deck, index `2` is the
steel material envelope, so the CP3 runner appends a local displacement-target
function and uses `ic_function 4`.

Pass criteria:

- IDC benchmark still passes.
- TS-N65 IDC reaches at least the same physical target as the strict 8-quarter
  baseline.
- Report `idc_time`, controlled displacement, load gauge, rows, cutbacks, and
  curve comparability.

Current checkpoint result:

- Report: `results/tsn65-indirect-control-20260520-cp3/report.md`.
- Script: `scripts/run_tsn65_indirect_control.py`.
- The historical `src/benchmark/Indirect_Control` benchmark did not pass on the
  current checkpoint tree: it exited nonzero after `146` steps, `21910`
  nonlinear rows, and `138` fallback acceptances under its legacy benchmark
  tolerance setup.
- TS-N65 sign calibration selected `ic_displ_weights -1 1`; the flipped sign
  failed immediately by driving the function query below the time range.
- Raw strict TS-N65 IDC completed `8/8` fixed quarter steps with no fallback,
  no NaNs, and no cutbacks, but undershot the physical target: final `v01 =
  1.233328e-5`, `18.84%` below the strict baseline `1.519697e-5`.
- A final-v01 calibrated IDC target, scaled by `1.23219208515`, was tried to
  reach the same physical `v01` target. It stalled in step `7`; it was stopped
  at iteration `69` with residual `2.352053e-3`, displacement `4.419994e-4`,
  and energy `4.186416e-3`.
- Verdict: existing IDC is not a baseline replacement for TS-N65. It can follow
  a lower-amplitude controlled path, but reaching the strict baseline gauge
  target still exposes late-step nonlinear instability. Proceed to CP4
  arc-length on a small benchmark before applying continuation to TS-N65.

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

Current checkpoint result:

- Report: `results/arc-length-cp4-20260521-cp4/report.md`.
- Summary TSV: `results/arc-length-cp4-20260521-cp4/cp4_arc_length.tsv`.
- Script: `scripts/run_arc_length_checkpoint.py`.
- Implemented disabled-by-default controls:
  - `nonlinear_control load|indirect|arc_length`;
  - `arc_length_radius_initial`, `arc_length_radius_min`,
    `arc_length_radius_max`, `arc_length_psi`, `arc_length_shrink`,
    `arc_length_expand`, `arc_length_target_iterations`,
    `arc_length_max_iterations`, `arc_length_constraint`,
    `arc_length_sign_strategy`.
- V1 arc-length is limited to proportional nodal loading with reference load
  `f_ext(load=1)-f_ext(load=0)`. Because OAS boundary conditions are active on
  `[begin, end)`, the implementation samples the left limit at `1.0`.
- Small proportional-load Timoshenko benchmark results:
  - no-keyword legacy load-control: passed, `3` steps, `6` nonlinear rows;
  - explicit `nonlinear_control load`: passed, `3` steps, `6` rows;
  - arc-length prototype: passed, `3` steps, `6` rows, final lambda `0.3`;
  - arc-length with line-search/snapshot rollback: passed, `4` steps, `18`
    rows, final lambda `0.3`, minimum alpha `0.5`, rejected alpha trials logged
    and restored.
- Strict TS-N65 baseline gate after the CP4 code change:
  - Report: `results/tsn65-cp4-baseline-gate-20260521-0058/report.md`.
  - Exact sequence matched: `6, 6, 10, 13, 17, 183, 187, 163`.
  - Total nonlinear rows: `585`.
  - Verdict: neutral/no regression, no fallback acceptance, no NaNs, no cutbacks.

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

Current checkpoint result:

- Reports:
  - `results/tsn65-arc-length-cp5-20260521/report.md`;
  - `results/tsn65-arc-length-cp5-radcal-20260521/report.md`;
  - `results/tsn65-cp5-baseline-gate-20260521/report.md`.
- Script: `scripts/run_tsn65_arc_length_cp5.py`.
- Added optional arc-length reference controls:
  - `arc_length_reference proportional_load|finite_difference`;
  - `arc_length_reference_delta`.
- Defaults preserve CP4 behavior: `arc_length_reference proportional_load`,
  `arc_length_reference_delta 1`.
- Direct proportional-load arc-length on TS-N65 is unsupported because the
  TS-N65 deck is prescribed-displacement driven and has no nonzero reduced
  nodal reference load `f_ext(load=1)-f_ext(load=0)`.
- The finite-difference reference mode can derive a direction from the residual
  change caused by prescribed-BC/load parameter perturbation, but the tested
  radii were not competitive:
  - radius `1.25e-3`: failed in the first arc step after `3` nonlinear rows;
    final lambda drifted back to `6.433438e-6`;
  - radius `1.35e-2`: calibrated the first predictor to `lambda=0.00125`, but
    the corrector oscillated and needed `39` rows before the first strict
    baseline quarter step was stably reached; strict baseline step 1 takes `6`
    rows.
- Verdict: CP5 arc-length is not a speed improvement for the strict TS-N65
  eight-quarter baseline. It should remain experimental until the continuation
  formulation is improved, likely by adding a monotone displacement-control
  constraint or a better scaled corrector instead of the current spherical norm
  on the full displacement field.

## CP5b: Gauge-Constrained Arc-Length Continuation

CP5b adds a TS-N65-specific arc-length continuation mode that keeps the
arc-length load factor as the equilibrium load parameter but replaces the
spherical corrector constraint with an indirect-control displacement gauge:

```text
nonlinear_control arc_length
arc_length_constraint gauge
arc_length_reference finite_difference
arc_length_sign_strategy monotone_lambda
arc_length_auto_radius 1
arc_length_gauge_tolerance 1e-3
indirect_control 2
ic_xcoords 2.1 2.1
ic_ycoords 0 0
ic_zcoords 0 -0.4
ic_directions 2 2
ic_displ_weights -1 1
ic_force_weights 0 0
```

New disabled-by-default controls added in this checkpoint:

- `arc_length_constraint gauge`;
- `arc_length_sign_strategy monotone_lambda`;
- `arc_length_auto_radius`;
- `arc_length_gauge_tolerance`.

The gauge constraint currently supports displacement-only indirect-control
blocks. It rejects force-weighted IDC blocks because the arc corrector does not
yet consistently linearize reaction-force gauges.

Current checkpoint result:

- Report: `results/tsn65-gauge-arc-cp5b-20260521-094117/report.md`.
- Strict baseline gate after CP5b code:
  `results/tsn65-cp5b-baseline-gate-20260521/report.md`.
- Summary TSV:
  `results/tsn65-gauge-arc-cp5b-20260521-094117/cp5b_gauge_arc_length.tsv`.
- Script: `scripts/run_tsn65_gauge_arc_cp5b.py`.
- First-quarter uncalibrated gauge arc-length completed, but the internal IDC
  gauge and exported `LD.out` interpolation gauge differ. The run hit the
  internal target exactly, while exported `v01` was `17.07%` low.
- First-quarter calibrated gauge arc-length matched the exported strict
  baseline `v01 = 2.613423e-6`, with no fallback acceptance or NaNs. It needed
  `25` nonlinear rows and about `4:08`, compared with strict baseline step 1:
  `6` rows and about `0:35`. The existing line-search option was neutral for
  rows and slightly slower in wall time.
- A full calibrated gauge-arc run was attempted. It reached reported step `7`
  and `235` nonlinear rows before being manually stopped because step `7`
  drifted away from convergence. The last parsed row was step `7`, iteration
  `89`, with residual `3.963647e-3`, displacement `2.671126e-3`, and energy
  `4.556947e-3`.
- The strict baseline gate still matches exactly after the CP5b code changes:
  `8/8` steps, sequence `6,6,10,13,17,183,187,163`, total `585` rows, no
  warnings, no NaNs, no cutbacks, and no fallback acceptance.
- The full calibrated run followed a softer linear gauge path than the strict
  baseline until the final target. It is therefore a continuation experiment,
  not an exact loading-history replacement.
- Verdict: gauge-constrained arc-length is now mechanically usable enough to
  hit prescribed displacement gauges, but it is not a strict-baseline speed
  replacement. The next usable direction is either a baseline-gauge PWL control
  path with stagnation cutback/adaptive matrix rebuild, or a deeper
  path-following formulation that constrains the same interpolated exporter
  gauge rather than nearest-node IDC coordinates.

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
