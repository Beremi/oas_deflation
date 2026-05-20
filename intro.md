# Nonlinear solver handoff intro

This repository is currently being worked on in branch `nonlinear_solver_testing`.

Read this file first when picking up the nonlinear solver work. It summarizes what has been implemented, what has been tested, what did and did not help, and the next useful checkpoint.

## Read order

1. `intro.md`
2. `docs/nonlinear_solver.md`
3. `docs/nonlinear-solver-stabilization-handoff.md`
4. `docs/nonlinear_solver_next_steps_workplan.md`
5. `results/nonlinear-solver-stabilization-20260518-093732/report.md`
6. `results/nonlinear-solver-next-steps-20260518-180841/report.md`
7. `results/nonlinear-solver-tangent-attribution-20260519-113815/report.md`
8. `results/nonlinear-solver-material-tangent-audit-20260519-143422/report.md`
9. `results/tsn65-zip-runs-20260519-190000/report.md`
10. `results/tsn65-csl-analytic-tangent-20260519-230004/report.md`

The `results/` directory is mostly ignored because it contains copied TS-N65 decks and large runtime files. The important Markdown reports should be kept or explicitly copied/force-added when making checkpoint commits.

## Main problem being solved

The target problem is the TS-N65 nonlinear mechanics benchmark. The problematic behavior is nonlinear Newton-like solve instability or very slow convergence in fixed load steps. Some full steps need hundreds of nonlinear iterations, and some hit the maximum iteration count.

The nonlinear solver solves a residual equilibrium problem of the form:

```text
R(u, history, load) = f_ext(load) - f_int(u, material_history) = 0
```

The local constitutive model is path/state dependent. This matters because trial evaluations can mutate material status unless state is restored correctly.

The key practical question is whether poor convergence is caused by:

- poor Newton step globalization,
- stale or wrong stiffness/tangent matrix,
- bad/inconsistent material tangent,
- load-control/path-following instability,
- local constitutive snap/localization that needs material or model-level changes.

The current evidence points most strongly at a material/tangent mismatch in LDPM/CSL-related `LDPMTetra` elements.

## Benchmark settings to preserve

Strict TS-N65 quarter-step phase closure:

```text
time_step = 1.25e-3
min_time_step = 1.25e-3
max_time_step = 1.25e-3
total_time = 1e-2
limit_tolerance = 0
max_iterations = 300
j = 16
solver_type = DeflatedFGMRES
dfgmres_tolerance = 1e-1
dfgmres_true_tolerance = 1e-1
dfgmres_deflation_vectors = 20
stiff_matrix_type tangent
stiffness_matrix_iter_update = 10
```

Run shape:

```bash
OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores OAS -j 16 master.inp
```

Saved baseline target from earlier notes:

```text
strict quarter steps: 8/8 converged
total nonlinear iterations: 585
total duration: 26:53
```

Every future checkpoint must finish a strict 8-quarter-step run with diagnostics disabled and compare to this target. If this cannot be reproduced, report that explicitly before claiming an improvement.

## Build used during latest checkpoint

The in-tree `build/` directory was not usable. An external build was used:

```bash
cmake -S . -B /tmp/oas_deflation_build_hypre_novtk -DCMAKE_BUILD_TYPE=Release -DUSE_VTK=OFF -DUSE_HYPRE=ON
cmake --build /tmp/oas_deflation_build_hypre_novtk -j 16
```

Executable:

```text
/tmp/oas_deflation_build_hypre_novtk/bin/OAS
```

VTK was disabled because system VTK required a missing `CLI11` package.

## Current implemented solver controls

Several local nonlinear stabilization tools have already been added behind disabled-by-default input controls:

- fixed/adaptive damping,
- rollback-adaptive damping,
- backtracking line search,
- bisection-style line search,
- line-search trial evaluation modes,
- stagnation/growth cutback,
- adaptive matrix rebuild triggers,
- step-norm trust region,
- material snapshot/rollback API for safe trial restore,
- global tangent finite-difference diagnostic,
- element/material tangent-attribution diagnostic.

The intended default remains legacy behavior when new nonlinear globalization or diagnostic controls are omitted or disabled.

## What did not materially improve TS-N65

The following were tested and did not give a convincing first-phase improvement:

- fixed damping,
- adaptive damping,
- rollback-adaptive damping,
- backtracking line search,
- bisection line search,
- step-norm trust region,
- alpha caps or alpha starts,
- damping plus adaptive rebuild in tested settings.

Main reason: changing the Newton increment length often either did almost the same as baseline Newton or failed/stalled. The hard steps appear more consistent with a tangent/model issue than a simple globalization issue.

## Best diagnostic evidence so far

Global tangent check at TS-N65 hard state:

```text
case: results/nonlinear-solver-tangent-attribution-20260519-113815/diag-step6-it10-localrestore
step: 6
iteration: 10
eps: 1e-6
matrix type: tangent
```

Result:

```text
direction   relative error   cosine
newton      5.320964         0.707442
random_1    0.032998         0.999465
```

Interpretation:

- The assembled tangent is very poor along the actual Newton direction.
- The same tangent is much better for a deterministic random direction.
- This suggests a path/state/localization direction issue rather than a simple global sign error.

Element attribution result:

- top mismatch elements are `LDPMTetra`,
- not beam elements,
- not obvious boundary/control elements.

The next useful work is therefore a targeted LDPM/CSL tangent audit.

## Latest checkpoint: LDPM/CSL material tangent audit

Checkpoint 2 is partially implemented:

- `OAS_material_tangent_audit` runs local one-point CSL and LDPM stress/tangent finite differences.
- `MaterialTangentAudit` is registered with CTest and passes expected linear/frozen checks.
- LDPM-family and coupled CSL material statuses now have typed snapshot/restore/hash support.
- Element attribution output now includes `material_id`, `material_name`, and `material_status_names`.
- CSL and LDPM now support explicit `consistent` branches backed by state-restored numerical material derivatives.

Key local audit result:

- CSL elastic, frozen, and damaged-unloading tangent paths close to finite differences.
- CSL active `damage_growth` with actual state update does not close: secant/tangent relative error is about `3.4089` with cosine about `-0.9963`.
- CSL active `damage_growth` with `consistent` closes: relative error is about `3.01e-5` with cosine near `1`.
- LDPM `damage_growth` and unloading paths close in the local one-point path with `secant`.
- The selected LDPM `elastic_loading` probe shows a large `secant` mismatch, while LDPM `consistent` closes it with relative error about `2.32e-6`.

Strict TS-N65 phase closure was not rerun during this material-audit checkpoint because the existing CSL `tangent`/`secant` and LDPM `secant` branches still had selected active-path mismatches, and the TS-N65 hard-state deck/state was not yet present in this checkout. The later zip benchmark checkpoint below completed the hard-state `element_top` rerun with material/status columns.

## Latest checkpoint: TS-N65 benchmark zip run

The newly added `TS-N_65.zip` and `Dogbone.zip` benchmark archives were inspected and extracted under:

```text
results/tsn65-zip-runs-20260519-190000/input/
```

`Dogbone` is a small deck with previous archived outputs. `TS-N_65` is the current large nonlinear target deck and was run for this checkpoint.

The original TS-N65 deck cannot be run as-is in this build because its `solver.inp` requests:

```text
solver_type PardisoLDLT
```

and the local executable reports:

```text
Solver type PardisoLDLT is not implemented
```

The strict TS-N65 hard-state rerun with `DeflatedFGMRES`, `stiff_matrix_type tangent`, and `stiffness_matrix_iter_update 10` reached step 6, iteration 10 and reproduced the global Newton-direction tangent mismatch:

```text
relative_error = 5.320964
cosine = 0.707442
```

The new material/status columns identify the top element-local mismatch rows as:

```text
element_name = LDPMTetra
material_id = 0
material_name = CSL material
material_status_names = CSL mat. statusx12
```

Top-50 element-local comparison at the same hard state:

| matrix type | mean relative error | max relative error | mean cosine | mean mismatch norm | max mismatch norm |
|---|---:|---:|---:|---:|---:|
| `tangent` | 3.75268 | 7.00578 | -0.308147 | 213602 | 358831 |
| diagnostic `consistent` | 0.133558 | 0.466171 | 0.986366 | 939.346 | 10748.7 |

This ties the benchmark hard-state mismatch to CSL active damage tangent behavior. However, a production smoke run with global `stiff_matrix_type consistent` did not complete the first quarter step: the first DFGMRES solve hit 500 iterations with true relative residual `0.843808` against the requested `0.1`, and the process was stopped after reaching about 27 GB RSS.

Current implication: use the numerical `consistent` branch as a reference for deriving a cheaper analytical CSL active-damage tangent, but do not switch TS-N65 production runs globally to `consistent` yet.

## Latest checkpoint: archived analytical CSL active-damage tangent

An analytical CSL active-damage tangent was tested by temporarily wiring it into the CSL `tangent` branch:

```text
D_tangent = (1 - damage) * D_elastic - (D_elastic * strain) outer d_damage/d_strain
```

It improved the local diagnostics, but it did not close the strict TS-N65 phase benchmark. The code is now archived behind the explicit branch name `archived_csl_damage_tangent`; the default CSL `tangent` path is back to the legacy degraded elastic stiffness.

Local audit result:

| CSL damage-growth branch | relative error | cosine |
|---|---:|---:|
| `secant` | 3.4089098751 | -0.9963323184 |
| `tangent` | 3.4089098751 | -0.9963323184 |
| `archived_csl_damage_tangent` | 2.3769400931e-5 | 0.99999999997 |
| `consistent` | 3.0108885587e-5 | 0.99999999955 |

TS-N65 hard-state element-top result with step-start quasi-Newton matrix updates:

| matrix type | mean relative error | max relative error | mean cosine | mean mismatch norm | max mismatch norm |
|---|---:|---:|---:|---:|---:|
| legacy `tangent` | 3.75268 | 7.00578 | -0.308147 | 213602 | 358831 |
| archived analytical branch | 0.0182187 | 0.638459 | 0.995489 | 563.047 | 22422.7 |

The archived branch fixes the local and element-local CSL active-damage tangent mismatch. The strict phase closure is not fixed:

- `stiffness_matrix_iter_update 10` with the analytical tangent destabilized step 3 after the iteration-10 matrix rebuild.
- Step-start analytical tangent with diagnostics disabled converged 5/8 steps and failed at step 6 after 300 nonlinear iterations.
- Step-start analytical tangent plus actual-state backtracking converged 6/8 steps, including step 6, but stalled in step 7 and was stopped.

Do not claim solver improvement from the analytical tangent alone. It is now an archived reproduction branch, not mainline behavior. The current evidence says the defect moved from a local CSL tangent mismatch to a global softening/load-control/path-following problem.

## Latest checkpoint: element/material tangent attribution

New input keywords:

```text
nonlinear_tangent_check_scope global|element_top
nonlinear_tangent_check_top_elements 20
nonlinear_tangent_check_element_output tangent_check_elements.tsv
nonlinear_tangent_check_matrix_type current
```

Behavior:

- `scope=global` keeps the existing global tangent check behavior.
- `scope=element_top` writes top local element mismatches ranked by `||Ke p_e - fd_e||`.
- Local finite difference compares local element internal force changes against local analytical stiffness contribution.
- Element-local material snapshots are used inside the attribution loop to avoid full-model restore per element.

Important caveat:

- Use `nonlinear_tangent_check_stop_after 1` for hard-state `element_top` runs unless specifically auditing non-intrusiveness.
- Pair every diagnostic-enabled run with a separate diagnostic-disabled phase-closure run.

## Latest phase closure result

A diagnostic-disabled strict 8-quarter-step run was launched after the packet 1 implementation:

```text
case: results/nonlinear-solver-tangent-attribution-20260519-113815/phase-close-disabled
nonlinear_tangent_check 0
```

It did not reproduce the saved baseline target. It stalled in step 6 after iteration 132 with residual near 1, repeated DFGMRES warnings, and no new nonlinear output for more than 10 minutes. The run was killed and reported as non-converged/stalled.

This means:

- do not claim the latest checkpoint improved solver convergence;
- the diagnostic code built and produced useful data;
- the current local executable/settings need baseline reproduction checked before performance claims.

## Files changed in latest checkpoint

Code:

```text
src/solver/src/solver_implicit.h
src/solver/src/solver_implicit.cpp
```

Docs/reports:

```text
docs/nonlinear_solver.md
docs/nonlinear_solver_next_steps_workplan.md
results/nonlinear-solver-tangent-attribution-20260519-113815/report.md
intro.md
```

## Blind or diagnostic-only code paths

Diagnostic-only:

- `nonlinear_tangent_check`
- `nonlinear_tangent_check_scope element_top`
- `nonlinear_tangent_check_element_output`
- `nonlinear_tangent_check_matrix_type`

These should not affect production behavior when `nonlinear_tangent_check 0`.

Blind/experimental solver-globalization paths already in the branch:

- damping variants,
- line search variants,
- trust-region variants,
- adaptive rebuild variants.

These are disabled by default and should not be treated as proven TS-N65 improvements.

## Next checkpoint

The local CSL tangent defect is understood and the analytical fix is archived. The next useful work should test solver/model strategies that can handle the softening path rather than simply replacing the default tangent:

1. Try path-following or indirect displacement control on the strict 8-quarter TS-N65 phase run.
2. Try a controlled positive-definite tangent blend/limiter using the archived branch as a diagnostic reference.
3. Consider material-level softening regularization or substepping if solver controls alone do not close the phase.
4. Keep every experiment gated by the diagnostic-disabled strict `8/8`, `585 iteration`, `26:53` baseline comparison.

## Practical warnings for next agent

- Do not commit copied TS-N65 result decks unless explicitly requested; they can be large.
- Do force-add compact Markdown reports if they live under ignored `results/`.
- Always separate diagnostic runs from phase-closure runs.
- If a run reaches huge DFGMRES norms or stalls with no nonlinear output for more than about 10 minutes in a clearly blown-up state, kill it and mark the phase failed/stalled.
- Do not interpret `merit < 1` as convergence unless residual, displacement, and energy convergence criteria are all satisfied.
- Do not use line-search acceptance as a replacement for normal convergence.
- The next likely high-value fix is not another alpha policy; it is local constitutive tangent correctness for LDPM/CSL.
