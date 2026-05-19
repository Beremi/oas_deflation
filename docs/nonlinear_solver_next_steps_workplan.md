# Next Nonlinear Solver Workplan: Matrix/Tangent Diagnosis and Safe Trial State

## Scope

This workplan covers the next local nonlinear-solver diagnostics and infrastructure after damping, backtracking, and trust-region step scaling showed limited improvement on TS-N65. The focus is to decide whether the main problem is the Newton matrix/tangent model, matrix rebuild cadence, or unsafe trial material-state mutation.

Large path-following or constitutive-model changes are explicitly out of scope for this pass: arc-length, PETSc SNES, material-law substepping, viscosity/regularization, and nonlocal damage regularization.

## Benchmark

Use TS-N65 strict 8 quarter steps:

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
```

All Newton globalization controls remain disabled for points 1 and 2.

## Result location

Primary report:

```text
results/nonlinear-solver-next-steps-<timestamp>/report.md
```

Machine-readable files:

```text
results/nonlinear-solver-next-steps-<timestamp>/stiffness_sweep.tsv
results/nonlinear-solver-next-steps-<timestamp>/matrix_type_sweep.tsv
results/nonlinear-solver-next-steps-<timestamp>/tangent_check.tsv
```

## Point 1: stiffness rebuild cadence sweep

Run strict 8 quarter steps with:

```text
stiffness_matrix_iter_update 1
stiffness_matrix_iter_update 3
stiffness_matrix_iter_update 10
```

Fixed settings:

```text
stiff_matrix_type secant
first_iteration_stiff_matrix_type void
```

Selection rule:

1. Prefer strict convergence of all 8 steps.
2. Then prefer lowest summed loading-step wall time.
3. If summed times are within 10%, prefer fewer matrix rebuilds.

Report per variant:

| variant | converged steps | iterations by step | summed step time | total OAS time | final errors | notes |
|---|---:|---|---:|---:|---|---|

## Point 2: matrix/tangent type sweep

Use the best cadence from point 1. Run:

```text
stiff_matrix_type secant
stiff_matrix_type elastic
stiff_matrix_type tangent
stiff_matrix_type consistent
first_iteration_stiff_matrix_type elastic + stiff_matrix_type secant
first_iteration_stiff_matrix_type elastic + stiff_matrix_type tangent
```

Unsupported or erroring variants are recorded explicitly and not silently modified.

Report per variant:

| variant | supported | converged steps | iterations by step | summed step time | total OAS time | final errors | notes |
|---|---|---:|---|---:|---:|---|---|

## Point 3: finite-difference tangent diagnostic

Add disabled-by-default input controls:

```text
nonlinear_tangent_check 0|1
nonlinear_tangent_check_step <int>
nonlinear_tangent_check_iteration <int>
nonlinear_tangent_check_eps <double>
nonlinear_tangent_check_random_vectors <int>
nonlinear_tangent_check_include_newton 0|1
nonlinear_tangent_check_stop_after 0|1
nonlinear_tangent_check_output tangent_check.tsv
```

At the selected nonlinear state, compare:

```text
Keff * p
```

against:

```text
[R(u) - R(u + eps*p)] / eps
```

where `R = f_ext - f_int`. The load is held fixed and solver/material state is restored after each perturbation.

Diagnostic states:

```text
step 6, iteration 10
step 6, iteration 80
step 7, iteration 10
```

Epsilon values by separate runs:

```text
1e-4
1e-5
1e-6
```

Directions:

```text
current Newton increment
deterministic random reduced vector(s)
```

Report:

| step | iteration | eps | direction | relative error | cosine | ||Kp|| | ||fd|| |
|---:|---:|---:|---|---:|---:|---:|---:|

## Point 4: material-status snapshot and rollback

Add snapshot/restore infrastructure:

```cpp
MaterialStatus::cloneState()
MaterialStatus::restoreStateFrom(...)
Element::createMaterialStatusSnapshot(...)
Element::restoreMaterialStatusSnapshot(...)
ElementContainer::createMaterialStatusSnapshot(...)
ElementContainer::restoreMaterialStatusSnapshot(...)
```

Implement first for TS-N65-relevant statuses:

```text
CSLMaterialStatus
NormalPlasticBeamMaterialStatus
required mechanical base status classes
```

Add disabled-by-default safe rollback control:

```text
nonlinear_material_snapshot_rollback 0|1
```

Acceptance for point 4:

1. Snapshot API builds cleanly.
2. Roundtrip restore works in a small CSL/beam-containing run or diagnostic trial.
3. Future line-search/trust-region rejected trials can use full material-state restore rather than only `resetMaterialStatuses()`.

## Final conclusions to report

The final report must state:

- Whether frequent stiffness rebuilds improve TS-N65 strict quarter-step convergence enough to justify cost.
- Which matrix/tangent type is best among supported local options.
- Whether current `Keff` is a good derivative of the residual at hard states.
- Whether safe material snapshot/rollback is available for future globalization work.
- Which remaining blockers are local solver issues versus constitutive/path-following issues.
