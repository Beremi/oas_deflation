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

Mandatory phase closure benchmark:

- Every diagnostic or implementation phase must finish with the full strict 8 quarter-step TS-N65 run.
- Compare one-to-one against the current best local baseline: `stiff_matrix_type tangent`, `stiffness_matrix_iter_update 10`, strict convergence `8/8`, total nonlinear iterations `585`, total OAS duration `26:53`.
- If a phase only adds diagnostics and should not alter solver behavior, the strict run must either match this baseline within normal runtime noise or clearly explain any overhead from enabled diagnostics.
- A phase is not considered finished until the report includes the 8-step table, total iterations, total wall time, convergence count, final errors, and a short baseline delta.

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
results/nonlinear-solver-next-steps-<timestamp>/tangent_check_elements.tsv
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
nonlinear_tangent_check_scope global|element_top
nonlinear_tangent_check_top_elements 20
nonlinear_tangent_check_element_output tangent_check_elements.tsv
nonlinear_tangent_check_matrix_type current
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

For `nonlinear_tangent_check_scope element_top`, also report the top local contributors:

| step | iteration | eps | direction | element id | element name | material id | material name | material status names | matrix type | relative error | cosine | ||Kp|| | ||fd|| | mismatch |
|---:|---:|---:|---|---:|---|---:|---|---|---|---:|---:|---:|---:|---:|

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

## Packet 1 execution note

- Every future checkpoint remains incomplete until it has a strict 8 quarter-step TS-N65 phase-closure run with diagnostics disabled and a table against the saved baseline target (`8/8`, `585` nonlinear iterations, `26:53` total duration).
- Element-level tangent attribution is diagnostic-only. If it is run at a hard state, prefer `nonlinear_tangent_check_stop_after 1` unless a separate non-intrusiveness audit proves that continuing after the local element FD sweep leaves all material and element-side state unchanged.
- If a diagnostic-enabled run is needed at a hard state, pair it with a separate diagnostic-disabled phase-closure run. The diagnostic run identifies the tangent mismatch source; the disabled run verifies production behavior.

## Packet 2 execution note

- `OAS_material_tangent_audit` now provides local one-point CSL/LDPM stress-update finite differences and is registered as `MaterialTangentAudit`.
- The local audit currently passes expected elastic/frozen checks and shows that numerical `consistent` branches close the selected active CSL and LDPM probes while the existing CSL `tangent`/`secant` and LDPM `secant` branches do not in those paths. Do not claim TS-N65 improvement until a hard-state `element_top` rerun with material/status columns and a diagnostic-disabled phase closure are complete.

## Packet 3 execution note

- The newly added `TS-N_65.zip` and `Dogbone.zip` archives were inspected. `TS-N_65` is the large target deck and `Dogbone` is a smaller deck with archived result files.
- The original TS-N65 deck requests `solver_type PardisoLDLT`, which is not implemented in the current build, so strict `DeflatedFGMRES` settings were used for the diagnostic run.
- The TS-N65 step-6/iteration-10 `element_top` rerun now ties the previous top `LDPMTetra` rows to `CSL material` with `CSL mat. statusx12`.
- At that hard state, the top-50 element-local `tangent` rows have mean relative error `3.75268`, mean cosine `-0.308147`, and max mismatch norm `358831`.
- Running the same element-top diagnostic with `nonlinear_tangent_check_matrix_type consistent` drops the top-50 mean relative error to `0.133558`, raises mean cosine to `0.986366`, and drops max mismatch norm to `10748.7`.
- A production smoke run with global `stiff_matrix_type consistent` assembled but the first DFGMRES solve hit 500 iterations with true relative residual `0.843808`; it was stopped after reaching about 27 GB RSS. Treat numerical `consistent` as the reference for deriving a cheaper CSL active-damage tangent, not as a production global matrix option yet.

## Packet 4 execution note

- An analytical CSL active-damage derivative was tested by temporarily wiring it into the `tangent` branch. It adds the rank-one damage term `-(D_elastic * strain) outer d_damage/d_strain` when CSL damage grows.
- That derivative is now archived as `archived_csl_damage_tangent`; the default CSL `tangent` branch is back to the legacy degraded elastic stiffness.
- Local `OAS_material_tangent_audit` shows CSL `damage_growth` / `archived_csl_damage_tangent` closes with relative error `2.37694e-5` and cosine approximately `1`, matching the numerical `consistent` reference. The default `tangent` and `secant` rows remain bad with relative error `3.4089`.
- TS-N65 hard-state `element_top` with the archived analytical branch reduces the top-50 mean relative error from `3.75268` to `0.0182187` and raises mean cosine from `-0.308147` to `0.995489`.
- Strict phase closure still fails. With `stiffness_matrix_iter_update -1`, the diagnostic-disabled 8-quarter run converged steps 1-5 and failed at step 6 after 300 nonlinear iterations. With actual-state backtracking line search it converged through step 6 but stalled in step 7 and was stopped.
- Therefore the analytical tangent is diagnostically valuable but not a standalone solver improvement. Future work should target path-following/globalization or tangent limiting/blending, and every candidate still must pass the strict `8/8`, `585 iteration`, `26:53` baseline comparison gate.
