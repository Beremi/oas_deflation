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

Checkpoint 2: targeted LDPM/CSL material tangent audit tests.

Tasks:

1. Identify the material/status classes used by the top `LDPMTetra` mismatch rows.
2. Add local material-level finite-difference tests for the stress update and tangent tensor.
3. Test loading, unloading, damage growth, frozen-state evaluation, and actual-state evaluation.
4. Compare tangent branches:
   - `tangent`,
   - `consistent`,
   - secant-like fallback if available.
5. Determine whether the mismatch is due to:
   - missing damage derivative,
   - wrong branch,
   - sign/scaling error,
   - state mutation during tangent/stress evaluation,
   - mismatch between residual update and stiffness update.
6. Only after local material tests pass, rerun strict 8 quarter steps and compare to the saved baseline.

## Practical warnings for next agent

- Do not commit copied TS-N65 result decks unless explicitly requested; they can be large.
- Do force-add compact Markdown reports if they live under ignored `results/`.
- Always separate diagnostic runs from phase-closure runs.
- If a run reaches huge DFGMRES norms or stalls with no nonlinear output for more than about 10 minutes in a clearly blown-up state, kill it and mark the phase failed/stalled.
- Do not interpret `merit < 1` as convergence unless residual, displacement, and energy convergence criteria are all satisfied.
- Do not use line-search acceptance as a replacement for normal convergence.
- The next likely high-value fix is not another alpha policy; it is local constitutive tangent correctness for LDPM/CSL.
