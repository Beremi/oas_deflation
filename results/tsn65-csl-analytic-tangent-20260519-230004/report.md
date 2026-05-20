# TS-N65 CSL Analytical Damage-Tangent Experiment

Date: 2026-05-19
Branch: `nonlinear_solver_testing`
Build: `/tmp/oas_deflation_build_material_audit`

## Scope

This checkpoint tested an analytical CSL active-damage tangent by temporarily wiring it into the production `tangent` branch:

```text
sigma = (1 - damage(eps)) * D_elastic * eps
D_tangent = (1 - damage) * D_elastic - (D_elastic * eps) outer d_damage/d_eps
```

The derivative is only active when CSL damage grows. Unloading and frozen-damage paths keep the previous degraded elastic stiffness.

Checkpoint update on 2026-05-20: because the strict TS-N65 phase benchmark did not close, this code has been moved out of the normal path and archived as `archived_csl_damage_tangent`. The default CSL `tangent` branch is back to the legacy degraded elastic stiffness. Historical run directories below still say `stiff_matrix_type tangent` because those runs were made before the branch was archived.

## Local Material Audit

Commands:

```bash
cmake --build /tmp/oas_deflation_build_material_audit --target OAS_material_tangent_audit OAS -j 16
/tmp/oas_deflation_build_material_audit/bin/OAS_material_tangent_audit --output results/nonlinear-solver-material-tangent-audit-20260519-143422/material_tangent_audit.tsv
ctest --test-dir /tmp/oas_deflation_build_material_audit -R MaterialTangentAudit --output-on-failure
```

Result:

```text
MaterialTangentAudit: Passed
```

Key CSL active damage-growth rows:

| branch | relative error | cosine | note |
|---|---:|---:|---|
| `secant` | 3.4089098751 | -0.9963323184 | old mismatch remains |
| `tangent` | 3.4089098751 | -0.9963323184 | legacy default mismatch restored |
| `archived_csl_damage_tangent` | 2.3769400931e-5 | 0.99999999997 | archived analytical tangent closes |
| `consistent` | 3.0108885587e-5 | 0.99999999955 | numerical reference |

## TS-N65 Hard-State Diagnostics

### Mid-Step Rebuild

Run directory:

```text
results/tsn65-csl-analytic-tangent-20260519-230004/runs/strict_tangent_element_top_step6_it10/
```

Settings:

```text
stiff_matrix_type tangent
stiffness_matrix_iter_update 10
nonlinear_tangent_check step 6 iteration 10
```

This was the temporary analytical `tangent` branch. To reproduce from the archived checkpoint, use `stiff_matrix_type archived_csl_damage_tangent`.

Result:

- The run did not reach step 6.
- It became unstable in step 3 immediately after the iteration-10 matrix rebuild.
- Last useful nonlinear row before stopping:

```text
step 3 iteration 15 residual 9.897274e-01 displacement 9.449099e-01 energy 1.552062e+00
```

Interpretation: the analytically correct local tangent is too aggressive/indefinite for the current load-control Newton path when rebuilt mid-step.

### Step-Start Quasi-Newton Diagnostic

Run directory:

```text
results/tsn65-csl-analytic-tangent-20260519-230004/runs/analytic_tangent_stepstart_element_top_step6_it10/
```

Settings:

```text
stiff_matrix_type tangent
stiffness_matrix_iter_update -1
nonlinear_tangent_check step 6 iteration 10
```

This was the temporary analytical `tangent` branch. To reproduce from the archived checkpoint, use `stiff_matrix_type archived_csl_damage_tangent`.

This variant reached the hard state in `00:04:35.559`.

Global Newton-direction check:

| step | iteration | eps | relative error | cosine | `||Kp||` | `||fd||` |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 10 | 1e-6 | 6.392516 | 0.7353447 | 2.036287e7 | 2.871311e6 |

Element-local top-50 summary:

| matrix type | rows | mean relative error | max relative error | mean cosine | mean mismatch norm | max mismatch norm |
|---|---:|---:|---:|---:|---:|---:|
| legacy `tangent` | 50 | 3.75268 | 7.00578 | -0.308147 | 213602 | 358831 |
| archived analytical branch | 50 | 0.0182187 | 0.638459 | 0.995489 | 563.047 | 22422.7 |

Top rows are still:

```text
element_name = LDPMTetra
material_id = 0
material_name = CSL material
material_status_names = CSL mat. statusx12
```

Interpretation: the archived analytical CSL tangent fixes the element-level CSL mismatch at the benchmark hard state, but the global Newton-direction check is still poor. The remaining global mismatch is no longer explained by the top CSL element-local tangent rows alone.

## Required Strict 8-Quarter Phase Closure

Saved baseline target from `intro.md`:

```text
converged steps: 8/8
total nonlinear iterations: 585
total OAS duration: 26:53
settings: tangent, stiffness_matrix_iter_update 10
```

The variant runs below used the temporary analytical `tangent` branch. From the archived checkpoint, the equivalent explicit setting is `stiff_matrix_type archived_csl_damage_tangent`.

### Analytical Tangent, Step-Start Matrix Only

Run directory:

```text
results/tsn65-csl-analytic-tangent-20260519-230004/runs/phase8_analytic_tangent_stepstart/
```

Settings:

```text
nonlinear_tangent_check 0
stiff_matrix_type tangent
stiffness_matrix_iter_update -1
```

Result:

| step | nonlinear rows | result |
|---:|---:|---|
| 1 | 8 | converged |
| 2 | 9 | converged |
| 3 | 12 | converged |
| 4 | 14 | converged |
| 5 | 18 | converged |
| 6 | 300 | failed max iterations |

Total OAS duration:

```text
00:10:56.446
```

Failure:

```text
Error: SteadyStateNonLinearSolver did not converge to the solution
```

Comparison to baseline:

```text
baseline: 8/8 converged, 585 iterations, 26:53
variant:  5/8 converged, failed at step 6, 361 nonlinear rows before failure, 10:56
```

### Analytical Tangent, Step-Start Matrix, Actual-State Backtracking

Run directory:

```text
results/tsn65-csl-analytic-tangent-20260519-230004/runs/phase8_analytic_tangent_stepstart_linesearch/
```

Added settings:

```text
nonlinear_material_snapshot_rollback 1
nonlinear_line_search backtracking
nonlinear_line_search_merit mixed
nonlinear_line_search_reduction 0.5
nonlinear_line_search_min_alpha 0.015625
nonlinear_line_search_max_trials 8
nonlinear_line_search_accept_any_decrease 1
nonlinear_line_search_evaluation actual
nonlinear_line_search_cutback_on_fail 0
```

Result:

| step | nonlinear rows | result |
|---:|---:|---|
| 1 | 8 | converged |
| 2 | 9 | converged |
| 3 | 12 | converged |
| 4 | 14 | converged |
| 5 | 18 | converged |
| 6 | 54 | converged with repeated `alpha=0.5` |
| 7 | 76 | stalled; manually stopped |

The run was stopped after about 23 minutes with step 7 stuck near:

```text
iteration 75 residual 5.843750e-03 displacement 3.003728e-04 energy 5.399725e-03
```

Comparison to baseline:

```text
baseline: 8/8 converged, 585 iterations, 26:53
variant:  6/8 converged, stalled at step 7, no valid 8-step phase closure
```

## Conclusion

The CSL active-damage tangent defect is real and the analytical derivative fixes it locally:

- local CSL active damage-growth `archived_csl_damage_tangent` matches finite differences,
- TS-N65 element-local top-row mismatch drops by roughly two orders of magnitude,
- the top rows remain `LDPMTetra` elements backed by `CSL material`.

However, the strict TS-N65 8-quarter phase closure does not improve:

- exact mid-step tangent rebuild destabilizes step 3,
- step-start analytical tangent fails at step 6,
- state-safe backtracking gets through step 6 but stalls in step 7.

This means the current blocker has moved from "CSL tangent is locally wrong" to "the locally correct softening tangent is not enough for strict load-control phase closure." The analytical branch is archived for reproduction, not enabled as default behavior. The next solver experiment should not claim improvement until it finishes the full `8/8` phase run. Good next candidates are path-following/indirect control, a controlled positive-definite tangent blend/limiter, or a model-level softening regularization/substepping strategy.
