# Runtime Accounting After Timing Instrumentation

Run:
`results/tsn65-force-accounting-firststep-timed-20260507-j16/dfgmres-hypre-tol1em1-N20`

Command:
`python3 scripts/run-tsn65-dfgmres-hypre-deflation.py --oas-bin ../oas_deflation-build/release/bin/OAS --threads 16 --total-time 5.000000e-03 --linear-tol 1e-1 --true-tol 1e-1 --nvecs 20 --only dfgmres-hypre-tol1em1-N20 --out-dir results/tsn65-force-accounting-firststep-timed-20260507-j16`

Status: partial first-step run, 1 loading step, 19 nonlinear iterations, 19 linear solves, 37 DFGMRES outer iterations.

## Top-Level Accounting

Total runtime is `run.json:wall_seconds = 122.670 s`. These rows are non-overlapping, so the sum including `rest` is the total runtime.

| Rank | Runtime bucket | Seconds | Share of total | Notes |
| --- | --- | ---: | ---: | --- |
| 1 | `model.init_total` | 68.775 | 56.1% | Full model initialization, including initial solver setup. |
| 2 | `model.solve_total` | 41.483 | 33.8% | Full solve loop, including the first loading step and post-step export. |
| 3 | `model.read_input` | 8.836 | 7.2% | Input parsing and model construction from input files. |
| 4 | `rest` | 3.576 | 2.9% | Process startup, final console reporting, destruction, and runner overhead not timed inside OAS. |
| 5 | `model.result_setup` | 0.000 | 0.0% | Result directory and `version.txt` setup. |
|  | **Total** | **122.670** | **100.0%** |  |

The missing time is now below the requested 5% threshold.

## Init Breakdown

These rows are nested inside `model.init_total`; they are not added again to the top-level total.

| Rank | Init phase | Seconds | Share of init |
| --- | --- | ---: | ---: |
| 1 | `model.init.elements` | 32.916 | 47.9% |
| 2 | `model.init.solver` | 22.670 | 33.0% |
| 3 | `model.init.preprocessing_blocks` | 12.616 | 18.3% |
| 4 | `model.init.exporters` | 0.430 | 0.6% |
| 5 | `model.init.nodes` | 0.128 | 0.2% |
| 6 | `model.init.constraints` | 0.007 | 0.0% |
| 7 | `model.init.fiber_assignment` | 0.006 | 0.0% |
| 8 | `model.init.rest` | 0.000 | 0.0% |
| 9 | `model.init.boundary_conditions` | 0.000 | 0.0% |
| 10 | `model.init.materials` | 0.000 | 0.0% |
|  | **Init total** | **68.775** | **100.0%** |

`model.init.solver` includes the initial stiffness assembly/factorization and initial force computations.

## First-Step Solve Detail

These rows are a non-overlapping accounting of `model.solve.solve_step = 41.481 s`. Force subphases are excluded here because `forces.total` already contains them. `forces.step_end` is also excluded to avoid double-counting the final converged force evaluation, which is already included in `forces.total active`.

| Rank | Solve-step bucket | Count | Seconds | Share of solve step |
| --- | --- | ---: | ---: | ---: |
| 1 | `forces.total active` | 20 | 16.417 | 39.6% |
| 2 | `linear.solve_total direct_residual` | 19 | 15.357 | 37.0% |
| 3 | `matrix.factorize_linear_system_total` | 1 | 5.330 | 12.8% |
| 4 | `matrix.constraint_transform` | 1 | 1.328 | 3.2% |
| 5 | `matrix.stiffness_update secant` | 1 | 1.069 | 2.6% |
| 6 | `forces.total frozen` | 1 | 0.794 | 1.9% |
| 7 | `solve_step.rest` | 1 | 0.779 | 1.9% |
| 8 | `step.run_after_each_step` | 1 | 0.227 | 0.5% |
| 9 | `field_update.expand_reduced_dofs` | 20 | 0.039 | 0.1% |
| 10 | `residual.reduce_force_array` | 19 | 0.036 | 0.1% |
| 11 | `matrix.compute_keff_copy` | 1 | 0.035 | 0.1% |
| 12 | `error.energy_update` | 19 | 0.031 | 0.1% |
| 13 | `error.accumulate_norms` | 19 | 0.029 | 0.1% |
| 14 | `field_update.apply_increment` | 20 | 0.004 | 0.0% |
| 15 | `field_update.dependency_update` | 20 | 0.003 | 0.0% |
| 16 | `error.prepare_buffers` | 19 | 0.001 | 0.0% |
| 17 | `load.direct_step_setup` | 1 | 0.001 | 0.0% |
| 18 | `export.iteration_data` | 19 | 0.001 | 0.0% |
| 19 | `step.run_before_each_step` | 1 | 0.000 | 0.0% |
| 20 | `error.finalize` | 19 | 0.000 | 0.0% |
|  | **Solve-step total** |  | **41.481** | **100.0%** |
