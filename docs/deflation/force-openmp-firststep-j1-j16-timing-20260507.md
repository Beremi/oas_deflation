# TS-N_65 First-Step Runtime Dissemination: j=1 vs j=16

Date: 2026-05-07

Artifacts compared:

- `j=1`: `results/tsn65-force-scaling-firststep-current-20260507-j1/dfgmres-hypre-tol1em1-N20`
- `j=16`: `results/tsn65-force-scaling-firststep-current-20260507-j16/dfgmres-hypre-tol1em1-N20`

Config: `dfgmres-hypre-tol1em1-N20`, `total_time=0.005`, `tol=1e-1`, `N=20`. Both runs completed the first loading step with identical nonlinear and linear iteration counts.

## Run-Level Summary

| metric | j=1 | j=16 | speedup j1/j16 |
| --- | --- | --- | --- |
| wall_seconds | 204.999 | 122.318 | 1.67595 |
| solver elapsed step 1 | 201.356 | 118.528 | 1.6988 |
| nonlinear_iterations | 19 | 19 | 1 |
| linear_solves | 19 | 19 | 1 |
| DFGMRES outer_iterations | 37 | 37 | 1 |
| max linear iterations | 4 | 4 | 1 |
| final true relres | 0.0622639 | 0.0622639 | 1 |
| final nonlinear residual | 7.876e-04 | 7.876e-04 | 1 |

## Runtime Phase Summary

Sorted by `j=16` total time. Speedup is `j=1 total / j=16 total`; values near 1 are serial or weakly scaling bottlenecks.

Accounting note:

- Literal sum of all `j=16` runtime-profile rows below: `65.23 s`.
- Solver elapsed time for step 1: `118.53 s`; runner wall time: `122.32 s`.
- Literal row-sum missing time: `53.30 s` vs solver elapsed, or `57.09 s` vs runner wall.
- The literal row sum is not a clean non-overlapping total because `forces.total:*` includes force subphases, and `forces.step_end:converged` includes the frozen force pass. A cleaner non-overlap estimate keeps `forces.total:*`, excludes force subphases and `forces.step_end:converged`, and adds the linear solve time from `run.json`: `59.82 s` accounted.
- By that cleaner estimate, missing/unprofiled time is about `58.70 s` vs solver elapsed, or `62.49 s` vs runner wall.

| phase:detail | j=1 count | j=1 total s | j=1 mean s | j=16 count | j=16 total s | j=16 mean s | speedup | classification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| forces.total:active | 21 | 82.8224 | 3.94392 | 21 | 17.3993 | 0.828539 | 4.76009 | partial |
| forces.integrate_internal_forces:active | 21 | 80.5082 | 3.83373 | 21 | 15.1516 | 0.721503 | 5.31353 | partial |
| matrix.factorize_linear_system_total:- | 2 | 13.5193 | 6.75963 | 2 | 13.2268 | 6.61338 | 1.02211 | serial/weak |
| matrix.stiffness_update:secant | 2 | 18.5701 | 9.28504 | 2 | 9.81607 | 4.90803 | 1.89181 | weak |
| matrix.constraint_transform:- | 2 | 2.58281 | 1.2914 | 2 | 2.50722 | 1.25361 | 1.03015 | serial/weak |
| forces.simplex_volumetric_strains:active | 21 | 2.19647 | 0.104594 | 21 | 2.1438 | 0.102086 | 1.02457 | serial/weak |
| forces.total:frozen | 2 | 6.80956 | 3.40478 | 2 | 1.66092 | 0.830461 | 4.09987 | partial |
| forces.integrate_internal_forces:frozen | 2 | 6.58058 | 3.29029 | 2 | 1.44665 | 0.723325 | 4.54885 | partial |
| forces.step_end:converged | 1 | 3.89773 | 3.89773 | 1 | 0.820393 | 0.820393 | 4.75105 | partial |
| step.run_after_each_step:- | 1 | 0.548295 | 0.548295 | 1 | 0.235965 | 0.235965 | 2.32363 | weak |
| step_finalize.material_status_update:- | 1 | 0.546278 | 0.546278 | 1 | 0.233424 | 0.233424 | 2.34029 | weak |
| forces.simplex_volumetric_strains:frozen | 2 | 0.217687 | 0.108844 | 2 | 0.198443 | 0.0992216 | 1.09697 | serial/weak |
| forces.external_reactions:active | 21 | 0.108698 | 0.00517611 | 21 | 0.0964187 | 0.00459137 | 1.12736 | serial/weak |
| matrix.compute_keff_copy:- | 2 | 0.0849292 | 0.0424646 | 2 | 0.0901122 | 0.0450561 | 0.942483 | slower at 16t |
| field_update.expand_reduced_dofs:- | 20 | 0.0390362 | 0.00195181 | 20 | 0.0390771 | 0.00195385 | 0.998955 | slower at 16t |
| residual.reduce_force_array:- | 19 | 0.0367733 | 0.00193543 | 19 | 0.036625 | 0.00192763 | 1.00405 | serial/weak |
| error.energy_update:- | 19 | 0.0319188 | 0.00167994 | 19 | 0.0312751 | 0.00164606 | 1.02058 | serial/weak |
| energy.accumulate_internal_external:- | 20 | 0.0301101 | 0.00150551 | 20 | 0.029644 | 0.0014822 | 1.01572 | serial/weak |
| error.accumulate_norms:- | 19 | 0.0295422 | 0.00155485 | 19 | 0.0275484 | 0.00144991 | 1.07237 | serial/weak |
| forces.external_reactions:frozen | 2 | 0.01033 | 0.00516501 | 2 | 0.0147298 | 0.00736491 | 0.7013 | slower at 16t |
| forces.residual_vector:active | 21 | 0.00684313 | 3.259e-04 | 21 | 0.00445202 | 2.120e-04 | 1.53708 | weak |
| field_update.apply_increment:- | 20 | 0.00404327 | 2.022e-04 | 20 | 0.00413949 | 2.070e-04 | 0.976755 | slower at 16t |
| energy.physical_field_lookup:- | 20 | 0.00299281 | 1.496e-04 | 20 | 0.00289396 | 1.447e-04 | 1.03416 | serial/weak |
| field_update.dependency_update:- | 20 | 0.00263663 | 1.318e-04 | 20 | 0.00238523 | 1.193e-04 | 1.1054 | serial/weak |
| step_finalize.energy:- | 1 | 0.00164131 | 0.00164131 | 1 | 0.0018627 | 0.0018627 | 0.881144 | slower at 16t |
| error.prepare_buffers:- | 19 | 0.00161841 | 8.518e-05 | 19 | 0.00161116 | 8.480e-05 | 1.0045 | serial/weak |
| load.direct_step_setup:- | 1 | 0.00159329 | 0.00159329 | 1 | 0.00119844 | 0.00119844 | 1.32947 | serial/weak |
| forces.residual_vector:frozen | 2 | 7.672e-04 | 3.836e-04 | 2 | 8.424e-04 | 4.212e-04 | 0.910773 | slower at 16t |
| export.iteration_data:- | 19 | 6.311e-04 | 3.322e-05 | 19 | 6.768e-04 | 3.562e-05 | 0.932612 | slower at 16t |
| step_finalize.state_copy:- | 1 | 3.040e-04 | 3.040e-04 | 1 | 5.844e-04 | 5.844e-04 | 0.520154 | slower at 16t |
| step.run_before_each_step:- | 1 | 6.537e-04 | 6.537e-04 | 1 | 5.594e-04 | 5.594e-04 | 1.16856 | serial/weak |
| energy.kinetic:- | 20 | 8.907e-06 | 4.454e-07 | 20 | 7.140e-06 | 3.570e-07 | 1.24748 | serial/weak |
| error.finalize:- | 19 | 5.847e-06 | 3.077e-07 | 19 | 4.761e-06 | 2.506e-07 | 1.2281 | serial/weak |
| step_finalize.perturbations:- | 1 | 2.200e-07 | 2.200e-07 | 1 | 5.720e-07 | 5.720e-07 | 0.384615 | slower at 16t |

## Linear Solver Subtimings

These are measured by the linear profiler/run manifest, not by `runtime_profile_summary.tsv`. Hypre preconditioner apply is intentionally scoped to one OpenMP thread in the corrected version.

| linear part | j=1 total s | j=16 total s | speedup | classification |
| --- | --- | --- | --- | --- |
| linear setup/analyze total | 11.9379 | 11.7005 | 1.02029 | serial/weak |
| linear solve total | 15.8399 | 15.2949 | 1.03563 | serial/weak |
| preconditioner apply | 12.8459 | 12.3867 | 1.03708 | serial/weak |
| matvec | 1.60796 | 1.5537 | 1.03492 | serial/weak |
| orthogonalization | 0.0224785 | 0.0216129 | 1.04005 | serial/weak |
| deflation | 0.103611 | 0.0954242 | 1.0858 | serial/weak |
| least squares | 8.769e-04 | 8.291e-04 | 1.05763 | serial/weak |

## Bottleneck Readout

- `forces.integrate_internal_forces:active` is the largest OpenMP-scaled part: 80.51 s -> 15.15 s, 5.31x. This is useful but only about 33% efficiency on 16 threads.
- The biggest j=16 residual costs are now `forces.total:active`, `forces.integrate_internal_forces:active`, linear solve/setup, stiffness update, and matrix transform.
- `forces.simplex_volumetric_strains:active`, `matrix.factorize_linear_system_total`, `matrix.constraint_transform`, and linear preconditioner apply are serial or weakly scaling bottlenecks in this first-step sample.
- Whole-step speedup is only 1.70x because the force kernel is not the only large phase and several important phases remain serial/weakly threaded.
