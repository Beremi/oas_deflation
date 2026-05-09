# Granular Runtime Comparison

Runs:
- `j=1`: `results/tsn65-granular-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20`
- `j=16`: `results/tsn65-granular-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

Both runs use the same first-step TS-N_65 DFGMRES/Hypre configuration: `total_time=5e-3`, `linear_tol=1e-1`, `true_tol=1e-1`, `N=20`.

The runtime profile summary `count` is now sample-aware. Normal phases count the number of events; DFGMRES subphases count actual inner actions such as matvecs and preconditioner applications; element assembly phases count element samples.

## Run Metrics

| Metric | j=1 | j=16 | j1/j16 |
| --- | --- | --- | --- |
| status | partial | partial |  |
| wall seconds | 197.338947 | 121.667372 | 1.62 |
| steps | 1 | 1 |  |
| nonlinear iterations | 19.000000 | 19.000000 | 1.00 |
| linear solves | 19.000000 | 19.000000 | 1.00 |
| DFGMRES outer iterations | 37.000000 | 37.000000 | 1.00 |
| final true residual | 0.062264 | 0.062264 | 1.00 |
| final nonlinear residual | 7.876e-04 | 7.876e-04 | 1.00 |

## Top-Level Accounting

These rows are non-overlapping. The `rest` row is the wall time not directly measured inside OAS.

| Bucket | j=1 s | j=1 share | j=16 s | j=16 share | speedup |
| --- | --- | --- | --- | --- | --- |
| `model.init_total` | 72.062 | 36.5% | 68.182 | 56.0% | 1.06 |
| `model.solve_total` | 113.252 | 57.4% | 41.366 | 34.0% | 2.74 |
| `model.read_input` | 8.583 | 4.3% | 8.540 | 7.0% | 1.01 |
| `rest` | 3.442 | 1.7% | 3.579 | 2.9% | 0.96 |
| `model.result_setup` | 1.244e-04 | 6.303e-05% | 1.801e-04 | 1.480e-04% | 0.69 |

## Force Assembly Detail

Rows are nested where appropriate and should not be summed unless explicitly noted by the phase name.

| Phase | Detail | j=1 count | j=1 total s | j=1 mean s | j=16 count | j=16 total s | j=16 mean s | speedup |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `forces.total` | `active` | 21 | 80.985 | 3.856448 | 21 | 17.208 | 0.819445 | 4.71 |
| `forces.integrate_internal_forces` | `active` | 21 | 78.754 | 3.750173 | 21 | 15.025 | 0.715486 | 5.24 |
| `forces.integrate.strain_stress` | `active` | 10,913,049 | 48.545 | 4.448e-06 | 10,913,049 | 8.295 | 7.601e-07 | 5.85 |
| `forces.integrate.element_force_scatter` | `active` | 10,913,049 | 28.379 | 2.600e-06 | 10,913,049 | 4.803 | 4.401e-07 | 5.91 |
| `forces.simplex_volumetric_strains` | `active` | 21 | 2.121 | 0.100976 | 21 | 2.080 | 0.099066 | 1.02 |
| `forces.integrate.reset_eigenstrain` | `active` | 21 | 1.826 | 0.086930 | 21 | 1.806 | 0.086018 | 1.01 |
| `forces.total` | `frozen` | 2 | 6.058 | 3.028812 | 2 | 1.528 | 0.764159 | 3.96 |
| `forces.integrate_internal_forces` | `frozen` | 2 | 5.848 | 2.924117 | 2 | 1.326 | 0.662844 | 4.41 |
| `forces.integrate.strain_stress` | `frozen` | 1,039,338 | 2.986 | 2.873e-06 | 1,039,338 | 0.708 | 6.814e-07 | 4.22 |
| `forces.integrate.element_force_scatter` | `frozen` | 1,039,338 | 2.692 | 2.590e-06 | 1,039,338 | 0.436 | 4.199e-07 | 6.17 |
| `forces.simplex_volumetric_strains` | `frozen` | 2 | 0.199 | 0.099611 | 2 | 0.193 | 0.096716 | 1.03 |
| `forces.integrate.reset_eigenstrain` | `frozen` | 2 | 0.170 | 0.084988 | 2 | 0.169 | 0.084584 | 1.00 |
| `forces.external_reactions` | `active` | 21 | 0.103 | 0.004915 | 21 | 0.096 | 0.004583 | 1.07 |
| `forces.integrate.thread_force_reduce` | `active` |  |  |  | 336 | 0.045 | 1.353e-04 |  |
| `forces.external_reactions` | `frozen` | 2 | 0.010 | 0.004752 | 2 | 0.009 | 0.004375 | 1.09 |
| `forces.residual_vector` | `active` | 21 | 0.007 | 3.248e-04 | 21 | 0.005 | 2.341e-04 | 1.39 |
| `forces.integrate.thread_force_reduce` | `frozen` |  |  |  | 32 | 0.004 | 1.391e-04 |  |
| `forces.integrate.zero_output` | `active` | 21 | 0.003 | 1.203e-04 | 21 | 0.002 | 9.445e-05 | 1.27 |
| `forces.residual_vector` | `frozen` | 2 | 5.549e-04 | 2.774e-04 | 2 | 3.211e-04 | 1.606e-04 | 1.73 |
| `forces.integrate.zero_output` | `frozen` | 2 | 1.946e-04 | 9.728e-05 | 2 | 1.770e-04 | 8.852e-05 | 1.10 |
| `forces.integrate.material_preparation` | `active` | 21 | 2.672e-05 | 1.272e-06 | 21 | 2.699e-05 | 1.285e-06 | 0.99 |
| `forces.integrate.material_preparation` | `frozen` | 2 | 2.385e-06 | 1.193e-06 | 2 | 2.034e-06 | 1.017e-06 | 1.17 |

## Matrix Assembly And Setup Detail

Rows are nested where appropriate and should not be summed unless explicitly noted by the phase name.

| Phase | Detail | j=1 count | j=1 total s | j=1 mean s | j=16 count | j=16 total s | j=16 mean s | speedup |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `matrix.factorize_linear_system_total` | `-` | 2 | 13.226 | 6.613187 | 2 | 13.038 | 6.518914 | 1.01 |
| `linear.preconditioner_setup` | `DeflatedFGMRES` | 2 | 9.848 | 4.923996 | 2 | 9.687 | 4.843466 | 1.02 |
| `matrix.stiffness_update` | `secant` | 2 | 17.217 | 8.608303 | 2 | 9.650 | 4.825128 | 1.78 |
| `matrix.assembly.entry_map` | `secant` |  |  |  | 1,039,338 | 7.443 | 7.161e-06 |  |
| `matrix.constraint_transform` | `-` | 2 | 2.498 | 1.248901 | 2 | 2.507 | 1.253585 | 1.00 |
| `linear.factorize_other` | `DeflatedFGMRES` | 2 | 1.849 | 0.924488 | 2 | 1.834 | 0.917191 | 1.01 |
| `matrix.prepare_pattern.compress` | `all` | 1 | 1.541 | 1.540728 | 1 | 1.542 | 1.542402 | 1.00 |
| `matrix.assembly.scatter` | `secant` | 1,039,338 | 10.716 | 1.031e-05 | 1,039,338 | 1.280 | 1.232e-06 | 8.37 |
| `matrix.prepare_pattern.element_triplets` | `all` | 519,669 | 1.130 | 2.175e-06 | 519,669 | 1.063 | 2.046e-06 | 1.06 |
| `matrix.assembly.local_matrix` | `secant` | 1,039,338 | 6.315 | 6.076e-06 | 1,039,338 | 0.867 | 8.337e-07 | 7.29 |
| `matrix.compute_keff_copy` | `-` | 2 | 0.081 | 0.040264 | 2 | 0.085 | 0.042492 | 0.95 |
| `matrix.assembly.zero_values` | `secant` | 2 | 0.148 | 0.074203 | 2 | 0.036 | 0.017912 | 4.14 |
| `matrix.prepare_pattern.clear_entry_map` | `all` | 1 | 5.710e-07 | 5.710e-07 | 1 | 7.410e-07 | 7.410e-07 | 0.77 |
| `matrix.assembly.lagrange_entries` | `secant` | 2 | 3.310e-07 | 1.655e-07 | 2 | 5.720e-07 | 2.860e-07 | 0.58 |
| `matrix.prepare_pattern.lagrange_triplets` | `stiffness` | 1 | 3.010e-07 | 3.010e-07 | 1 | 3.700e-07 | 3.700e-07 | 0.81 |

## Linear Solve Detail

Rows are nested where appropriate and should not be summed unless explicitly noted by the phase name.

| Phase | Detail | j=1 count | j=1 total s | j=1 mean s | j=16 count | j=16 total s | j=16 mean s | speedup |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `linear.solve_total` | `direct_residual` | 19 | 15.620 | 0.822129 | 19 | 15.408 | 0.810948 | 1.01 |
| `linear.dfgmres.preconditioner_apply` | `direct_residual` | 37 | 12.540 | 0.338930 | 37 | 12.435 | 0.336085 | 1.01 |
| `linear.preconditioner_setup` | `DeflatedFGMRES` | 2 | 9.848 | 4.923996 | 2 | 9.687 | 4.843466 | 1.02 |
| `linear.factorize_other` | `DeflatedFGMRES` | 2 | 1.849 | 0.924488 | 2 | 1.834 | 0.917191 | 1.01 |
| `linear.dfgmres.matvec` | `direct_residual` | 74 | 1.636 | 0.022111 | 74 | 1.607 | 0.021710 | 1.02 |
| `linear.dfgmres.solve_other` | `direct_residual` | 19 | 1.321 | 0.069516 | 19 | 1.244 | 0.065489 | 1.06 |
| `linear.deflation_basis_orthogonalization` | `direct_residual` | 19 | 0.782 | 0.041184 | 19 | 0.788 | 0.041484 | 0.99 |
| `linear.dfgmres.deflation_projection` | `direct_residual` | 55 | 0.100 | 0.001810 | 55 | 0.099 | 0.001797 | 1.01 |
| `linear.dfgmres.krylov_orthogonalization` | `direct_residual` | 37 | 0.023 | 6.107e-04 | 37 | 0.022 | 6.048e-04 | 1.01 |
| `linear.dfgmres.least_squares` | `direct_residual` | 37 | 8.296e-04 | 2.242e-05 | 37 | 8.432e-04 | 2.279e-05 | 0.98 |

## Full Runtime Profile

All measured runtime summary rows, sorted by `j=16 total s`. Nested rows are intentionally included for diagnosis, so this table should not be summed.

| Phase | Detail | j=1 count | j=1 total s | j=1 mean s | j=16 count | j=16 total s | j=16 mean s | speedup |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `model.init_total` | `-` | 1 | 72.062 | 72.061736 | 1 | 68.182 | 68.181744 | 1.06 |
| `model.solve_total` | `-` | 1 | 113.252 | 113.251772 | 1 | 41.366 | 41.365836 | 2.74 |
| `model.solve.solve_step` | `-` | 1 | 113.250 | 113.249588 | 1 | 41.364 | 41.363514 | 2.74 |
| `model.init.elements` | `-` | 1 | 32.426 | 32.426068 | 1 | 32.753 | 32.753072 | 0.99 |
| `model.init.solver` | `-` | 1 | 27.383 | 27.382841 | 1 | 22.140 | 22.140306 | 1.24 |
| `forces.total` | `active` | 21 | 80.985 | 3.856448 | 21 | 17.208 | 0.819445 | 4.71 |
| `linear.solve_total` | `direct_residual` | 19 | 15.620 | 0.822129 | 19 | 15.408 | 0.810948 | 1.01 |
| `forces.integrate_internal_forces` | `active` | 21 | 78.754 | 3.750173 | 21 | 15.025 | 0.715486 | 5.24 |
| `matrix.factorize_linear_system_total` | `-` | 2 | 13.226 | 6.613187 | 2 | 13.038 | 6.518914 | 1.01 |
| `model.init.preprocessing_blocks` | `-` | 1 | 11.710 | 11.709912 | 1 | 12.768 | 12.768285 | 0.92 |
| `linear.dfgmres.preconditioner_apply` | `direct_residual` | 37 | 12.540 | 0.338930 | 37 | 12.435 | 0.336085 | 1.01 |
| `linear.preconditioner_setup` | `DeflatedFGMRES` | 2 | 9.848 | 4.923996 | 2 | 9.687 | 4.843466 | 1.02 |
| `matrix.stiffness_update` | `secant` | 2 | 17.217 | 8.608303 | 2 | 9.650 | 4.825128 | 1.78 |
| `model.read_input` | `-` | 1 | 8.583 | 8.582998 | 1 | 8.540 | 8.540231 | 1.01 |
| `forces.integrate.strain_stress` | `active` | 10,913,049 | 48.545 | 4.448e-06 | 10,913,049 | 8.295 | 7.601e-07 | 5.85 |
| `matrix.assembly.entry_map` | `secant` |  |  |  | 1,039,338 | 7.443 | 7.161e-06 |  |
| `forces.integrate.element_force_scatter` | `active` | 10,913,049 | 28.379 | 2.600e-06 | 10,913,049 | 4.803 | 4.401e-07 | 5.91 |
| `matrix.constraint_transform` | `-` | 2 | 2.498 | 1.248901 | 2 | 2.507 | 1.253585 | 1.00 |
| `forces.simplex_volumetric_strains` | `active` | 21 | 2.121 | 0.100976 | 21 | 2.080 | 0.099066 | 1.02 |
| `linear.factorize_other` | `DeflatedFGMRES` | 2 | 1.849 | 0.924488 | 2 | 1.834 | 0.917191 | 1.01 |
| `forces.integrate.reset_eigenstrain` | `active` | 21 | 1.826 | 0.086930 | 21 | 1.806 | 0.086018 | 1.01 |
| `linear.dfgmres.matvec` | `direct_residual` | 74 | 1.636 | 0.022111 | 74 | 1.607 | 0.021710 | 1.02 |
| `matrix.prepare_pattern.compress` | `all` | 1 | 1.541 | 1.540728 | 1 | 1.542 | 1.542402 | 1.00 |
| `forces.total` | `frozen` | 2 | 6.058 | 3.028812 | 2 | 1.528 | 0.764159 | 3.96 |
| `forces.integrate_internal_forces` | `frozen` | 2 | 5.848 | 2.924117 | 2 | 1.326 | 0.662844 | 4.41 |
| `matrix.assembly.scatter` | `secant` | 1,039,338 | 10.716 | 1.031e-05 | 1,039,338 | 1.280 | 1.232e-06 | 8.37 |
| `linear.dfgmres.solve_other` | `direct_residual` | 19 | 1.321 | 0.069516 | 19 | 1.244 | 0.065489 | 1.06 |
| `matrix.prepare_pattern.element_triplets` | `all` | 519,669 | 1.130 | 2.175e-06 | 519,669 | 1.063 | 2.046e-06 | 1.06 |
| `matrix.assembly.local_matrix` | `secant` | 1,039,338 | 6.315 | 6.076e-06 | 1,039,338 | 0.867 | 8.337e-07 | 7.29 |
| `forces.step_end` | `converged` | 1 | 3.881 | 3.880784 | 1 | 0.806 | 0.806327 | 4.81 |
| `linear.deflation_basis_orthogonalization` | `direct_residual` | 19 | 0.782 | 0.041184 | 19 | 0.788 | 0.041484 | 0.99 |
| `forces.integrate.strain_stress` | `frozen` | 1,039,338 | 2.986 | 2.873e-06 | 1,039,338 | 0.708 | 6.814e-07 | 4.22 |
| `forces.integrate.element_force_scatter` | `frozen` | 1,039,338 | 2.692 | 2.590e-06 | 1,039,338 | 0.436 | 4.199e-07 | 6.17 |
| `model.init.exporters` | `-` | 1 | 0.395 | 0.394587 | 1 | 0.378 | 0.378166 | 1.04 |
| `step.run_after_each_step` | `-` | 1 | 0.489 | 0.488762 | 1 | 0.244 | 0.243832 | 2.00 |
| `step_finalize.material_status_update` | `-` | 1 | 0.487 | 0.486753 | 1 | 0.242 | 0.241787 | 2.01 |
| `forces.simplex_volumetric_strains` | `frozen` | 2 | 0.199 | 0.099611 | 2 | 0.193 | 0.096716 | 1.03 |
| `forces.integrate.reset_eigenstrain` | `frozen` | 2 | 0.170 | 0.084988 | 2 | 0.169 | 0.084584 | 1.00 |
| `model.init.nodes` | `-` | 1 | 0.133 | 0.133067 | 1 | 0.128 | 0.128202 | 1.04 |
| `linear.dfgmres.deflation_projection` | `direct_residual` | 55 | 0.100 | 0.001810 | 55 | 0.099 | 0.001797 | 1.01 |
| `forces.external_reactions` | `active` | 21 | 0.103 | 0.004915 | 21 | 0.096 | 0.004583 | 1.07 |
| `matrix.compute_keff_copy` | `-` | 2 | 0.081 | 0.040264 | 2 | 0.085 | 0.042492 | 0.95 |
| `forces.integrate.thread_force_reduce` | `active` |  |  |  | 336 | 0.045 | 1.353e-04 |  |
| `field_update.expand_reduced_dofs` | `-` | 20 | 0.039 | 0.001957 | 20 | 0.040 | 0.002001 | 0.98 |
| `residual.reduce_force_array` | `-` | 19 | 0.036 | 0.001870 | 19 | 0.036 | 0.001901 | 0.98 |
| `matrix.assembly.zero_values` | `secant` | 2 | 0.148 | 0.074203 | 2 | 0.036 | 0.017912 | 4.14 |
| `error.energy_update` | `-` | 19 | 0.032 | 0.001667 | 19 | 0.031 | 0.001645 | 1.01 |
| `error.accumulate_norms` | `-` | 19 | 0.031 | 0.001617 | 19 | 0.030 | 0.001602 | 1.01 |
| `energy.accumulate_internal_external` | `-` | 20 | 0.030 | 0.001510 | 20 | 0.030 | 0.001484 | 1.02 |
| `linear.dfgmres.krylov_orthogonalization` | `direct_residual` | 37 | 0.023 | 6.107e-04 | 37 | 0.022 | 6.048e-04 | 1.01 |
| `forces.external_reactions` | `frozen` | 2 | 0.010 | 0.004752 | 2 | 0.009 | 0.004375 | 1.09 |
| `model.init.constraints` | `-` | 1 | 0.009 | 0.008546 | 1 | 0.007 | 0.007363 | 1.16 |
| `model.init.fiber_assignment` | `-` | 1 | 0.007 | 0.006515 | 1 | 0.006 | 0.006163 | 1.06 |
| `forces.residual_vector` | `active` | 21 | 0.007 | 3.248e-04 | 21 | 0.005 | 2.341e-04 | 1.39 |
| `forces.integrate.thread_force_reduce` | `frozen` |  |  |  | 32 | 0.004 | 1.391e-04 |  |
| `field_update.apply_increment` | `-` | 20 | 0.004 | 2.037e-04 | 20 | 0.004 | 2.072e-04 | 0.98 |
| `energy.physical_field_lookup` | `-` | 20 | 0.003 | 1.366e-04 | 20 | 0.003 | 1.386e-04 | 0.99 |
| `field_update.dependency_update` | `-` | 20 | 0.002 | 1.198e-04 | 20 | 0.002 | 1.246e-04 | 0.96 |
| `model.solve.step_export` | `-` | 1 | 0.002 | 0.002161 | 1 | 0.002 | 0.002297 | 0.94 |
| `forces.integrate.zero_output` | `active` | 21 | 0.003 | 1.203e-04 | 21 | 0.002 | 9.445e-05 | 1.27 |
| `step_finalize.energy` | `-` | 1 | 0.002 | 0.001652 | 1 | 0.002 | 0.001660 | 1.00 |
| `load.direct_step_setup` | `-` | 1 | 0.002 | 0.001512 | 1 | 0.001 | 0.001359 | 1.11 |
| `error.prepare_buffers` | `-` | 19 | 0.001 | 7.044e-05 | 19 | 0.001 | 6.897e-05 | 1.02 |
| `linear.dfgmres.least_squares` | `direct_residual` | 37 | 8.296e-04 | 2.242e-05 | 37 | 8.432e-04 | 2.279e-05 | 0.98 |
| `export.iteration_data` | `-` | 19 | 5.596e-04 | 2.945e-05 | 19 | 5.586e-04 | 2.940e-05 | 1.00 |
| `step.run_before_each_step` | `-` | 1 | 5.086e-04 | 5.086e-04 | 1 | 5.356e-04 | 5.356e-04 | 0.95 |
| `step_finalize.state_copy` | `-` | 1 | 3.002e-04 | 3.002e-04 | 1 | 3.220e-04 | 3.220e-04 | 0.93 |
| `forces.residual_vector` | `frozen` | 2 | 5.549e-04 | 2.774e-04 | 2 | 3.211e-04 | 1.606e-04 | 1.73 |
| `model.result_setup` | `-` | 1 | 1.244e-04 | 1.244e-04 | 1 | 1.801e-04 | 1.801e-04 | 0.69 |
| `forces.integrate.zero_output` | `frozen` | 2 | 1.946e-04 | 9.728e-05 | 2 | 1.770e-04 | 8.852e-05 | 1.10 |
| `forces.integrate.material_preparation` | `active` | 21 | 2.672e-05 | 1.272e-06 | 21 | 2.699e-05 | 1.285e-06 | 0.99 |
| `energy.kinetic` | `-` | 20 | 8.395e-06 | 4.197e-07 | 20 | 7.494e-06 | 3.747e-07 | 1.12 |
| `error.finalize` | `-` | 19 | 5.721e-06 | 3.011e-07 | 19 | 5.238e-06 | 2.757e-07 | 1.09 |
| `model.solve.initial_export` | `-` | 1 | 4.198e-06 | 4.198e-06 | 1 | 4.668e-06 | 4.668e-06 | 0.90 |
| `model.init.materials` | `-` | 1 | 1.973e-06 | 1.973e-06 | 1 | 3.467e-06 | 3.467e-06 | 0.57 |
| `model.solve.print_step_time` | `-` | 1 | 1.973e-06 | 1.973e-06 | 1 | 2.245e-06 | 2.245e-06 | 0.88 |
| `forces.integrate.material_preparation` | `frozen` | 2 | 2.385e-06 | 1.193e-06 | 2 | 2.034e-06 | 1.017e-06 | 1.17 |
| `model.init.boundary_conditions` | `-` | 1 | 1.594e-06 | 1.594e-06 | 1 | 1.904e-06 | 1.904e-06 | 0.84 |
| `matrix.prepare_pattern.clear_entry_map` | `all` | 1 | 5.710e-07 | 5.710e-07 | 1 | 7.410e-07 | 7.410e-07 | 0.77 |
| `matrix.assembly.lagrange_entries` | `secant` | 2 | 3.310e-07 | 1.655e-07 | 2 | 5.720e-07 | 2.860e-07 | 0.58 |
| `matrix.prepare_pattern.lagrange_triplets` | `stiffness` | 1 | 3.010e-07 | 3.010e-07 | 1 | 3.700e-07 | 3.700e-07 | 0.81 |
| `step_finalize.perturbations` | `-` | 1 | 3.510e-07 | 3.510e-07 | 1 | 3.300e-07 | 3.300e-07 | 1.06 |
