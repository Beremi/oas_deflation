# TS-N_65 Two-Step DFGMRES/HYPRE j16 vs Pardiso Artifact

Generated: 2026-05-07

This report records the current 16-thread `DeflatedFGMRES` + hypre BoomerAMG-preconditioned run for the first two TS-N_65 loading steps and compares it with the saved 16-thread PardisoLDLT artifact. Pardiso was not rerun.

## Artifacts

| item | path |
| --- | --- |
| current DFGMRES/HYPRE run | `results/tsn65-dfgmres-hypre-cheby3-sweeps3-twostep-20260507-j16/dfgmres-hypre-tol1em1-N20` |
| current summary | `results/tsn65-dfgmres-hypre-cheby3-sweeps3-twostep-20260507-j16/summary.md` |
| current full runtime phase timings | `results/tsn65-dfgmres-hypre-cheby3-sweeps3-twostep-20260507-j16/dfgmres-hypre-tol1em1-N20/runtime_profile_summary.tsv` |
| current linear event timings | `results/tsn65-dfgmres-hypre-cheby3-sweeps3-twostep-20260507-j16/dfgmres-hypre-tol1em1-N20/linear_profile_events.tsv` |
| Pardiso 16-thread artifact | `results/tsn65-two-step-comparison-20260503-093031/pardisoldlt-reference` |

The current run used:

```text
threads = 16
total_time = 1.000000e-02
dfgmres_tolerance = 1e-1
dfgmres_true_tolerance = 1e-1
dfgmres_deflation_vectors = 20
hypre_relax_type = 16
hypre_cheby_order = 3
hypre_num_sweeps = 3
hypre_boomer_max_iterations = 1
dfgmres_elastic_reorder = 2
```

The OAS log confirms `hypre_threads=16`, `hypre_matrix_threads=16`, `active_elastic_reorder=coordinate_node_major`, `active_permutation=active`, `matvec_backend=librsb`, and `librsb_threads=8`.

## Main Comparison

`current/Pardiso` below is a time ratio for timing rows, so values below `1.0` mean the current run is faster. The Pardiso artifact was produced by an older code version, so the wall-time comparison includes both linear-solver differences and later OpenMP force/matrix changes.

| metric | current DFGMRES/HYPRE j16 | Pardiso j16 artifact | current/Pardiso |
| --- | ---: | ---: | ---: |
| wall s | 1200.121 | 4596.646 | 0.261 |
| OAS elapsed s | 1196.362 | 4589.158 | 0.261 |
| linear total s | 255.984 | 786.674 | 0.325 |
| analyze/setup/factor s | 11.259 | 55.649 | 0.202 |
| solve s | 244.725 | 731.026 | 0.335 |
| other/nonlinear s | 944.137 | 3802.484 | 0.248 |
| nonlinear iterations | 961 | 945 | 1.017 |
| linear solves | 961 | 945 | 1.017 |
| factorizations | 3 | 3 | 1.000 |
| outer Krylov iterations | 434 | 0 |  |
| median linear iterations | 0 | 0 |  |
| max linear iterations | 5 | 0 |  |

## Per-Step Timing

| solver | step | step wall s | linear total s | setup/factor s | solve s | solves | outer iters | max iter |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| current DFGMRES/HYPRE | 1 | 115.180 | 20.247 | 3.029 | 17.219 | 19 | 37 | 5 |
| current DFGMRES/HYPRE | 2 | 1081.180 | 230.786 | 3.280 | 227.506 | 942 | 397 | 2 |
| Pardiso artifact | 1 | 228.823 | 31.871 | 17.684 | 14.187 | 18 | 0 | 0 |
| Pardiso artifact | 2 | 4360.335 | 734.439 | 17.600 | 716.839 | 927 | 0 | 0 |

## Current DFGMRES Linear Internals

| phase | detail | count | total s | mean s | max s |
| --- | --- | ---: | ---: | ---: | ---: |
| linear.solve_total | direct_residual | 961 | 244.725 | 0.255 | 2.127 |
| linear.dfgmres.preconditioner_apply | direct_residual | 434 | 174.461 | 0.402 | 0.449 |
| linear.deflation_basis_orthogonalization | direct_residual | 961 | 47.686 | 0.050 | 0.064 |
| linear.dfgmres.solve_other | direct_residual | 961 | 41.936 | 0.044 | 0.202 |
| linear.dfgmres.matvec | direct_residual | 1807 | 23.544 | 0.013 | 0.017 |
| linear.preconditioner_setup | DeflatedFGMRES | 3 | 5.755 | 1.918 | 3.162 |
| linear.factorize_other | DeflatedFGMRES | 3 | 5.209 | 1.736 | 1.789 |
| linear.dfgmres.deflation_projection | direct_residual | 1394 | 4.637 | 0.003 | 0.014 |
| linear.deflation_basis_reorthogonalization | DeflatedFGMRES | 1 | 0.295 | 0.295 | 0.295 |
| linear.dfgmres.krylov_orthogonalization | direct_residual | 434 | 0.138 | 3.185e-04 | 0.004 |
| linear.dfgmres.least_squares | direct_residual | 434 | 0.009 | 2.168e-05 | 8.847e-05 |

## Current Runtime Hotspots

The full 83-row runtime timing table is in `runtime_profile_summary.tsv`. Largest rows:

| phase | detail | count | total s | mean s | max s |
| --- | --- | ---: | ---: | ---: | ---: |
| model.solve_total | - | 1 | 1122.539 | 1122.539 | 1122.539 |
| model.solve.solve_step | - | 2 | 1122.534 | 561.267 | 1081.180 |
| forces.total | active | 964 | 808.570 | 0.839 | 0.943 |
| forces.integrate_internal_forces | active | 964 | 704.731 | 0.731 | 0.832 |
| forces.integrate.strain_stress | active | 500960916 | 384.920 | 7.684e-07 | 8.879e-07 |
| linear.solve_total | direct_residual | 961 | 244.725 | 0.255 | 2.127 |
| forces.integrate.element_force_scatter | active | 500960916 | 227.513 | 4.542e-07 | 5.356e-07 |
| linear.dfgmres.preconditioner_apply | direct_residual | 434 | 174.461 | 0.402 | 0.449 |
| forces.simplex_volumetric_strains | active | 964 | 98.916 | 0.103 | 0.129 |
| forces.integrate.reset_eigenstrain | active | 964 | 86.528 | 0.090 | 0.108 |
| model.init_total | - | 1 | 65.479 | 65.479 | 65.479 |
| linear.deflation_basis_orthogonalization | direct_residual | 961 | 47.686 | 0.050 | 0.064 |
| linear.dfgmres.solve_other | direct_residual | 961 | 41.936 | 0.044 | 0.202 |
| model.init.elements | - | 1 | 32.536 | 32.536 | 32.536 |
| linear.dfgmres.matvec | direct_residual | 1807 | 23.544 | 0.013 | 0.017 |
| model.init.solver | - | 1 | 20.479 | 20.479 | 20.479 |
| matrix.factorize_linear_system_total | - | 3 | 13.255 | 4.418 | 6.018 |
| model.init.preprocessing_blocks | - | 1 | 11.908 | 11.908 | 11.908 |
| matrix.stiffness_update | secant | 3 | 10.697 | 3.566 | 8.500 |
| model.read_input | - | 1 | 8.344 | 8.344 | 8.344 |

## Previous Corrected DFGMRES Context

Compared with the earlier corrected two-step DFGMRES/HYPRE artifact at `results/tsn65-openmp-corrected-two-step-20260507/dfgmres-hypre-tol1em1-N20`, the current HYPRE-threaded Chebyshev/sweeps configuration is slightly faster in wall time, has similar linear time, and preserves the same order of Krylov work.

| metric | previous corrected DFGMRES j16 | current DFGMRES/HYPRE j16 |
| --- | ---: | ---: |
| wall s | 1222.163 | 1200.121 |
| linear total s | 252.407 | 255.984 |
| solve s | 234.834 | 244.725 |
| setup/analyze s | 17.574 | 11.259 |
| preconditioner apply s | 151.034 | 174.461 |
| matvec s | 38.294 | 23.544 |
| nonlinear iterations | 964 | 961 |
| linear solves | 964 | 961 |
| outer iterations | 438 | 434 |
| max linear iterations | 4 | 5 |
| final true residual | 0.007050 | 0.009181 |
| final nonlinear residual | 9.929e-04 | 9.995e-04 |
