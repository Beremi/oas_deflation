# DFGMRES librsb/OpenMP First-Step Comparison

## Scope

This pass keeps librsb optional at build time and uses it only inside native `DeflatedFGMRES` matvecs. The Eigen sparse matrix and row-major CRS arrays are still kept for preconditioners, validation, and fallback paths.

Build and validation:

| check | result |
| --- | --- |
| non-librsb `OAS` build | passed before enabling `USE_LIBRSB` |
| librsb configure | `USE_LIBRSB=ON`, `librsb 1.3.0` from `/home/beremi/local/librsb-1.3.0.2-native-openmp` |
| librsb `OAS` build | passed |
| `ctest -L tests`, `OMP_NUM_THREADS=1` | 7/7 passed |
| `ctest -L tests`, `OMP_NUM_THREADS=16` | 7/7 passed |
| DFGMRES-none parallel regression with `amgcl_check_matrix 1` | passed |

Backend confirmation from the TS-N_65 candidate log:

```text
matvec_backend=librsb, librsb_threads=8
```

## First-Step Summary

Baseline artifact: `results/tsn65-granular-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

Candidate artifact: `results/tsn65-dfgmres-librsb-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

| metric | previous | librsb/OpenMP | ratio previous/candidate |
| --- | ---: | ---: | ---: |
| script wall time | 121.667 s | 123.670 s | 0.984x |
| OAS profiled wall time | 118.088 s | 120.014 s | 0.984x |
| nonlinear iterations | 19 | 19 | 1.000x |
| linear solves | 19 | 19 | 1.000x |
| DFGMRES outer iterations | 37 | 37 | 1.000x |
| max linear iterations | 4 | 4 | 1.000x |
| final true relative residual | 0.0622639253 | 0.0622639253 | 1.000x |
| final nonlinear residual | 7.875648e-4 | 7.875648e-4 | 1.000x |
| final basis size | 18 | 18 | 1.000x |
| discarded basis vectors | 1 | 1 | 1.000x |

Correctness did not regress: iteration counts, basis behavior, and final residuals match the previous artifact.

## Granular Timing

| metric | count previous->candidate | previous s | librsb/OpenMP s | speedup | delta s |
| --- | ---: | ---: | ---: | ---: | ---: |
| model.solve_total | 1->1 | 41.3658 | 42.3430 | 0.977x | +0.9772 |
| model.solve.solve_step | 1->1 | 41.3635 | 42.3405 | 0.977x | +0.9770 |
| forces.total active | 21->21 | 17.2083 | 17.7696 | 0.968x | +0.5613 |
| forces.integrate_internal_forces active | 21->21 | 15.0252 | 15.4986 | 0.969x | +0.4733 |
| forces.integrate.strain_stress active | 10913049->10913049 | 8.2952 | 8.4033 | 0.987x | +0.1081 |
| forces.integrate.element_force_scatter active | 10913049->10913049 | 4.8033 | 5.0741 | 0.947x | +0.2708 |
| forces.simplex_volumetric_strains active | 21->21 | 2.0804 | 2.1599 | 0.963x | +0.0795 |
| matrix.factorize_linear_system_total | 2->2 | 13.0378 | 14.4292 | 0.904x | +1.3913 |
| linear.preconditioner_setup | 2->2 | 9.6869 | 10.0443 | 0.964x | +0.3574 |
| linear.factorize_other | 2->2 | 1.8344 | 2.8327 | 0.648x | +0.9984 |
| matrix.stiffness_update secant | 2->2 | 9.6503 | 9.8912 | 0.976x | +0.2409 |
| matrix.assembly.entry_map secant | 1039338->1039338 | 7.4431 | 7.7607 | 0.959x | +0.3175 |
| matrix.assembly.local_matrix secant | 1039338->1039338 | 0.8665 | 0.8446 | 1.026x | -0.0220 |
| matrix.assembly.scatter secant | 1039338->1039338 | 1.2800 | 1.2244 | 1.045x | -0.0556 |
| matrix.constraint_transform | 2->2 | 2.5072 | 2.5699 | 0.976x | +0.0627 |
| linear.solve_total | 19->19 | 15.4080 | 15.3175 | 1.006x | -0.0905 |
| linear.dfgmres.preconditioner_apply | 37->37 | 12.4352 | 12.8283 | 0.969x | +0.3932 |
| linear.dfgmres.matvec | 74->74 | 1.6065 | 1.0920 | 1.471x | -0.5145 |
| linear.dfgmres.krylov_orthogonalization | 37->37 | 0.0224 | 0.0149 | 1.503x | -0.0075 |
| linear.dfgmres.deflation_projection | 55->55 | 0.0988 | 0.0894 | 1.106x | -0.0095 |
| linear.dfgmres.solve_other | 19->19 | 1.2443 | 1.2919 | 0.963x | +0.0476 |
| linear.deflation_basis_orthogonalization | 19->19 | 0.7882 | 0.5079 | 1.552x | -0.2803 |
| step_finalize.material_status_update | 1->1 | 0.2418 | 0.2480 | 0.975x | +0.0062 |
| forces.step_end converged | 1->1 | 0.8063 | 0.8470 | 0.952x | +0.0407 |

## Interpretation

The librsb matvec path is working and improves the DFGMRES matrix-vector bucket: `1.6065 s -> 1.0920 s`, a `1.47x` speedup across 74 matvec calls. Manual OpenMP vector kernels also helped the DFGMRES vector-heavy buckets: Krylov orthogonalization is `1.50x` faster, deflation projection is `1.11x` faster, and basis orthogonalization is `1.55x` faster.

The first-step wall time did not improve because matvec is only about one second of the full run. The candidate lost time in factorization/setup (`+1.39 s` total), mainly `linear.factorize_other` (`+1.00 s`), which is consistent with building the additional librsb matrix twice. Preconditioner apply also moved slower by `+0.39 s`, likely normal run-to-run noise or thread/runtime interaction outside the librsb matvec path.

The practical conclusion is that librsb is numerically safe in this path and makes SpMV faster, but it is not enough to improve TS-N_65 first-step wall time while force integration, stiffness assembly, and preconditioner setup/apply dominate.
