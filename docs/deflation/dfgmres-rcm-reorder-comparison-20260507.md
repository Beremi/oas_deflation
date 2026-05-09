# DFGMRES RCM Reordering Check, 2026-05-07

## Change

Added DFGMRES elastic reorder mode `3`, a node-block Reverse Cuthill-McKee ordering:

```text
active matrix pattern -> node graph -> RCM node order -> keep local DOFs contiguous
```

This keeps the HYPRE system-AMG-friendly node-major layout while testing a bandwidth-reducing graph ordering. It is different from mode `2`, which is coordinate-sorted node-major ordering.

Implementation:

- `dfgmres_elastic_reorder 3`
- `hypre_elastic_reorder 3` remains a compatibility path through the existing campaign script
- rows/columns are symmetrically permuted before diagonal scaling, CRS extraction, librsb construction, and HYPRE setup
- RHS, deflation vectors, and solution are mapped at the DFGMRES boundary

## Spy Images

Artifact: `results/tsn65-rcm-active-matrix-spy-20260507`

Source replay: `results/amgcl-tsn65-first-solve-20260503-201037/002-hybrid-aggr-blocksize1/linear_replay/solve_1`

| metric | value |
| --- | ---: |
| reduced rows | 567,923 |
| active rows | 587,190 |
| active block size | 6 |
| active nodes | 97,865 |
| nodes with matrix entries | 97,379 |
| isolated/empty nodes in node graph | 486 |
| nnz | 49,429,959 |
| active density | 1.433616e-04 |

![RCM overview](../../results/tsn65-rcm-active-matrix-spy-20260507/tsn65_rcm_active_spy_overview_log_density.png)

![RCM zoom 3000](../../results/tsn65-rcm-active-matrix-spy-20260507/tsn65_rcm_active_spy_zoom_first_3000.png)

![RCM zoom 12000](../../results/tsn65-rcm-active-matrix-spy-20260507/tsn65_rcm_active_spy_zoom_first_12000.png)

![RCM density](../../results/tsn65-rcm-active-matrix-spy-20260507/tsn65_rcm_active_nnz_density_by_index.png)

## Correctness

Build:

```bash
PKG_CONFIG_PATH=/home/beremi/local/librsb-1.3.0.2-native-openmp/lib/pkgconfig \
cmake --build ../oas_deflation-build/release --target OAS -j16
```

CTest:

| command | result |
| --- | --- |
| `OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |
| `OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |

TS-N_65 first-step RCM runs both completed with status `0`.

## First-Step Timing

Baseline is the previously saved global coordinate node-major reorder artifact. RCM was run fresh with the same first-step setup.

| ordering | threads | wall s | solve s | setup/analyze s | nonlinear iters | linear solves | outer iters | max iter | final true relres |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| coordinate node-major | 1 | 205.693 | 18.407 | 12.434 | 19 | 19 | 52 | 8 | 0.0866123 |
| RCM node-major | 1 | 205.709 | 16.474 | 13.522 | 19 | 19 | 53 | 8 | 0.0739647 |
| coordinate node-major | 16 | 114.340 | 12.772 | 7.796 | 19 | 19 | 52 | 8 | 0.0866123 |
| RCM node-major | 16 | 115.998 | 13.175 | 9.006 | 19 | 19 | 53 | 8 | 0.0739647 |

## Requested Kernel Comparison

| ordering | threads | HYPRE setup count | HYPRE setup s | HYPRE apply count | HYPRE apply s | matvec count | matvec s | factorize other s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| coordinate node-major | 1 | 2 | 8.671 | 52 | 15.031 | 90 | 1.829 | 3.763 |
| RCM node-major | 1 | 2 | 8.520 | 53 | 13.417 | 91 | 1.532 | 5.002 |
| coordinate node-major | 16 | 2 | 4.331 | 52 | 10.069 | 90 | 1.148 | 3.465 |
| RCM node-major | 16 | 2 | 4.277 | 53 | 10.477 | 91 | 1.179 | 4.729 |

Interpretation:

- Single-thread RCM improves the repeated solve kernels: matvec and HYPRE apply are faster, but setup/analyze grows because RCM graph construction and permuted matrix rebuild cost more.
- At 16 threads, RCM does not help the solve path: matvec and HYPRE apply are slightly slower, and setup/analyze is worse.
- HYPRE preconditioner setup alone is slightly faster with RCM, but that does not offset the extra non-HYPRE factorization work.
- RCM increases DFGMRES outer iterations by one (`52 -> 53`) while keeping the same max iteration count (`8`) and a better final true residual in this run.

## Artifacts

| artifact | path |
| --- | --- |
| RCM spy | `results/tsn65-rcm-active-matrix-spy-20260507` |
| RCM j=1 | `results/tsn65-dfgmres-hypre-rcm-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20` |
| RCM j=16 | `results/tsn65-dfgmres-hypre-rcm-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20` |
| coordinate baseline j=1 | `results/tsn65-dfgmres-hypre-global-reorder-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20` |
| coordinate baseline j=16 | `results/tsn65-dfgmres-hypre-global-reorder-firststep-20260507-j16-nobind/dfgmres-hypre-tol1em1-N20` |

## Decision

Do not switch the TS-N_65 campaign default from coordinate node-major to RCM based on this test. RCM is numerically safe in the first-step run, but the 16-thread path is not faster and it adds setup overhead.

## Original HYPRE Smoother With RCM

After the Chebyshev test above, the original guarded HYPRE smoother was also tested with RCM:

```text
hypre_relax_type 6
hypre_cheby_order 0
hypre_boomer_max_iterations 1
hypre_threads 0
dfgmres_elastic_reorder 3
```

This uses the current automatic HYPRE policy: HYPRE matrix insertion can use OpenMP, but BoomerAMG setup/apply is clamped to one HYPRE thread because relax type `6` changes behavior when run with OpenMP solve threads.

| ordering / smoother | threads | wall s | linear total s | solve s | setup/analyze s | HYPRE setup s | HYPRE apply s | matvec s | DFGMRES outer | max iter | final true relres |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| original HYPRE, no active reorder | 16 | 120.415 | 27.227 | 14.888 | 12.339 | 9.590 | 12.461 | 1.048 | 37 | 4 | 0.0622639 |
| original HYPRE, RCM node-major | 16 | 123.018 | 26.647 | 14.213 | 12.433 | 7.606 | 11.790 | 0.969 | 37 | 4 | 0.0880175 |
| Cheby4, coordinate node-major | 16 | 114.340 | 20.568 | 12.772 | 7.796 | 4.331 | 10.069 | 1.148 | 52 | 8 | 0.0866123 |

Interpretation:

- RCM with the original smoother keeps the original iteration count: `37` DFGMRES outer iterations and max `4`.
- RCM improves the repeated original-HYPRE kernels: setup, apply, and matvec all drop.
- Total linear time improves slightly (`27.227 s -> 26.647 s`), but wall time is worse in this one run because non-linear-solver phases outside the linear profile were slower.
- Cheby4 coordinate ordering is still faster in linear time and first-step wall time, but it gets there by accepting more Krylov work (`52` outer iterations instead of `37`).

Artifact:

- `results/tsn65-dfgmres-hypre-original-rcm-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`
