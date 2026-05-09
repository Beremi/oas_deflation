# libRSB SpMV Scaling Diagnostics

## Scope

This diagnostic isolates pure SpMV on the dumped TS-N_65 replay matrix:

- Replay dump: `results/amgcl-tsn65-first-solve-20260503-201037/002-hybrid-aggr-blocksize1/linear_replay/solve_1`
- Matrix: `567923 x 567923`, `49429959` nonzeros
- Vector: `solution.tsv`
- Timed operation: repeated `y = A*x`, `beta = 0`, separate `x/y`, setup excluded
- Validation: every backend compared against Eigen row-major SpMV with relative tolerance `1e-10`
- Sweep artifact: `results/librsb-spmv-scaling-diagnostics-20260507-tsn65`
- Targeted autotune artifact: `results/librsb-spmv-scaling-diagnostics-20260507-tsn65-tune-targeted`

The production DFGMRES path was not changed by this pass.

## Main Results

All 100 untuned sweep rows passed validation. The largest observed relative error was around `1.1e-14`, so the tested backends are numerically consistent with Eigen for this matrix/vector.

| Case | Threads | Binding | Leaves | Time / SpMV | Speedup vs own 1-thread | Est. bandwidth |
|---|---:|---|---:|---:|---:|---:|
| Best overall, libRSB noflags subdivision 2 | 8 | close | 226 | 0.012127 s | 2.07x | 76.44 GiB/s |
| Best 16-thread, libRSB noflags subdivision 1 | 16 | close | 226 | 0.012972 s | 1.73x | 71.46 GiB/s |
| Eigen row-major | 16 | close | n/a | 0.013557 s | 1.44x | 68.38 GiB/s |
| Manual CRS OpenMP | 16 | close | n/a | 0.013305 s | 1.49x | 69.68 GiB/s |

The suspicious `1 -> 16` libRSB speedup is reproducible as a kernel-level limit, not an OAS timing artifact. libRSB reported the requested executing-thread count for every thread setting, including `16`.

## Threading And Partitioning

For the current production-like construction (`noflags`, subdivision `1`), libRSB reported:

| Threads | Actual RSB threads | Leaves | Time / SpMV, close |
|---:|---:|---:|---:|
| 1 | 1 | 16 | 0.022491 s |
| 2 | 2 | 34 | 0.015228 s |
| 4 | 4 | 52 | 0.013751 s |
| 8 | 8 | 160 | 0.013222 s |
| 16 | 16 | 226 | 0.012972 s |

Increasing the subdivision multiplier created more leaves, but did not monotonically improve time:

| Threads | Subdivision | Leaves | Time / SpMV, close noflags |
|---:|---:|---:|---:|
| 16 | 0.5 | 160 | 0.013252 s |
| 16 | 1 | 226 | 0.012972 s |
| 16 | 2 | 601 | 0.013943 s |
| 16 | 4 | 911 | 0.014898 s |

So the problem is not simply "too few leaves". More aggressive partitioning increases overhead or hurts locality once the matrix is already decomposed enough for this memory-bound operation.

## Autotune Check

A targeted `rsb_tune_spmm` run was tested for 8 and 16 threads, close/spread binding, subdivisions `1` and `2`.

- Best tuned 8-thread row: `0.012103 s`, close, subdivision `2`
- Best untuned 8-thread row: `0.012127 s`, close, subdivision `2`
- Best tuned 16-thread row: `0.012988 s`, close, subdivision `2`
- Best untuned 16-thread row: `0.012972 s`, close, subdivision `1`

Autotuning did not produce a meaningful 16-thread improvement. It should stay out of the production path.

## Interpretation

The low 16-thread speedup is not caused by librsb silently using one thread, because `RSB_IO_WANT_EXECUTING_THREADS` matched the requested count. It is also not caused by setup/tuning being included in the timed region; the timed loop contains only SpMV calls.

The likely limiting factor is memory bandwidth and locality. All backends flatten around `68-76 GiB/s`, and both Eigen row-major and manual CRS OpenMP show similarly poor scaling from 1 to 16 threads. On this matrix, 8 threads is usually the practical sweet spot for libRSB SpMV.

`RSB_FLAG_NOFLAGS` and `RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS` both produced matrix flags `0x2046186` in this build, and there was no clear performance difference between them.

## Verification

Build:

```bash
PKG_CONFIG_PATH=/home/beremi/local/librsb-1.3.0.2-native-openmp/lib/pkgconfig \
cmake --build ../oas_deflation-build/release --target OAS -j16
```

Regression tests:

```bash
OMP_NUM_THREADS=1  ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure
OMP_NUM_THREADS=16 ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure
```

Both test runs passed, including `2DPatchTestMechanics_DFGMRES_parallel_16`, which enables `amgcl_check_matrix 1`.

## Recommendation

Do not change production DFGMRES settings based on this diagnostic. The current cap around 8 librsb threads remains defensible: the pure SpMV benchmark shows 8 threads can be faster than 16 on the actual TS_N65 matrix.

The next useful performance work is outside single-vector SpMV: either reduce the number of DFGMRES vector operations and projections, batch operations where possible, or continue optimizing the force integration path, which dominates the first-step runtime more than SpMV.
