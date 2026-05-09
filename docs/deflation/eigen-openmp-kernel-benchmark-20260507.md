# Eigen OpenMP Kernel Microbenchmark, 2026-05-07

Purpose: isolate Eigen's native OpenMP behavior for dense vector operations,
dense dot products, and sparse matrix-vector products with RowMajor and
ColMajor storage. A second pass adds manual OpenMP reference kernels to verify
whether poor scaling is caused by OpenMP itself or by Eigen not parallelizing a
specific operation.

## Setup

- Machine: AMD Ryzen 9 9950X, 16 physical cores, 32 hardware threads.
- Compiler flags:
  `-std=c++17 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG -fopenmp`
- Eigen include: `../oas_deflation-build/release/_deps/eigen-src`
- Runtime binding:
  `OMP_PROC_BIND=close OMP_PLACES=cores`
- Benchmark source:
  `scripts/benchmark-eigen-openmp-kernels.cpp`
- Raw results, native Eigen only:
  `results/eigen-openmp-kernels-20260507/`
- Raw results, native Eigen plus manual OpenMP reference kernels:
  `results/eigen-openmp-kernels-20260507-v2/`
- Raw results, corrected metadata/thread-observation pass:
  `results/eigen-openmp-kernels-20260507-v3/`
- Raw results, CBLAS linked to OpenBLAS:
  `results/eigen-openmp-kernels-20260507-blas/`
- Raw results, CBLAS linked to Intel MKL:
  `results/eigen-openmp-kernels-20260507-mkl-blas/`
- Raw results, system OpenBLAS with explicit `openblas_set_num_threads()`:
  `results/eigen-openmp-kernels-20260507-system-openblas-config/`
- Raw results, local OpenBLAS `ZEN` OpenMP build:
  `results/eigen-openmp-kernels-20260507-local-openblas/`
- Raw results, local OpenBLAS `ZEN` pthread build:
  `results/eigen-openmp-kernels-20260507-local-openblas-pthread/`
- Raw results, local OpenBLAS `COOPERLAKE` OpenMP build:
  `results/eigen-openmp-kernels-20260507-local-openblas-cooperlake/`
- Raw results, local `librsb` SpMV build, narrow band:
  `results/eigen-openmp-kernels-20260507-librsb/`
- Raw results, local `librsb` SpMV build, wide band:
  `results/eigen-openmp-kernels-20260507-librsb-wideband/`
- Raw results, local `librsb` tuned SpMV build, narrow band:
  `results/eigen-openmp-kernels-20260507-librsb-tuned/`
- Raw results, local `librsb` tuned SpMV build, wide band:
  `results/eigen-openmp-kernels-20260507-librsb-wideband-tuned/`
- Raw results, actual TS-N_65 replay matrix/vectors:
  `results/eigen-openmp-kernels-20260507-replay-tsn65/`

Benchmark sizes:

- Dense vectors: 8,000,000 doubles.
- Vector add/sub/axpby repeats: 80.
- Dot repeats: 120.
- Sparse matrix: 1,000,000 x 1,000,000, band radius 3.
- Sparse nonzeros: 6,999,988.
- Sparse matvec repeats: 30.

The corrected binary reports both Eigen's requested thread count and the number
of threads observed inside an actual OpenMP parallel region. For the v3
`OMP_NUM_THREADS=16` run it reported:

```text
eigen_version=3.5.0 openmp_enabled=1 openmp_macro=201511
openmp_max_threads=16 openmp_observed_threads=16 eigen_threads=16
```

`Eigen::initParallel()` is not used in the benchmark because the vendored Eigen
3.5 headers mark it deprecated; thread control is done through
`Eigen::setNbThreads()` and the OpenMP runtime.

## Initial Timing Summary

Average of three runs.

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup |
|---|---:|---:|---:|---:|
| vector add | `VectorXd` | 0.453868 | 0.441048 | 1.03x |
| vector subtract | `VectorXd` | 0.439876 | 0.429765 | 1.02x |
| vector axpby | `VectorXd` | 0.466118 | 0.434281 | 1.07x |
| dot | `VectorXd` | 0.310233 | 0.291252 | 1.07x |
| sparse matvec | `RowMajor` | 0.090420 | 0.057718 | 1.57x |
| sparse matvec | `ColMajor` | 0.097868 | 0.098286 | 1.00x |

This first table was incomplete: it measured Eigen's native kernels only. It
did not prove that OpenMP itself was ineffective.

## Repaired Timing Summary

Average of three runs. Manual OpenMP rows are direct loop implementations using
`#pragma omp parallel for`, with static scheduling for vector loops and RowMajor
CRS matvec.

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup |
|---|---:|---:|---:|---:|
| vector add | `VectorXd` | 0.448659 | 0.429973 | 1.04x |
| vector add | `manual_openmp` | 0.439109 | 0.369730 | 1.19x |
| vector subtract | `VectorXd` | 0.444627 | 0.449708 | 0.99x |
| vector subtract | `manual_openmp` | 0.437184 | 0.364857 | 1.20x |
| vector axpby | `VectorXd` | 0.430099 | 0.443266 | 0.97x |
| vector axpby | `manual_openmp` | 0.436260 | 0.365185 | 1.19x |
| dot | `VectorXd` | 0.291157 | 0.300666 | 0.97x |
| dot | `manual_openmp` | 0.378213 | 0.183362 | 2.06x |
| sparse matvec | `RowMajor` | 0.090930 | 0.057965 | 1.57x |
| sparse matvec | `manual_rowmajor_openmp` | 0.091156 | 0.028146 | 3.24x |
| sparse matvec | `ColMajor` | 0.099585 | 0.098981 | 1.01x |

## Corrected Metadata Pass

Average of three runs from `results/eigen-openmp-kernels-20260507-v3/`. This
pass keeps the same kernels but adds explicit reporting that OpenMP was compiled
in and that a parallel region actually launched 16 threads.

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup | j=1 GiB/s | j=16 GiB/s |
|---|---:|---:|---:|---:|---:|---:|
| vector add | `VectorXd` | 0.443022 | 0.439513 | 1.01x | 32.3 | 32.5 |
| vector add | `manual_openmp` | 0.434147 | 0.375023 | 1.16x | 33.0 | 38.2 |
| vector subtract | `VectorXd` | 0.449786 | 0.444978 | 1.01x | 31.8 | 32.1 |
| vector subtract | `manual_openmp` | 0.450133 | 0.385753 | 1.17x | 31.8 | 37.2 |
| vector axpby | `VectorXd` | 0.440117 | 0.463387 | 0.95x | 32.5 | 30.9 |
| vector axpby | `manual_openmp` | 0.447587 | 0.374326 | 1.20x | 32.0 | 38.2 |
| dot | `VectorXd` | 0.291312 | 0.310926 | 0.94x | 49.1 | 46.0 |
| dot | `manual_openmp` | 0.388456 | 0.199160 | 1.95x | 36.8 | 72.1 |
| sparse matvec | `RowMajor` | 0.095484 | 0.062421 | 1.53x | 29.3 | 44.8 |
| sparse matvec | `manual_rowmajor_openmp` | 0.097372 | 0.030926 | 3.15x | 28.9 | 90.7 |
| sparse matvec | `ColMajor` | 0.101007 | 0.102363 | 0.99x | 27.7 | 27.3 |

## BLAS Level-1 Comparison

This pass adds direct CBLAS calls for vector addition and dot products. BLAS
does not have a fused three-vector `z = x + y` operation, so `vector_add` uses
the standard two-call sequence:

```text
dcopy(x, z)
daxpy(1.0, y, z)
```

That is not memory-equivalent to Eigen's fused `z = x + y`: it touches about
five vector streams rather than three. Its bandwidth numbers therefore use a
five-stream byte model and should not be compared directly with the fused
Eigen/manual rows as if all rows were doing the same operation.

The sparse part of the benchmark was reduced to a tiny matrix for these runs
because this section only compares dense vector add and dot.

OpenBLAS-linked CBLAS, average of three runs:

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup | j=1 GiB/s | j=16 GiB/s |
|---|---:|---:|---:|---:|---:|---:|
| vector add | `VectorXd` | 0.460145 | 0.446814 | 1.03x | 31.1 | 32.0 |
| vector add | `manual_openmp` | 0.452146 | 0.382184 | 1.18x | 31.6 | 37.4 |
| vector add | `cblas_dcopy_daxpy` | 0.587934 | 0.494521 | 1.19x | 40.7 | 48.4 |
| dot | `VectorXd` | 0.293792 | 0.289972 | 1.01x | 48.7 | 49.3 |
| dot | `manual_openmp` | 0.380414 | 0.175584 | 2.17x | 37.6 | 81.6 |
| dot | `cblas_ddot` | 0.295758 | 0.193759 | 1.53x | 48.4 | 74.0 |

MKL-linked CBLAS, average of three runs:

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup | j=1 GiB/s | j=16 GiB/s |
|---|---:|---:|---:|---:|---:|---:|
| vector add | `VectorXd` | 0.464448 | 0.441362 | 1.05x | 30.8 | 32.4 |
| vector add | `manual_openmp` | 0.467977 | 0.376537 | 1.24x | 30.6 | 38.0 |
| vector add | `cblas_dcopy_daxpy` | 0.585745 | 0.346867 | 1.69x | 40.7 | 69.0 |
| dot | `VectorXd` | 0.307613 | 0.293154 | 1.05x | 46.6 | 48.8 |
| dot | `manual_openmp` | 0.381828 | 0.180856 | 2.11x | 37.5 | 79.1 |
| dot | `cblas_ddot` | 0.309784 | 0.180577 | 1.72x | 46.2 | 79.2 |

Interpretation:

- Eigen-native vector add and dot remain effectively serial.
- OpenBLAS `ddot` does thread, but less effectively than the simple manual
  OpenMP reduction on this benchmark.
- MKL `ddot` and the manual OpenMP reduction are essentially tied at 16
  threads.
- MKL `dcopy+daxpy` scales better than OpenBLAS here, but it is still not a
  drop-in replacement for Eigen's fused `c = a + b`; it performs more memory
  traffic and needs two BLAS calls.
- For OAS, BLAS can help selected large dot/reduction kernels. For vector
  elementwise expressions, loop fusion or a small OAS-owned OpenMP vector helper
  is more appropriate than replacing fused Eigen expressions with BLAS
  `dcopy+daxpy`.

## Local OpenBLAS Builds

Built OpenBLAS 0.3.32 locally outside the repo from tag `v0.3.32`
(`8cecf89`) under `/home/beremi/build/OpenBLAS-0.3.32`.

The useful build was BLAS/CBLAS-only, `TARGET=ZEN`, OpenMP threaded:

```bash
make -j16 DYNAMIC_ARCH=0 TARGET=ZEN USE_OPENMP=1 NO_AFFINITY=1 \
  ONLY_CBLAS=1 NUM_THREADS=32 BINARY=64 libs shared
make PREFIX=/home/beremi/local/openblas-0.3.32-zen-openmp \
  ONLY_CBLAS=1 install
```

The benchmark confirmed it was linked against the local library:

```text
libopenblas.so.0 => /home/beremi/local/openblas-0.3.32-zen-openmp/lib/libopenblas.so.0
openblas_config="OpenBLAS 0.3.32 NO_LAPACK NO_LAPACKE NO_AFFINITY USE_OPENMP ZEN MAX_THREADS=32"
openblas_runtime_threads=16 openblas_parallel_mode=2
```

I also tested a local `TARGET=ZEN` pthread build and a local
`TARGET=COOPERLAKE USE_OPENMP=1` build. The pthread build was much slower at
16 threads for these level-1 kernels. The `COOPERLAKE` build did not beat the
local `ZEN` OpenMP build, despite the system dynamic OpenBLAS reporting a
Cooperlake dispatch on this CPU.

Final BLAS rows only, average of three runs:

| Backend | Kernel | j=1 avg s | j=16 avg s | Speedup | j=1 GiB/s | j=16 GiB/s |
|---|---|---:|---:|---:|---:|---:|
| system OpenBLAS configured | `cblas_dcopy_daxpy` | 0.560776 | 0.493502 | 1.14x | 42.6 | 48.3 |
| system OpenBLAS configured | `cblas_ddot` | 0.293145 | 0.192185 | 1.53x | 48.8 | 74.4 |
| local OpenBLAS `ZEN` OpenMP | `cblas_dcopy_daxpy` | 0.558127 | 0.479279 | 1.16x | 42.7 | 49.8 |
| local OpenBLAS `ZEN` OpenMP | `cblas_ddot` | 0.289480 | 0.182222 | 1.59x | 49.4 | 78.6 |
| local OpenBLAS `COOPERLAKE` OpenMP | `cblas_dcopy_daxpy` | 0.561329 | 0.491216 | 1.14x | 42.5 | 48.5 |
| local OpenBLAS `COOPERLAKE` OpenMP | `cblas_ddot` | 0.289244 | 0.199666 | 1.45x | 49.5 | 71.7 |
| MKL | `cblas_dcopy_daxpy` | 0.585745 | 0.346867 | 1.69x | 40.7 | 69.0 |
| MKL | `cblas_ddot` | 0.309784 | 0.180577 | 1.72x | 46.2 | 79.2 |

Conclusion from the local OpenBLAS pass:

- The local `ZEN` OpenMP build is slightly faster than the system OpenBLAS for
  `ddot` and slightly faster for `dcopy+daxpy`, but the difference is small.
- MKL remains better for `dcopy+daxpy` at 16 threads and essentially tied with
  the best local OpenBLAS/manual OpenMP result for `ddot`.
- Replacing fused Eigen vector expressions with BLAS remains unattractive:
  BLAS `dcopy+daxpy` is not the same operation as `c = a + b` and moves more
  data.

## librsb SpMV Comparison

`librsb` is a shared-memory OpenMP sparse matrix library using the Recursive
Sparse Blocks format. Its upstream documentation describes support for SpMV,
SpMM, triangular solves, row/column scaling, norms, and CSR/CSC/COO
interoperability; it is especially aimed at multicore sparse kernels and
symmetric/transposed multiplication variants.

No saved OAS stiffness/active matrices were found in the repo results tree, so
this section uses the existing synthetic sparse benchmark matrices:

- narrow band: `1,000,000 x 1,000,000`, radius 3, `6,999,988` nonzeros.
- wide band: `200,000 x 200,000`, radius 31, `12,599,008` nonzeros.

Built locally from `librsb-1.3.0.2.tar.gz` under
`/home/beremi/build/librsb-1.3.0.2` and installed to
`/home/beremi/local/librsb-1.3.0.2-native-openmp`:

```bash
./configure \
  --prefix=/home/beremi/local/librsb-1.3.0.2-native-openmp \
  --enable-shared --disable-static --enable-pkg-config-install \
  --disable-programs --disable-c-examples --disable-c++-examples \
  --disable-fortran-examples --disable-sparse-blas-interface \
  --without-zlib --without-xdr --with-max-threads=32 \
  CC=gcc CXX=g++ FC=' ' \
  CFLAGS='-O3 -march=native -DNDEBUG' \
  CXXFLAGS='-include omp.h -O3 -march=native -DNDEBUG -std=c++17' \
  OPENMP_CFLAGS='-fopenmp' OPENMP_CXXFLAGS='-fopenmp'
make -j16
make install
```

The `-include omp.h` workaround is needed with GCC 15 here because this librsb
release otherwise includes `omp.h` inside an `extern "C"` block in one C++
translation unit, while GCC 15's OpenMP header contains C++ templates.

The benchmark binary confirmed the local library and 16-thread runtime:

```text
librsb.so.0 => /home/beremi/local/librsb-1.3.0.2-native-openmp/lib/librsb.so.0
librsb_runtime_threads=16 rsb_num_threads_env=16
```

Narrow band, average of three runs except tuned rows, which are one run:

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup | j=16 GiB/s |
|---|---:|---:|---:|---:|---:|
| SpMV | `ColMajor` | 0.098684 | 0.102933 | 0.96x | 27.2 |
| SpMV | `RowMajor` | 0.090527 | 0.061416 | 1.47x | 45.5 |
| SpMV | `manual_rowmajor_openmp` | 0.092984 | 0.032904 | 2.83x | 85.1 |
| SpMV | `librsb_rsb_spmv` | 0.093685 | 0.051270 | 1.83x | 54.5 |
| setup | `librsb_from_csr` | 0.048953 | 0.050753 | 0.96x | 0.0 |
| SpMV | `librsb_tuned_rsb_spmv` | 0.096362 | 0.048778 | 1.98x | 57.3 |
| setup | `librsb_tuned_from_csr` | 0.112918 | 0.130137 | 0.87x | 0.0 |

Wide band, average of three runs except tuned rows, which are one run:

| Kernel | Variant | j=1 avg s | j=16 avg s | Speedup | j=16 GiB/s |
|---|---:|---:|---:|---:|---:|
| SpMV | `ColMajor` | 0.102768 | 0.104631 | 0.98x | 41.2 |
| SpMV | `RowMajor` | 0.104492 | 0.085360 | 1.22x | 50.6 |
| SpMV | `manual_rowmajor_openmp` | 0.102914 | 0.067129 | 1.53x | 64.4 |
| SpMV | `librsb_rsb_spmv` | 0.097803 | 0.070874 | 1.38x | 60.9 |
| setup | `librsb_from_csr` | 0.073713 | 0.067666 | 1.09x | 0.0 |
| SpMV | `librsb_tuned_rsb_spmv` | 0.093578 | 0.077143 | 1.21x | 55.9 |
| setup | `librsb_tuned_from_csr` | 0.490004 | 0.280691 | 1.75x | 0.0 |

Interpretation:

- `librsb` works and uses OpenMP in this local build.
- It beats Eigen's native RowMajor sparse matvec on both synthetic matrices.
- It does not beat the simple static CRS OpenMP loop on these two matrices.
  These synthetic matrices are very regular, so static row scheduling is close
  to ideal.
- `rsb_tune_spmm()` slightly improves the narrow-band case, but the extra setup
  cost is larger than the gain for short repeated solves. On the wide-band case
  it made the measured SpMV slower.
- For OAS DFGMRES SpMV specifically, `librsb` is not obviously better than the
  current OAS-owned CRS loop. It would add an external dependency and a matrix
  conversion/setup step. It may still be worth retesting on real OAS matrices,
  especially if we need transposed, symmetric-storage, or less uniformly
  balanced sparse products.

## Actual TS-N_65 Replay Matrix

This pass uses the larger dumped linear replay from:

```text
results/amgcl-tsn65-first-solve-20260503-200641/001-hybrid-nullspace-fixed/linear_replay/solve_1/
```

Inputs:

- Matrix: `matrix.mtx`, `567923 x 567923`, `49429959` nonzeros.
- Vectors: actual replay `rhs.tsv` and `solution.tsv`.
- Vector add/axpby/dot: `rhs` with `solution`.
- SpMV: `matrix * solution`.
- Benchmark source: `scripts/benchmark-replay-linear-kernels.cpp`.
- Load/build time is excluded from the kernel rows. It was `8.83 s` at `j=1`
  and `8.69 s` at `j=16`.

All `j=1` and `j=16` checksums matched for the same kernel/variant, up to the
expected last-bit differences in threaded reductions and sparse accumulation.

| Kernel | Variant | j=1 s | j=16 s | Speedup | j=16 GiB/s | Repeats |
|---|---:|---:|---:|---:|---:|---:|
| vector add | `Eigen_VectorXd_rhs_plus_solution` | 0.165009 | 0.158963 | 1.04x | 159.7 | 2000 |
| vector add | `manual_openmp_rhs_plus_solution` | 0.165736 | 0.018131 | 9.14x | 1400.3 | 2000 |
| vector add | `openblas_dcopy_daxpy_rhs_plus_solution` | 0.249563 | 0.388854 | 0.64x | 108.8 | 2000 |
| vector axpby | `Eigen_VectorXd_rhs_solution` | 0.166553 | 0.165391 | 1.01x | 153.5 | 2000 |
| vector axpby | `manual_openmp_rhs_solution` | 0.166879 | 0.018308 | 9.11x | 1386.7 | 2000 |
| dot | `Eigen_VectorXd_rhs_dot_solution` | 0.164347 | 0.165374 | 0.99x | 153.5 | 3000 |
| dot | `manual_openmp_rhs_dot_solution` | 0.630441 | 0.047195 | 13.36x | 537.9 | 3000 |
| dot | `openblas_ddot_rhs_dot_solution` | 0.191035 | 0.020224 | 9.45x | 1255.3 | 3000 |
| SpMV | `Eigen_RowMajor_matrix_times_solution` | 1.777139 | 1.432769 | 1.24x | 39.1 | 100 |
| SpMV | `manual_rowmajor_openmp_matrix_times_solution` | 2.075918 | 1.426758 | 1.45x | 39.3 | 100 |
| SpMV | `Eigen_ColMajor_matrix_times_solution` | 1.899200 | 2.063405 | 0.92x | 27.2 | 100 |
| setup | `librsb_from_rowmajor_csr` | 0.464640 | 0.463479 | 1.00x | 0.0 | 1 |
| SpMV | `librsb_rsb_spmv_matrix_times_solution` | 2.432017 | 1.357681 | 1.79x | 41.3 | 100 |

Interpretation:

- Eigen dense vector expressions and Eigen dot remain effectively serial on the
  actual OAS vector size.
- The OAS-size vectors are only `567923` doubles, so repeated vector operations
  are warm-cache microbenchmarks. Manual OpenMP and OpenBLAS `ddot` scale well
  here, but this does not imply the same speedup for every solver vector phase.
- OpenBLAS `ddot` is the best dot kernel in this replay test.
- OpenBLAS `dcopy+daxpy` is not useful for `c = a + b` here; it is slower than
  serial Eigen because it performs two BLAS calls and more memory traffic than a
  fused vector expression.
- On the real matrix, SpMV scaling is weak. RowMajor Eigen and the direct
  OpenMP CRS loop both land around `39 GiB/s` at 16 threads. `librsb` is the
  fastest measured SpMV variant, but only by about `5%` versus the simple CRS
  loop at `j=16`, and it has a one-time CSR-to-RSB setup cost.
- ColMajor Eigen SpMV is serial and gets slightly slower at `j=16`; this
  matches the earlier Eigen source-level finding.
- The actual OAS matrix is much less friendly than the synthetic banded test.
  The bottleneck is sparse memory locality and irregular row structure, not
  merely whether OpenMP threads are launched.

### librsb Option Sweep On TS-N_65

I checked whether the first replay test left obvious librsb performance on the
table.

Build/runtime settings already used correctly:

- Local librsb `1.3.0.2`, built with `-O3 -march=native -DNDEBUG`.
- OpenMP recursive kernels enabled.
- Maximum supported librsb threads: `32`.
- Detected memory hierarchy: `L3:16/64/32768K,L2:16/64/1024K,L1:12/64/48K`.
- Runtime executing thread count explicitly set through
  `RSB_IO_WANT_EXECUTING_THREADS`.

Raw sweep files:

- `j16-librsb-flags-sweep-100.csv`
- `j1/j2/j4/j8/j16-librsb-thread-sweep.csv`
- `j16-librsb-confirm-plain.csv`
- `j8-librsb-confirm-plain.csv`
- `j16-librsb-confirm-tune.csv`
- `j16-librsb-confirm-tune-threads.csv`

Storage flag sweep, `j=16`, 100 SpMV repeats:

| Variant | setup s | SpMV s / 100 | GiB/s |
|---|---:|---:|---:|
| `RSB_FLAG_NOFLAGS` | 0.424 | 1.287102 | 43.6 |
| `RSB_FLAG_SORTED_INPUT` | 0.427 | 1.337069 | 41.9 |
| `RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS` | 0.422 | 1.300932 | 43.1 |
| `DEFAULT_RSB | SORTED_INPUT` | 0.414 | 1.338842 | 41.9 |
| `DEFAULT_RSB | SORTED_INPUT | MORE_LEAVES` | 0.419 | 1.296247 | 43.3 |
| `DEFAULT_RSB | SORTED_INPUT | MORE_DIAG` | 0.422 | 1.329519 | 42.2 |

Thread count sweep, untuned, 50 SpMV repeats:

| Threads | Eigen RowMajor s | manual CRS s | librsb s | fastest librsb GiB/s |
|---:|---:|---:|---:|---:|
| 1 | 0.735382 | 0.801596 | 1.121085 | 25.0 |
| 2 | 0.658006 | 0.690798 | 0.714195 | 39.3 |
| 4 | 0.667336 | 0.674041 | 0.679956 | 41.2 |
| 8 | 0.688695 | 0.674491 | 0.615893 | 45.5 |
| 16 | 0.656598 | 0.678399 | 0.653511 | 42.9 |

Longer confirmation, 100 SpMV repeats:

| Case | setup/tune s | SpMV s / 100 | GiB/s |
|---|---:|---:|---:|
| plain librsb, 16 threads | 0.409 | 1.290319 | 43.5 |
| plain librsb, 8 threads | 0.393 | 1.257839 | 44.6 |
| structure tuning, launched from 16 threads | 1.571 | 1.323950 | 42.4 |
| structure+thread tuning, launched from 16 threads | 4.150 | 1.315696 | 42.6 |

Interpretation:

- `RSB_FLAG_NOFLAGS` is the best setting for this CSR import. It already
  creates the expected internal RSB structure with COO+CSR storage and halfword
  local indices; explicit default flags did not improve it.
- `RSB_FLAG_SORTED_INPUT` did not help, even though the Eigen RowMajor CSR is
  sorted. It is slightly slower in the long sweep.
- `rsb_tune_spmm()` is not worthwhile here. It can find different structures,
  but the longer confirmation was slower and the setup cost is too high for
  repeated Newton matrices.
- librsb is memory/locality limited on this OAS matrix. The best measured
  setting is plain librsb at 8 threads, but the gain over plain 16-thread
  librsb is only about `2.6%`.
- If librsb were integrated into OAS, the sensible first setting would be
  `RSB_FLAG_NOFLAGS` plus an internal librsb thread cap around 8 for this
  matrix family, not autotuning on every matrix update.

## Effective Bandwidth

This is an approximate bytes-touched model, useful only for comparing the two
thread counts inside this benchmark. Values below are from the repaired run.

| Kernel | Variant | j=1 GiB/s | j=16 GiB/s | Ratio |
|---|---:|---:|---:|---:|
| vector add | `VectorXd` | 31.9 | 33.3 | 1.04x |
| vector add | `manual_openmp` | 32.6 | 38.7 | 1.19x |
| vector subtract | `VectorXd` | 32.2 | 31.8 | 0.99x |
| vector subtract | `manual_openmp` | 32.7 | 39.2 | 1.20x |
| vector axpby | `VectorXd` | 33.3 | 32.3 | 0.97x |
| vector axpby | `manual_openmp` | 32.8 | 39.2 | 1.19x |
| dot | `VectorXd` | 49.1 | 47.6 | 0.97x |
| dot | `manual_openmp` | 37.8 | 78.4 | 2.07x |
| sparse matvec | `RowMajor` | 30.7 | 48.2 | 1.57x |
| sparse matvec | `manual_rowmajor_openmp` | 30.7 | 99.3 | 3.24x |
| sparse matvec | `ColMajor` | 28.1 | 28.2 | 1.00x |

## Interpretation

- OpenMP is working in the benchmark. The manual OpenMP kernels do get faster at
  16 threads.
- Eigen dense vector add/sub/axpby does not use native OpenMP. The local Eigen
  multithreading documentation lists selected kernels such as dense
  matrix-matrix and RowMajor sparse-dense products; dense vector elementwise
  expressions are not on that list.
- Eigen `VectorXd::dot()` is implemented through Eigen's redux path and does not
  use OpenMP here.
- Eigen sparse matvec with `ColMajor` storage does not use OpenMP.
- Eigen sparse matvec with `RowMajor` storage does enter an OpenMP path, but it
  uses Eigen's generic scheduling. On this uniform CRS-like matrix, a direct
  static-scheduled RowMajor CRS loop is about 2x faster than Eigen's RowMajor
  OpenMP path at 16 threads.
- Dense vector operations remain memory-bandwidth bound. Even the manual OpenMP
  vector loops only improve by about 1.2x on this machine.
- For OAS DFGMRES, relying on Eigen-native OpenMP is not enough for dense vector
  reductions and is only moderately useful for RowMajor SpMV. A targeted
  OAS-owned CRS matvec with static scheduling is likely a better SpMV direction,
  provided it is guarded by strict consistency tests against `A * x`.

## Source-Level Cause

- `doc/TopicMultithreading.dox` in the vendored Eigen tree lists
  `row-major-sparse * dense vector/matrix products` as a parallelized case, but
  not dense vector elementwise operations or dense dot products.
- `Eigen/src/Core/Dot.h` delegates to `dot_impl`, and `Eigen/src/Core/Redux.h`
  implements reductions without OpenMP pragmas.
- `Eigen/src/SparseCore/SparseDenseProduct.h` contains OpenMP pragmas for the
  RowMajor sparse-dense path. The ColMajor sparse-dense path is serial.

## OAS Repair Applied

- DFGMRES now uses an OAS-owned row-major CRS matvec helper for active-matrix
  products instead of `impl->activeMatrix * v`.
- The helper uses `schedule(static)` and falls back to the serial loop for small
  matrices or one-thread runs.
- The existing CRS arrays are still built from `impl->activeMatrix`, so the
  preconditioner setup and DFGMRES matvec use the same active matrix state.
- When `amgcl_check_matrix 1` is set, DFGMRES compares the CRS matvec against
  Eigen sparse matvec during factorization and fails setup if the relative
  difference is larger than `1e-10`.

## Reproduction

Build:

```bash
c++ -std=c++17 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG -fopenmp \
  -I../oas_deflation-build/release/_deps/eigen-src \
  scripts/benchmark-eigen-openmp-kernels.cpp \
  -o /tmp/oas_eigen_openmp_bench
```

Build with OpenBLAS CBLAS rows:

```bash
c++ -std=c++17 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG -fopenmp \
  -DOAS_EIGEN_BENCH_USE_CBLAS \
  -I../oas_deflation-build/release/_deps/eigen-src \
  scripts/benchmark-eigen-openmp-kernels.cpp \
  $(pkg-config --cflags --libs cblas) \
  -o /tmp/oas_eigen_openmp_bench_blas
```

Build with MKL CBLAS rows:

```bash
c++ -std=c++17 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG -fopenmp \
  -DOAS_EIGEN_BENCH_USE_CBLAS -DOAS_EIGEN_BENCH_USE_MKL_CBLAS \
  -I../oas_deflation-build/release/_deps/eigen-src \
  scripts/benchmark-eigen-openmp-kernels.cpp \
  $(pkg-config --cflags --libs mkl-dynamic-lp64-gomp) \
  -o /tmp/oas_eigen_openmp_bench_mkl_blas
```

Build with local `librsb` SpMV rows:

```bash
c++ -std=c++17 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG -fopenmp \
  -DOAS_EIGEN_BENCH_USE_LIBRSB \
  -I../oas_deflation-build/release/_deps/eigen-src \
  scripts/benchmark-eigen-openmp-kernels.cpp \
  $(PKG_CONFIG_PATH=/home/beremi/local/librsb-1.3.0.2-native-openmp/lib/pkgconfig \
    pkg-config --cflags --libs librsb) \
  -Wl,-rpath,/home/beremi/local/librsb-1.3.0.2-native-openmp/lib \
  -o /tmp/oas_eigen_openmp_bench_librsb
```

Run:

```bash
mkdir -p results/eigen-openmp-kernels-20260507-v3
OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  /tmp/oas_eigen_openmp_bench \
  > results/eigen-openmp-kernels-20260507-v3/j1.csv \
  2> results/eigen-openmp-kernels-20260507-v3/j1.log

OMP_NUM_THREADS=16 OMP_PROC_BIND=close OMP_PLACES=cores \
  /tmp/oas_eigen_openmp_bench \
  > results/eigen-openmp-kernels-20260507-v3/j16.csv \
  2> results/eigen-openmp-kernels-20260507-v3/j16.log
```

For BLAS-only vector/dot runs, also set the BLAS thread controls and reduce the
sparse portion:

```bash
OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 \
  OMP_PROC_BIND=close OMP_PLACES=cores \
  OAS_EIGEN_BENCH_SPARSE_SIZE=16 OAS_EIGEN_BENCH_SPMV_REPEATS=1 \
  /tmp/oas_eigen_openmp_bench_mkl_blas \
  > results/eigen-openmp-kernels-20260507-mkl-blas/j16.csv \
  2> results/eigen-openmp-kernels-20260507-mkl-blas/j16.log
```

For `librsb` SpMV runs:

```bash
OMP_NUM_THREADS=16 RSB_NUM_THREADS=16 \
  OMP_PROC_BIND=close OMP_PLACES=cores \
  /tmp/oas_eigen_openmp_bench_librsb \
  > results/eigen-openmp-kernels-20260507-librsb/j16.csv \
  2> results/eigen-openmp-kernels-20260507-librsb/j16.log
```

Set `OAS_EIGEN_BENCH_LIBRSB_TUNE=1` to include `rsb_tune_spmm()` in the
`librsb` setup phase.
