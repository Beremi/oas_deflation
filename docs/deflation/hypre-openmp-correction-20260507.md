# HYPRE OpenMP Threading Correction - 2026-05-07

## Finding

The DFGMRES HYPRE preconditioner path was previously wrapped in an OpenMP
thread clamp of one thread for both `HYPRE_BoomerAMGSetup()` and
`HYPRE_BoomerAMGSolve()`. Removing that clamp blindly made HYPRE use 16 threads,
but it changed the effective BoomerAMG preconditioner for the TS-N_65 settings:

| run | HYPRE smoother threads | matrix insertion threads | DFGMRES outer iterations | step wall s | solve s |
| --- | ---: | ---: | ---: | ---: | ---: |
| previous good j16 | 1 | row-by-row serial | 37 | 120.014 | 15.3175 |
| unrestricted HYPRE j16 | 16 | 16 | 92 | 189.291 | 39.0196 |
| corrected policy j16 | 1 | 16 | 37 | 116.895 | 14.8880 |

The local HYPRE source explains why this happens. Relaxation type `6` is
hybrid symmetric Gauss-Seidel/SSOR. In `par_relax.c`, when
`hypre_NumThreads() > 1`, HYPRE partitions rows across OpenMP threads and runs a
different threaded smoother path. That is not numerically equivalent to the
serial smoother and was weaker for this benchmark.

## Code Change

- Added `hypre_threads` / `hypre_thread_count` solver input option.
  - `0` means automatic.
  - Positive values explicitly cap HYPRE's BoomerAMG setup/apply threads.
- Automatic HYPRE BoomerAMG setup/apply thread policy:
  - Use OpenMP threads for thread-stable smoother types: weighted Jacobi `0`,
    Jacobi `7`, Chebyshev `16`, and l1-Jacobi `18`.
  - Use one thread for hybrid Gauss-Seidel/SSOR-style smoothers by default,
    including the TS-N_65 `hypre_relax_type 6`.
- Kept HYPRE IJ matrix insertion parallel-capable:
  - Bulk row insertion uses `HYPRE_IJMatrixSetOMPFlag(..., 1)`.
  - The corrected TS-N_65 run logs `hypre_matrix_threads=16`.
- Added setup logs for `hypre_threads`, `hypre_matrix_threads`, and
  `hypre_relax_type`.

## Validation

CTest:

| command | result |
| --- | --- |
| `OMP_NUM_THREADS=1 ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |
| `OMP_NUM_THREADS=16 ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |

TS-N_65 first loading step, DFGMRES + HYPRE, 16 threads:

| metric | previous good j16 | corrected policy j16 |
| --- | ---: | ---: |
| nonlinear iterations | 19 | 19 |
| linear solves | 19 | 19 |
| DFGMRES outer iterations | 37 | 37 |
| max linear iterations | 4 | 4 |
| final true relative residual | 0.0622639 | 0.0622639 |
| step wall seconds | 120.014 | 116.895 |
| linear solve seconds | 15.3175 | 14.8880 |
| preconditioner setup seconds | 10.0443 | 9.59005 |
| preconditioner apply seconds | 12.8283 | 12.4608 |

The production behavior is now explicit: HYPRE OpenMP is used automatically only
when the selected smoother is suitable for it. For TS-N_65's current
`hypre_relax_type 6`, the matrix build can use OpenMP, but the BoomerAMG
preconditioner setup/apply remains serial to preserve convergence.

## Follow-Up: Threaded Smoother Sweep

The next question was whether HYPRE could scale without breaking the nonlinear
solve. The answer is yes, but not with the original smoother. The original
TS-N_65 `hypre_relax_type 6` path is stable only when BoomerAMG setup/apply is
kept effectively serial. Thread-friendly smoothers do scale per apply, but some
are much weaker as preconditioners and increase DFGMRES work too much.

I also fixed the DFGMRES HYPRE apply path to allocate and reuse IJ vectors during
factorization. That avoids per-apply vector construction and keeps the HYPRE
preconditioner path explicit. It did not materially change TS-N_65 timing by
itself, so the bottleneck is the BoomerAMG smoother/apply work, not wrapper
allocation.

TS-N_65 first loading step, 16 OpenMP threads:

| configuration | HYPRE solve threads | relax type | Cheby order | nonlinear iters | DFGMRES applies | max linear iter | precond setup s | precond apply s | mean apply s | solve s | profile wall s | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| broken forced OpenMP | 16 | 6 | - | 19 | 92 | 16 | 8.391 | 33.191 | 0.3608 | 39.020 | 189.291 | reject: iterations more than double |
| guarded original | 1 | 6 | - | 19 | 37 | 4 | 9.590 | 12.461 | 0.3368 | 14.888 | 116.895 | stable, but HYPRE apply is serial |
| reused vectors, guarded original | 1 | 6 | - | 19 | 37 | 4 | 9.580 | 12.452 | 0.3365 | 14.959 | 120.008 | allocation was not the bottleneck |
| l1-Jacobi | 16 | 18 | - | 28 | 367 | 45 | 4.267 | 23.623 | 0.0644 | 34.312 | 140.275 | reject: much weaker preconditioner |
| Chebyshev default | 16 | 16 | hypre default | 19 | 73 | 9 | 4.599 | 8.460 | 0.1159 | 11.618 | 110.374 | faster, but applications nearly double |
| Chebyshev order 3 | 16 | 16 | 3 | 19 | 58 | 9 | 4.649 | 9.266 | 0.1598 | 12.216 | 111.384 | usable, but worse than order 4 |
| Chebyshev order 4 | 16 | 16 | 4 | 19 | 52 | 8 | 4.636 | 10.470 | 0.2013 | 13.194 | 110.676 | first threaded candidate |
| Chebyshev order 4, 2 AMG iters | 16 | 16 | 4 | 18 | 51 | 8 | 4.762 | 12.118 | 0.2376 | 14.878 | 114.383 | no benefit over one AMG iter |
| relax 8 forced OpenMP | 16 | 8 | - | 20 | 140 | 15 | 4.211 | 16.987 | 0.1213 | 21.490 | 119.616 | reject |
| relax 6 forced to 4 threads | 4 | 6 | - | 19 | 55 | 8 | 5.483 | 12.089 | 0.2198 | 14.802 | 115.929 | little gain, more iterations |
| relax 6 forced to 2 threads | 2 | 6 | - | 19 | 54 | 8 | 7.153 | 15.730 | 0.2913 | 18.372 | 118.148 | reject |

The first threaded candidate was `hypre_relax_type 16` with
`hypre_cheby_order 4` and one BoomerAMG iteration. It kept the nonlinear solve
stable and kept DFGMRES applications below a 2x increase relative to the
guarded original path: 52 applies instead of 37.

Follow-up tuning added `hypre_num_sweeps` and found a stronger threaded
candidate:

| configuration | HYPRE solve threads | relax type | Cheby order | HYPRE sweeps | nonlinear iters | DFGMRES applies | max linear iter | precond setup s | precond apply s | mean apply s | solve s | linear total s | profile wall s | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| guarded original | 1 | 6 | - | default | 19 | 37 | 4 | 9.590 | 12.461 | 0.3368 | 14.888 | 27.227 | 120.415 | original iteration count |
| Cheby order 3, 3 sweeps | 16 | 16 | 3 | 3 | 19 | 37 | 5 | 4.330 | 14.830 | 0.4008 | 17.179 | 25.059 | 118.159 | selected 37-outer threaded candidate |
| Cheby order 4, 2 sweeps | 16 | 16 | 4 | 2 | 19 | 40 | 6 | 4.406 | 14.248 | 0.3562 | 16.722 | 24.586 | 118.007 | faster, but not 37 outer |
| Cheby order 4, 3 sweeps | 16 | 16 | 4 | 3 | 19 | 33 | 5 | 4.367 | 16.867 | 0.5111 | 19.068 | 26.901 | 119.870 | stronger than requested |

The selected 37-outer threaded candidate is therefore:

```text
hypre_relax_type 16
hypre_cheby_order 3
hypre_num_sweeps 3
hypre_boomer_max_iterations 1
hypre_threads 0
```

It preserves the original DFGMRES application count (`37`) while allowing HYPRE
BoomerAMG setup/apply to run with 16 threads. The apply itself is more expensive
than the guarded original because it performs three Chebyshev sweeps, but setup,
matvec, and total linear time improve enough to reduce first-step wall time in
this run.

Apples-to-apples selected-candidate comparison:

| metric | j=1 | j=16 | speedup / change |
| --- | ---: | ---: | ---: |
| run wall seconds | 208.662 | 118.159 | 1.77x |
| linear total seconds | 36.672 | 25.059 | 1.46x |
| solve seconds | 24.576 | 17.179 | 1.43x |
| setup/analyze seconds | 12.096 | 7.880 | 1.53x |
| nonlinear iterations | 19 | 19 | unchanged |
| linear solves | 19 | 19 | unchanged |
| DFGMRES applies | 37 | 37 | unchanged |
| max linear iterations | 5 | 5 | unchanged |
| final true relative residual | 0.0688078 | 0.0688078 | unchanged |
| preconditioner setup seconds | 8.484 | 4.330 | 1.96x |
| preconditioner apply seconds | 21.756 | 14.830 | 1.47x |
| mean preconditioner apply seconds | 0.5880 | 0.4008 | 1.47x |
| matvec seconds | 1.456 | 0.944 | 1.54x |

Apples-to-apples Chebyshev order 4 comparison:

| metric | j=1 | j=16 | speedup / change |
| --- | ---: | ---: | ---: |
| profile wall seconds | 199.301 | 110.676 | 1.80x |
| run wall seconds | 202.943 | 114.349 | 1.78x |
| nonlinear iterations | 19 | 19 | unchanged |
| linear solves | 19 | 19 | unchanged |
| DFGMRES applies | 52 | 52 | unchanged |
| max linear iterations | 8 | 8 | unchanged |
| final true relative residual | 0.0886686 | 0.0886686 | unchanged |
| preconditioner setup seconds | 10.417 | 4.636 | 2.25x |
| preconditioner apply seconds | 14.939 | 10.470 | 1.43x |
| mean preconditioner apply seconds | 0.2873 | 0.2013 | 1.43x |
| total linear solve seconds | 18.465 | 13.194 | 1.40x |

## Current Code State

- Generic HYPRE defaults remain conservative: `hypre_relax_type 6`,
  automatic `hypre_threads 0`, and no Chebyshev override.
- Automatic HYPRE solve threading now only enables OpenMP for thread-stable
  smoothers (`0`, `7`, `16`, `18`). Non-thread-stable smoothers use one HYPRE
  solve thread unless the user explicitly overrides `hypre_threads`.
- `hypre_cheby_order`, `hypre_cheby_fraction`, and `hypre_num_sweeps` are now solver-file options.
- The TS-N_65 DFGMRES/HYPRE campaign script now defaults to the validated
  37-outer threaded candidate: `hypre_relax_type 16`,
  `hypre_cheby_order 3`, `hypre_num_sweeps 3`,
  `hypre_boomer_max_iterations 1`, and `hypre_threads 0`.
- The setup log now includes HYPRE solve threads, matrix insertion threads,
  relax type, Chebyshev order, Chebyshev fraction, and BoomerAMG max iterations.

Validation after the final code changes:

| command | result |
| --- | --- |
| `OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |
| `OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |
