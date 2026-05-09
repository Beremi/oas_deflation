# DFGMRES Active Ordering Experiment, 2026-05-07

## What changed

The permutation is now applied to the whole native `DeflatedFGMRES` active system, not only to the HYPRE preconditioner matrix.

- Added solver input key `dfgmres_elastic_reorder`.
- Kept `hypre_elastic_reorder` as a compatibility fallback for existing TS-N_65 scripts.
- `DeflatedFGMRES::factorize()` now permutes the lifted active matrix before diagonal scaling, CRS extraction, librsb construction, and HYPRE setup.
- RHS vectors, collected deflation vectors, and final solutions are mapped at the DFGMRES boundary:
  - reduced OAS order -> lifted elastic order -> active permuted order,
  - active permuted order -> lifted elastic order -> reduced OAS order.
- HYPRE no longer has a private per-apply permutation path; its matrix and vectors are already in active DFGMRES order.
- AMGCL near-nullspace rows are permuted if the active matrix is permuted.

The implemented ordering mode used here is coordinate-sorted node-major ordering:

```text
old node order -> stable sort by x, y, z -> [ux, uy, uz] per node
```

This keeps the cyclic node-major vector-component convention expected by HYPRE system AMG. HYPRE's BoomerAMG documentation states that for systems of PDEs, `num_functions > 1` generates a default repeating `dof_func` sequence and supports nodal coarsening for unknowns at the same physical point. MFEM's HYPRE wrapper follows the same distinction by using system options and explicit `DofFunc` only when the ordering differs from the default node ordering.

Sources:

- HYPRE BoomerAMG system AMG notes: https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html
- MFEM `HypreBoomerAMG::SetSystemsOptions`: https://docs.mfem.org/html/hypre_8cpp_source.html

## Correctness

Build:

```bash
PKG_CONFIG_PATH=/home/beremi/local/librsb-1.3.0.2-native-openmp/lib/pkgconfig \
cmake --build ../oas_deflation-build/release --target OAS -j16
```

CTest:

| test command | result |
| --- | --- |
| `OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |
| `OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure` | 7/7 passed |

Note: `src/tests/check.py` was made more robust because the smoke check was trying to parse wrapped text headers such as `uy`/`fy` as numeric data when no `check_results` directory exists.

## TS-N_65 First-Step Runs

Command shape:

```bash
python3 scripts/run-tsn65-dfgmres-hypre-deflation.py \
  --oas-bin ../oas_deflation-build/release/bin/OAS \
  --total-time 5.000000e-03 \
  --linear-tol 1e-1 --true-tol 1e-1 \
  --nvecs 20 --only dfgmres-hypre-tol1em1-N20 \
  --hypre-relax-type 16 --hypre-cheby-order 4 \
  --hypre-boomer-max-iterations 1 \
  --hypre-elastic-reorder 2
```

The current script writes both `dfgmres_elastic_reorder 2` and `hypre_elastic_reorder 2`.

## Comparison

| run | threads | wall s | step s | linear solve s | setup s | precond apply s | matvec s | DFGMRES outer | max iter | final true relres |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| previous Cheby4 | 1 | 202.943 | 199.301 | 18.465 | 13.299 | 14.939 | 2.113 | 52 | 8 | 0.0886686 |
| HYPRE-only reorder | 1 | 204.703 | 201.193 | 18.306 | 12.890 | 14.732 | 2.108 | 52 | 8 | 0.0866123 |
| global DFGMRES reorder | 1 | 205.693 | 201.947 | 18.407 | 12.434 | 15.031 | 1.829 | 52 | 8 | 0.0866123 |
| previous Cheby4 | 16 | 114.349 | 110.676 | 13.194 | 7.368 | 10.470 | 1.285 | 52 | 8 | 0.0886686 |
| HYPRE-only reorder | 16 | 114.696 | 111.125 | 13.003 | 8.561 | 10.163 | 1.303 | 52 | 8 | 0.0866123 |
| global DFGMRES reorder | 16 | 114.340 | 110.671 | 12.772 | 7.796 | 10.069 | 1.148 | 52 | 8 | 0.0866123 |

Artifacts:

- Previous j=1: `results/tsn65-dfgmres-hypre-relax16-cheby4-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20`
- Previous j=16: `results/tsn65-dfgmres-hypre-relax16-cheby4-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`
- HYPRE-only reorder j=1: `results/tsn65-dfgmres-hypre-reorder-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20`
- HYPRE-only reorder j=16: `results/tsn65-dfgmres-hypre-reorder-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`
- Global reorder j=1: `results/tsn65-dfgmres-hypre-global-reorder-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20`
- Global reorder j=16: `results/tsn65-dfgmres-hypre-global-reorder-firststep-20260507-j16-nobind/dfgmres-hypre-tol1em1-N20`

## Interpretation

The global permutation is numerically safe in this first-step test:

- nonlinear iterations unchanged: 19,
- linear solves unchanged: 19,
- DFGMRES outer iterations unchanged: 52,
- max DFGMRES iterations unchanged: 8,
- final true residual stays within the same loose target.

For 16 threads, the active global permutation improved the linear solve components modestly:

- solve time: `13.194 s -> 12.772 s`,
- preconditioner apply: `10.470 s -> 10.069 s`,
- librsb matvec: `1.285 s -> 1.148 s`.

The whole first-step wall time is effectively unchanged because forces and model initialization dominate the run. For 1 thread, wall time is slightly worse, while matvec is faster but preconditioner apply is slightly slower.

One run with `OMP_PROC_BIND=close OMP_PLACES=cores` produced a much slower wall time (`170.731 s`) because force integration slowed dramatically. That run is kept as an artifact at `results/tsn65-dfgmres-hypre-global-reorder-firststep-20260507-j16`, but it is not the fair comparison against the previous no-affinity artifacts.

## Decision

Keep the implementation, but keep the core solver default as no active reorder (`dfgmres_elastic_reorder = 0`). The TS-N_65 campaign script currently opts into mode `2` so this experiment is reproducible.

The measured benefit is real but small in the first-step wall time. The most promising part is the lower 16-thread matvec and preconditioner apply time without any iteration increase.
