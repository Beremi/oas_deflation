# TS-N65 Baseline Replication

This note records the reproducible recipe for the strict 8-quarter-step TS-N65
baseline. The important lesson from the local reruns is that the baseline is not
defined by the short nonlinear solver block alone. It requires the full
DFGMRES/HYPRE control block below.

## Known Good Result

The target nonlinear iteration sequence is:

```text
6, 6, 10, 13, 17, 183, 187, 163
```

The total number of printed nonlinear iteration rows is `585`, counting
iteration `0` in each step. A successful run reaches `total_time = 1.0e-2` in
exactly `8/8` fixed quarter steps.

Known successful code states:

- `main` at `779608498250c976d807c7d8c8501e66ea48d4ae`
- `nonlinear_solver_testing` checkpoint
  `16c60d903808d69c7a06273d9f6d445f3603bae9`

Runtime varies by machine and build. The iteration sequence is the primary
replication criterion.

## Build

Build an out-of-tree Release executable with HYPRE enabled and VTK disabled:

```bash
cmake -S . -B /tmp/oas_tsn65_full_baseline_build \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_VTK=OFF \
  -DUSE_HYPRE=ON

cmake --build /tmp/oas_tsn65_full_baseline_build --target OAS -j 16
```

The executable banner should report a Release build with OpenMP, BLAS, and
LAPACK:

```bash
/tmp/oas_tsn65_full_baseline_build/bin/OAS
```

The command exits with a usage message when no input file is supplied, but the
version banner is still printed.

## Input Deck

Use the TS-N65 input deck. In a fully populated checkout this may be available
as `data/cases/TS-N_65`. In the local benchmark-import workflow used for the
2026-05-20 runs, the extracted deck was at:

```text
results/tsn65-zip-runs-20260519-190000/input/TS-N_65/
```

Copy the deck to a clean run directory, for example:

```bash
RUN_ROOT=results/tsn65-full-baseline-replication-$(date +%Y%m%d-%H%M%S)
RUN_DIR="$RUN_ROOT/runs/current-branch-full-baseline"

mkdir -p "$RUN_DIR"
cp -a results/tsn65-zip-runs-20260519-190000/input/TS-N_65/. "$RUN_DIR"/
```

Replace `solver.inp` in the run directory with the full baseline deck below.

## Solver Input

```text
SteadyStateNonLinearSolver
time_step 1.25e-3
total_time 1e-2

solver_type DeflatedFGMRES
max_iterations 300
stiffness_matrix_iter_update 10
stiffness_matrix_step_update 1
limit_tolerance 0
tolerance 1e-3

dfgmres_tolerance 1e-1
dfgmres_true_tolerance 1e-1
dfgmres_max_iterations 500
dfgmres_restart 80
dfgmres_deflation_vectors 20
dfgmres_deflation_eps 1e-15
dfgmres_collect_newton_steps 1
dfgmres_preconditioner hypre
dfgmres_reorthogonalize_on_matrix_change 1
dfgmres_elastic_reorder 2

hypre_coarsen_type 8
hypre_interp_type 6
hypre_strong_threshold 0.25
hypre_nodal 4
hypre_relax_type 16
hypre_num_sweeps 3
hypre_p_max 4
hypre_boomer_max_iterations 1
hypre_cheby_order 3
hypre_cheby_fraction -1
hypre_elastic_reorder 2
hypre_threads 0
hypre_use_dof_functions 0
hypre_use_interp_vectors 0

max_time_step 1.25e-3
min_time_step 1.25e-3
stiff_matrix_type tangent
```

`hypre_threads 0` means use the active OpenMP thread count. With
`OMP_NUM_THREADS=16`, the solver setup line should resolve to
`hypre_threads=16`.

## Run

Run from inside the copied TS-N65 deck directory:

```bash
OMP_NUM_THREADS=16 \
OMP_DYNAMIC=FALSE \
OMP_PROC_BIND=close \
OMP_PLACES=cores \
/tmp/oas_tsn65_full_baseline_build/bin/OAS -j 16 master.inp > solver.out 2>&1
```

The `-j 16` argument is kept for consistency with the baseline workflow, even
though the DFGMRES/HYPRE path is controlled by OpenMP and the solver deck.

## Replication Checkpoints

The first `DeflatedFGMRES setup complete` line must show the HYPRE
preconditioner, coordinate-node-major elastic reordering, and 16 HYPRE threads:

```text
preconditioner=hypre_boomeramg_apply
elastic_reorder=coordinate_node_major
hypre_threads=16
```

If the setup line instead reports `elastic_reorder=none` or `hypre_threads=1`,
the full baseline deck is not active.

The per-step iteration checkpoints are:

| step | time | target rows | final residual | final displacement | final energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `0.00125` | 6 | `7.705844e-06` | `6.743608e-04` | `7.134058e-05` |
| 2 | `2.500000e-03` | 6 | `4.198577e-04` | `1.021696e-04` | `6.255293e-04` |
| 3 | `3.750000e-03` | 10 | `7.308889e-04` | `1.346596e-04` | `9.259996e-04` |
| 4 | `5.000000e-03` | 13 | `6.064366e-04` | `1.451312e-04` | `8.371949e-04` |
| 5 | `6.250000e-03` | 17 | `6.189476e-04` | `1.136484e-04` | `9.322705e-04` |
| 6 | `7.500000e-03` | 183 | `6.428869e-04` | `8.777311e-05` | `9.986343e-04` |
| 7 | `8.750000e-03` | 187 | `6.721379e-04` | `1.022238e-04` | `9.840756e-04` |
| 8 | `1.000000e-02` | 163 | `5.838771e-04` | `1.075313e-04` | `9.817164e-04` |

During a healthy run, steps 6, 7, and 8 are long. Step 6 should not be judged a
failure merely because it passes 100 nonlinear rows. It should close at 183
rows without NaNs or DFGMRES true-tolerance warnings.

## Parsing The Log

Use this command to summarize the printed nonlinear rows and final error row for
each step:

```bash
awk '
/######### Solving step/ {
  if (step != "") {
    print step "\t" count "\t" last
  }
  step=$4
  count=0
  last=""
  next
}
/^[[:space:]]*[0-9]+[[:space:]]/ {
  count++
  last=$0
}
END {
  if (step != "") {
    print step "\t" count "\t" last
  }
}
' solver.out
```

Expected output shape:

```text
1  6    5   7.705844e-06   6.743608e-04   7.134058e-05
2  6    5   4.198577e-04   1.021696e-04   6.255293e-04
3  10   9   7.308889e-04   1.346596e-04   9.259996e-04
4  13   12  6.064366e-04   1.451312e-04   8.371949e-04
5  17   16  6.189476e-04   1.136484e-04   9.322705e-04
6  183  182 6.428869e-04   8.777311e-05   9.986343e-04
7  187  186 6.721379e-04   1.022238e-04   9.840756e-04
8  163  162 5.838771e-04   1.075313e-04   9.817164e-04
```

Also verify the run completed cleanly:

```bash
rg "did not converge|warning: performed|nan|total duration|end of calculation" solver.out
```

A successful baseline replication should contain `end of calculation` and
`total duration`, but no `did not converge`, no `warning: performed`, and no
`nan`.

## Reduced Deck Failure Mode

The earlier local failure used a shortened `solver.inp` that omitted the
DFGMRES/HYPRE controls. That changed the linear solve and produced this pattern:

- setup line used `elastic_reorder=none` and `hypre_threads=1`;
- steps 1-5 converged with a different sequence: `8, 9, 12, 14, 18`;
- step 6 eventually emitted DFGMRES true-tolerance warnings;
- NaNs began at step 6 iteration 134;
- the run stopped at `max_iterations=300` for step 6.

If that pattern appears, the run is not using the baseline configuration. Check
that the full block from `dfgmres_deflation_eps` through
`hypre_use_interp_vectors` is present in `solver.inp`.
