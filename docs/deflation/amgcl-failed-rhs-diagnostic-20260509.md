# AMGCL Failed-RHS Diagnostic, 2026-05-09

## Question

The direct AMGCL hybrid6 one-cycle preconditioner looked good on the first exported TS-N65 linear system, but failed inside the first loading step. The suspected failure modes were a bad near-nullspace map, DOF-order mismatch, or an active-matrix permutation bug.

## What Changed For Diagnostics

- Added `amgcl_check_matrix` diagnostics in `DeflatedFGMRES::factorize()`:
  - active backend matvec versus Eigen sparse matvec,
  - active near-nullspace row count, permutation mode, and `A * mode` scale,
  - one deterministic preconditioner probe `||A M b - b|| / ||b||`.
- Added `linear_solver_replay_start` so one late linear solve can be dumped without writing every previous matrix.

These diagnostics are inactive unless requested through solver input.

## OAS First-Step Diagnostics

Run:

`results/tsn65-amgcl-diagnostic-full-firststep-20260509-j16/dfgmres-amgcl-hybrid6-firststep-diagnostic`

Settings:

- `dfgmres_preconditioner=amgcl`
- `dfgmres_deflation_vectors=20`
- `amgcl_backend=hybrid`
- `amgcl_block_size=6`
- `amgcl_near_nullspace=1`
- `amgcl_npre=1`, `amgcl_npost=1`, `amgcl_ncycle=1`
- `dfgmres_elastic_reorder=2`

Key results:

| check | value |
| --- | ---: |
| active rows | 587190 |
| reduced rows | 567923 |
| lifted identity rows | 19267 |
| active ordering | coordinate node-major |
| near-nullspace modes | 6 |
| CRS/librsb matvec vs Eigen | `3.37e-16` relative difference |
| `A` action per unit generic vector norm | `2.633e9` |
| near-nullspace `A*mode` / reference action | `0.0065` to `0.0218` |
| preconditioner probe `||A M b - b|| / ||b||` | `5.22` |
| matrix factorizations in first step | 2 |
| linear solves in first step | 42 |
| successful linear solves | 41 |
| failed solve | solve 42 |
| failed solve DFGMRES iterations | 500 |
| failed solve true residual | `0.1216888` |

Interpretation:

- The active matrix permutation and matvec path are correct.
- The near-nullspace vector table has the expected size and is permuted consistently.
- The six rigid-body modes are not exact null modes because the constrained/lifted system has boundary constraints and identity rows, but they are still much lower energy than a generic vector.
- The matrix did not meaningfully change before the failure. The same active matrix handled early RHS vectors and then failed on a later RHS.

## Deflation Check

Run:

`results/tsn65-amgcl-no-deflation-firststep-20260509-j16/dfgmres-amgcl-hybrid6-nodefl-firststep`

Same AMGCL setup, but:

- `dfgmres_deflation_vectors=0`
- `dfgmres_collect_newton_steps=0`

Result:

| case | nonlinear iterations | linear solves | successes | max DFGMRES iters | final failed true residual |
| --- | ---: | ---: | ---: | ---: | ---: |
| AMGCL hybrid6 direct, N=20 | 41 | 42 | 41 | 500 | `0.1216888` |
| AMGCL hybrid6 direct, no deflation | 44 | 45 | 44 | 500 | `0.151146` |

Interpretation: deflation/recycling is not the root cause. Removing it makes the late solve no better.

## Failed RHS Replay

Dumped failed solve:

`results/tsn65-amgcl-failed-rhs-replaydump-20260509-j16/dfgmres-amgcl-hybrid6-failed-rhs-dump/linear_replay_failed/solve_42`

Summary:

| field | value |
| --- | ---: |
| step | 1 |
| nonlinear iteration | 41 |
| replay solve index | 42 |
| rows | 567923 |
| nnz | 49429959 |
| OAS solver success | 0 |
| OAS iterations | 500 |
| OAS true residual | `0.1216888` |

Standalone replay on this exact RHS:

| backend | representation | nullspace | solver | settings | iterations | true residual | solve time, j=16 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: |
| AMGCL hybrid6 | full-coordinate-node-major | RBM6 | FGMRES | `eps=0`, `1/1/1` | 500 | `0.206862` | 32.38 s |
| AMGCL hybrid6 | full-coordinate-node-major | off | FGMRES | `eps=0`, `1/1/1` | 500 | `0.182066` | 32.05 s |
| AMGCL hybrid6 | full-coordinate-node-major | RBM6 | FGMRES | `eps=0.25`, `1/1/1` | 500 | `0.158848` | 30.36 s |
| AMGCL hybrid6 | full-coordinate-node-major | RBM6 | FGMRES | `eps=0.25`, `1/1/1`, maxiter 1000 | 1000 | `0.111505` | 59.02 s |
| AMGCL hybrid6 | full-coordinate-node-major | RBM6 | FGMRES | `eps=0`, `3/3/1` | 500 | `0.164377` | 53.94 s |
| AMGCL hybrid6 | full-coordinate-node-major | RBM6 | FGMRES | `eps=0`, `1/1/2` | 500 | `0.178753` | 57.24 s |
| hypre BoomerAMG | full-coordinate-node-major | off | FGMRES | nodal 4, interp 6, strong 0.25, relax 16, sweeps 3 | 16 | `0.093010` | 6.76 s |

The previous "good" AMGCL replay result was from the first RHS, not this late RHS. The late RHS is the real failing case.

## Krylov Residual History

Residual-history plot:

`results/tsn65-amgcl-failed-rhs-krylov-history-20260509/krylov-history.md`

Key replay curves on the failed RHS:

| method | restart | first printed residual | last printed residual at 500 iters | behavior |
| --- | ---: | ---: | ---: | --- |
| AMGCL FGMRES hybrid6/RBM6 | 80 | `0.326509` | `0.169541` | monotone decreasing in printed samples, too slow |
| AMGCL FGMRES hybrid6/RBM6 | 30 | `0.326509` | `0.207129` | monotone decreasing in printed samples, slower than restart 80 |
| AMGCL CG hybrid6/RBM6 | - | `0.342840` | `0.351052` | residual jumps to `1.85` at iter 5 and oscillates |
| AMGCL CG scalar/reduced | - | `3.666568` | `652.997` | diverges badly |

One longer bounded check:

| method | restart | max iterations | reached iterations | final true residual | solve time, j=16 |
| --- | ---: | ---: | ---: | ---: | ---: |
| AMGCL FGMRES hybrid6/RBM6 | 80 | 1500 | 1352 | `0.0999394` | 86.23 s |

Interpretation:

- FGMRES is not decreasing and then increasing. It decreases steadily, but the AMGCL one-cycle preconditioner is too weak for this RHS at the `500` iteration cap.
- CG is not suitable with this preconditioner/operator combination. Its residual increases and oscillates, which points to a non-CG-compatible preconditioned operator. The matrix is symmetric by structure, but the AMGCL/SPAI0 hierarchy and lifted constrained system should not be assumed SPD in the PCG sense.
- I also fixed the replay benchmark to actually apply `--gmres-restart`; previous AMGCL FGMRES replay used AMGCL's default `M=30`, not the OAS restart of 80.

## Conclusion

I do not see evidence for a near-nullspace row mismatch or active DOF permutation bug. The failure reproduces in standalone replay on the dumped late RHS, and it also happens with nullspace disabled and with deflation disabled.

The concrete problem is that AMGCL hybrid6/SPAI0 one-cycle is too weak as a true DFGMRES preconditioner for this late nonlinear RHS. It can solve the first exported RHS, but not the hard RHS near nonlinear iteration 41 at the current 500-iteration cap. It does eventually cross `1e-1`, but only at about 1352 FGMRES iterations on the dumped RHS. Hypre on the same dumped RHS reaches the tolerance in 16 iterations.

Recommendation: do not switch production DFGMRES to AMGCL one-cycle hybrid6 for TS-N65 yet. Keep using hypre for this path, or use AMGCL only as an inner solve if we accept the much higher apply cost. The next useful AMGCL work would be testing a different relaxation/coarsening backend, not more near-nullspace permutation debugging.
