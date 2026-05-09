# hypre follow-up on the hard TS-N65 replay system

Date: 2026-05-09

Replay system:

- `results/tsn65-amgcl-failed-rhs-replaydump-20260509-j16/dfgmres-amgcl-hybrid6-failed-rhs-dump/linear_replay_failed/solve_42`
- Active full-coordinate system: `587190 x 587190`, `49449226` nonzeros.
- Symmetry probe: `1.36e-15`, so PCG was also tested as a diagnostic.

Relevant hypre documentation:

- BoomerAMG systems AMG uses `HYPRE_BoomerAMGSetNumFunctions`, `HYPRE_BoomerAMGSetDofFunc`, and `HYPRE_BoomerAMGSetNodal`.
- hypre documents rigid-body interpolation vectors through `HYPRE_BoomerAMGSetInterpVectors` and `HYPRE_BoomerAMGSetInterpVecVariant`, but only in the nodal systems-AMG path.
- The local hypre header documents `HYPRE_BoomerAMGSetNodalDiag`: `1` replaces the nodal diagonal by the negative off-diagonal row sum, and `2` inverts nodal diagonal signs.

## What was added for testing

- Added standalone replay options:
  - `--nullspace rbm6-normalized|rbm6-physical|rbm6-scale|rbm6-rot-neg|rbm6-rot-zero|transrot6`
  - `--nullspace-coordinate-scale`
  - `--interp-vec-qmax`
  - `--interp-vec-abs-qtrunc`
  - `--nodal-diag`
  - `--keep-same-sign`
- Added production OAS option:
  - `hypre_nodal_diag`
  - default is `0`, preserving existing behavior.

## Stable full-coordinate variants

These are all on the same hard RHS replay. Times are seconds.

| variant | threads | precond setup | AMG apply mean | Krylov solve | Krylov iters | true residual | solution error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| current cheby, interp 6, strong 0.25 | 1 | 2.969 | 0.563 | 9.591 | 16 | 9.301e-02 | 2.741e-01 |
| current cheby, interp 6, strong 0.25 | 16 | 0.937 | 0.388 | 6.346 | 16 | 9.301e-02 | 2.741e-01 |
| `hypre_nodal_diag=2`, same settings | 1 | 0.532 | 0.389 | 6.668 | 16 | 8.700e-02 | 2.813e-01 |
| `hypre_nodal_diag=2`, same settings | 16 | 0.442 | 0.265 | 4.544 | 16 | 8.700e-02 | 2.813e-01 |
| l1-Jacobi, interp 17, strong 0.25 | 16 | 0.680 | 0.062 | 4.686 | 69 | 9.975e-02 | 2.617e-01 |
| `keep_same_sign=1`, current cheby | 16 | 0.807 | 0.411 | 6.875 | 16 | 7.035e-02 | 2.752e-01 |

Best stable setting from this batch:

- `hypre_nodal_diag=2`
- It keeps the same Krylov iteration count as the current cheby setup: `16`.
- At 16 threads, compared with the current cheby setting:
  - preconditioner setup: `0.937 -> 0.442 s` (`2.12x` faster)
  - one-cycle AMG apply: `0.388 -> 0.265 s` (`1.46x` faster)
  - standalone Krylov solve: `6.346 -> 4.544 s` (`1.40x` faster)
  - true residual improves slightly: `9.30e-02 -> 8.70e-02`

The thread scaling of `hypre_nodal_diag=2` is still not excellent by itself:

| metric | 1 thread | 16 threads | speedup |
| --- | ---: | ---: | ---: |
| preconditioner setup | 0.532 | 0.442 | 1.20x |
| AMG apply mean | 0.389 | 0.265 | 1.47x |
| Krylov solve | 6.668 | 4.544 | 1.47x |

However, it is materially faster than the previous 16-thread setting, and it does not increase iterations on this replay.

## PCG diagnostic

PCG is not better here even though the matrix is symmetric.

| variant | threads | precond setup | AMG apply mean | Krylov solve | iters | true residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PCG + `hypre_nodal_diag=2` | 1 | 0.533 | 0.389 | 7.394 | 17 | 9.595e-02 |
| PCG + `hypre_nodal_diag=2` | 16 | 0.448 | 0.273 | 4.965 | 17 | 9.595e-02 |

Recommendation: keep DFGMRES/FGMRES-style use for this path.

## Rejected variants

Compressed-free reduced representations failed on this hard RHS:

- `compressed-free-node-major`, l1-Jacobi/interp 6 or 17: `hypre error 256 during HYPRE_ParCSRFlexGMRESSolve`.
- `compressed-free-node-major` + PCG: `hypre error 256 during HYPRE_ParCSRPCGSolve`.

Rigid-body interpolation vectors are still not production-safe for hypre on this system:

- `rbm6-scale` with scales `1.0`, `1.25`, `1.5`, `2.0` either timed out in setup or failed `HYPRE_BoomerAMGSetup`.
- `interp_vec_qmax = 1,2,4` did not fix setup.
- `interp_vec_abs_qtrunc = 1e-3, 1e-2, 1e-1` did not fix setup.
- `nodal_diag = 1,2` combined with interpolation vectors did not fix setup.

This is different from AMGCL, where physical/scaled rigid-body modes improved convergence. For hypre on this matrix, the stable improvement is nodal diagonal handling, not interpolation-vector fitting.

## Recommendation

The standalone replay result was strong enough to justify a full OAS check, but the nonlinear DFGMRES validation rejected it.

## First-step OAS validation

Both runs use current code, 16 threads, TS-N65 first loading step, DFGMRES, `N=20`, and `tol=true_tol=1e-1`.

| metric | current `hypre_nodal_diag=0` | candidate `hypre_nodal_diag=2` |
| --- | ---: | ---: |
| wall time | 117.114 s | 211.660 s |
| nonlinear iterations | 19 | 29 |
| linear solves | 19 | 29 |
| DFGMRES outer iterations | 37 | 351 |
| median DFGMRES iterations | 1 | 12 |
| max DFGMRES iterations | 5 | 23 |
| preconditioner setup | 4.303 s | 3.299 s |
| preconditioner apply | 14.458 s | 94.166 s |
| DFGMRES solve total | 16.684 s | 103.504 s |
| final true linear residual | 6.881e-02 | 8.882e-02 |
| final nonlinear residual | 8.138e-04 | 5.787e-04 |

Artifacts:

- Baseline: `results/tsn65-dfgmres-hypre-nodaldiag0-firststep-20260509-j16/dfgmres-hypre-tol1em1-N20`
- Candidate: `results/tsn65-dfgmres-hypre-nodaldiag2-firststep-20260509-j16/dfgmres-hypre-tol1em1-N20`

The candidate lowered setup time, but it made the preconditioner much worse for DFGMRES: the first step went from `37` to `351` outer iterations. That is the decisive metric for this solver path.

## Final Recommendation

Do not use `hypre_nodal_diag 2` for the current DFGMRES+hypre production path. Keep the current production setting:

```text
hypre_nodal_diag 0
hypre_use_interp_vectors 0
```

The new `hypre_nodal_diag` option can remain available for diagnostics and future sweeps, but the first-step validation shows it should not become the default. The replay-only result was misleading because standalone hypre Krylov solve time did not predict the quality of one BoomerAMG cycle as a right preconditioner inside DFGMRES.
