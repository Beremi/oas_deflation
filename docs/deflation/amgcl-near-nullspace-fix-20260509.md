# AMGCL Near-Nullspace Diagnostic And Fix

Dumped system: `results/tsn65-amgcl-failed-rhs-replaydump-20260509-j16/dfgmres-amgcl-hybrid6-failed-rhs-dump/linear_replay_failed/solve_42`

The replay tests isolate the failed TS-N65 RHS from nonlinear force timing. All rows below use `full-coordinate-node-major`, `hybrid6`, `spai0`, FGMRES restart 80, one AMG V-cycle (`npre=1`, `npost=1`, `ncycle=1`), and 16 OpenMP threads unless noted.

![nullspace-residual-history.png](../../results/tsn65-amgcl-nullspace-variants-20260509/nullspace-residual-history.png)

## Nullspace Variant Sweep

| near-nullspace | setup s | solve s | iterations | final true residual | solution error |
| --- | ---: | ---: | ---: | ---: | ---: |
| `rbm6-physical` | 0.930 | 31.877 | 500 | 0.132618 | 0.197469 |
| `rbm6-amgcl-physical` | 0.976 | 32.103 | 500 | 0.132627 | 0.197468 |
| `rbm6-rot-neg` | 0.968 | 31.905 | 500 | 0.155833 | 0.246822 |
| `rbm6-normalized` | 0.961 | 32.033 | 500 | 0.168996 | 0.248497 |
| `rbm6-amgcl-order` | 0.929 | 32.122 | 500 | 0.168996 | 0.248499 |
| `rbm6` | 0.969 | 32.042 | 500 | 0.169 | 0.248546 |
| `rbm6-rot-zero` | 0.945 | 31.958 | 500 | 0.201871 | 0.375444 |

`rbm6-physical` is the same six particle rigid-body modes as production, but with centered physical coordinates. The old/current exported `rbm6` normalized coordinates by the model RMS radius before adding unit rotational DOF rows; that changes the actual six-DOF rigid-body vector directions, not just column scaling.

## Convergence-To-Tolerance Check

| near-nullspace | max iterations | solve s | iterations used | true residual |
| --- | ---: | ---: | ---: | ---: |
| `rbm6-physical` | 800 | 49.648 | 779 | 0.0999311 |
| `rbm6-physical` | 1000 | 49.760 | 779 | 0.0999311 |
| `rbm6-physical` | 1500 | 49.913 | 779 | 0.0999311 |

The physical-coordinate modes reach the `1e-1` target in 779 iterations on the isolated failed RHS. The old normalized six-DOF modes required 1352 iterations in the earlier 1500-iteration replay.

## Krylov Method Check

| method | solve s | iterations | true residual | solution error | interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| `cg` | 29.908 | 500 | 0.333954 | 0.263609 | stalls; matrix/preconditioner path is not CG-safe |
| `bicgstab` | 38.887 | 327 | 0.0994349 | 0.554376 | fast residual drop, but large solution error vs dumped reference |
| `fgmres` | 31.873 | 500 | 0.132618 | 0.197469 | robust but still slow at one V-cycle |

## Code Fix

Changed `SteadyStateLinearSolver::buildElasticDofMap()` so coordinate normalization is skipped when the mechanical block contains rotational DOFs (`blockSize > dim`). For displacement-only elasticity, normalization is still kept because it only rescales rotational columns. For particle six-DOF mechanics, it changes the relative magnitude of translational motion `omega x r` against rotational rows `omega`, so it corrupts the rigid-body mode shape.

## Remaining Risk

This fixes a real near-nullspace bug and cuts the failed RHS convergence from ~1352 to 779 FGMRES iterations in replay, but it does not make AMGCL competitive with the best hypre replay on this RHS. The next tuning knob, if AMGCL remains interesting, is a stronger preconditioner cycle or a BiCGSTAB-style outer method with nonlinear validation.
