# AMGCL Near-Nullspace Enhancement For TS-N65 Particle System

Replay system: `results/tsn65-amgcl-failed-rhs-replaydump-20260509-j16/dfgmres-amgcl-hybrid6-failed-rhs-dump/linear_replay_failed/solve_42`

## Model-Specific Near-Nullspace

For a 6-DOF particle with translational unknown `u_i` and rotational unknown `theta_i`, a global rigid-body motion should be represented as

```text
u_i     = t + omega x (x_i - x0)
theta_i = omega
```

This is different from the standard displacement-only AMGCL/FEM elasticity example, where only the three displacement rows exist per node. OAS now keeps physical coordinates by default for rotational particle blocks. The diagnostic below asks whether a slightly modified translation/rotation ratio gives AMGCL a better coarse space on the damaged failed RHS.

![rbm6-scale-sweep.png](../../results/tsn65-amgcl-nullspace-enhanced-20260509/rbm6-scale-sweep.png)

## Six-Mode Scale Sweep

All rows use 16 threads, `full-coordinate-node-major`, AMGCL `hybrid6`, FGMRES restart 80, `eps_strong=0`, and one AMG V-cycle. Scale is the divisor used in the rotational displacement part `omega x ((x_i-x0)/scale)`; `scale=1` is the physical rigid-body mode and `scale=0.7296073` reproduces the old RMS-normalized OAS basis.

| scale | solve s | iterations | final true residual | note |
| ---: | ---: | ---: | ---: | --- |
| 0.5 | 32.000 | 500 | 0.171109 |  |
| 0.7296073 | 31.864 | 500 | 0.169 | old normalized scale |
| 0.85 | 31.853 | 500 | 0.16306 |  |
| 1 | 31.882 | 500 | 0.132618 | physical rigid-body scale |
| 1.1 | 31.659 | 500 | 0.106893 |  |
| 1.15 | 31.803 | 500 | 0.110264 |  |
| 1.2 | 31.732 | 500 | 0.0999973 |  |
| 1.25 | 31.846 | 499 | 0.0998527 | best tested crossing target |
| 1.3 | 31.719 | 500 | 0.105002 |  |
| 1.35 | 31.714 | 500 | 0.111224 |  |
| 1.4 | 31.775 | 500 | 0.116528 |  |
| 1.5 | 31.836 | 500 | 0.121923 |  |
| 2 | 31.903 | 500 | 0.136026 |  |
| 3 | 31.890 | 500 | 0.131858 |  |

Best tested six-mode basis: `scale=1.25`, which reached `9.985e-2` in 499 FGMRES iterations. This is an empirical near-nullspace weighting, not the exact physical rigid-body nullspace. It improves AMGCL on the failed RHS but is still much weaker than the best hypre replay.

## Augmented Mode Diagnostics

| variant | backend | block size option | solve/setup result | final residual | interpretation |
| --- | --- | ---: | --- | ---: | --- |
| `transrot6` | `hybrid6` | 6 | ok; iters 500 | 0.139995 | pure translation/rotation fields are worse than scaled rigid modes |
| `rbm9-physical` | `hybrid6` | 6 | Matrix size is not divisible by block size!; iters -1 | - | more than 6 modes conflicts with block grouping/hybrid hierarchy |
| `rbm9-physical` | `scalar` | 6 | Matrix size should be divisible by block_size; iters -1 | - | more than 6 modes conflicts with block grouping/hybrid hierarchy |
| `rbm6-scale` | `scalar` | 6 | ok; iters 499 | 0.0998527 | same convergence as hybrid6 but slower apply path |
| `rbm9-physical` | `scalar` | 1 | ok; iters 300 | 0.156659 | runs only with scalar/block-size workaround and is slower/weaker |
| `rbm9-physical` | `scalar` | 3 | ok; iters 300 | 0.201873 | runs only with scalar/block-size workaround and is slower/weaker |

## Code Changes

- Production default remains the physical six-DOF particle rigid-body basis for rotational blocks.
- Added optional solver-file control: `elastic_near_nullspace_coordinate_scale 1.25` (aliases: `amg_near_nullspace_coordinate_scale`, `amgcl_near_nullspace_coordinate_scale`) to reproduce the best replay weighting in full OAS runs.
- Replay benchmark now supports `rbm6-scale`, `transrot6`, and `rbm9-physical` variants.

## Recommendation

Use physical six rigid-particle modes as the correctness default. If we run another full TS-N65 AMGCL nonlinear test, use `elastic_near_nullspace_coordinate_scale 1.25` as the candidate enhanced near-nullspace. Do not use the 9-mode augmented basis with AMGCL `hybrid6`; it is structurally incompatible because the coarse level no longer aligns with 6x6 blocks.
