# TS_N65 First-Step Thread Comparison

## Runs

All three runs use the current native `DeflatedFGMRES` + hypre preconditioner path with libRSB matvec enabled where available.

- `j=1`: `results/tsn65-dfgmres-librsb-firststep-20260507-j1/dfgmres-hypre-tol1em1-N20`
- `j=8`: `results/tsn65-dfgmres-librsb-firststep-20260507-j8/dfgmres-hypre-tol1em1-N20`
- `j=16`: `results/tsn65-dfgmres-librsb-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

Settings:

- First TS_N65 load step only: `total_time=5.000000e-03`
- `dfgmres_preconditioner=hypre_boomeramg_apply`
- `dfgmres_tolerance=1e-1`
- `dfgmres_true_tolerance=1e-1`
- `dfgmres_deflation_vectors=20`

Note: with `j=16`, the current DFGMRES libRSB matvec still uses the production cap of `librsb_threads=8`.

## Outcome

| Metric | j=1 | j=8 | j=16 | j1/j8 | j1/j16 |
|---|---:|---:|---:|---:|---:|
| Status | partial | partial | partial |  |  |
| Steps | 1 | 1 | 1 | 1.00x | 1.00x |
| Step elapsed from `solver.out` | 201.960 s | 124.080 s | 120.015 s | 1.63x | 1.68x |
| Wall incl. campaign overhead | 205.657 s | 127.828 s | 123.670 s | 1.61x | 1.66x |
| Model init total | 75.302 s | 70.832 s | 68.598 s | 1.06x | 1.10x |
| Model solve total / first step | 117.815 s | 44.689 s | 42.343 s | 2.64x | 2.78x |
| Nonlinear iterations | 19 | 19 | 19 | 1.00x | 1.00x |
| Linear solves | 19 | 19 | 19 | 1.00x | 1.00x |
| DFGMRES outer iterations | 37 | 37 | 37 | 1.00x | 1.00x |
| Max outer iters per solve | 4 | 4 | 4 | 1.00x | 1.00x |
| Final true relative residual | 0.0622639253 | 0.0622639253 | 0.0622639253 | 1.00x | 1.00x |
| Final nonlinear residual | 0.000787565 | 0.000787565 | 0.000787565 | 1.00x | 1.00x |

The important correctness signal is clean: increasing threads did not change nonlinear iteration count, linear solve count, DFGMRES outer iterations, or final residuals.

## Phase Timing

| Phase | j=1 | j=8 | j=16 | j1/j8 | j1/j16 |
|---|---:|---:|---:|---:|---:|
| Init total | 75.302 s | 70.832 s | 68.598 s | 1.06x | 1.10x |
| Init elements | 33.154 s | 32.591 s | 32.827 s | 1.02x | 1.01x |
| Init solver | 28.930 s | 23.757 s | 23.231 s | 1.22x | 1.25x |
| Solve step total | 117.815 s | 44.689 s | 42.343 s | 2.64x | 2.78x |
| Forces total active | 83.992 s | 20.049 s | 17.770 s | 4.19x | 4.73x |
| Integrate internal forces active | 81.648 s | 17.818 s | 15.499 s | 4.58x | 5.27x |
| Strain/stress kernel | 47.612 s | 9.721 s | 8.403 s | 4.90x | 5.67x |
| Element force/scatter kernel | 32.047 s | 6.060 s | 5.074 s | 5.29x | 6.32x |
| Simplex volumetric strains | 2.229 s | 2.127 s | 2.160 s | 1.05x | 1.03x |
| Stiffness update total | 17.901 s | 10.170 s | 9.891 s | 1.76x | 1.81x |
| Matrix factorize total | 14.526 s | 14.398 s | 14.429 s | 1.01x | 1.01x |
| Constraint transform | 2.597 s | 2.608 s | 2.570 s | 1.00x | 1.01x |
| Linear solve total | 15.942 s | 15.350 s | 15.318 s | 1.04x | 1.04x |
| Preconditioner setup | 10.164 s | 10.109 s | 10.044 s | 1.01x | 1.01x |
| Preconditioner apply | 12.681 s | 12.778 s | 12.828 s | 0.99x | 0.99x |
| DFGMRES matvec | 1.830 s | 1.075 s | 1.092 s | 1.70x | 1.68x |
| DFGMRES orthogonalization | 0.030 s | 0.015 s | 0.015 s | 2.08x | 2.04x |
| DFGMRES deflation projection/update | 0.102 s | 0.090 s | 0.089 s | 1.13x | 1.15x |
| DFGMRES least squares | 0.0009 s | 0.0010 s | 0.0010 s | 0.93x | 0.85x |

## Interpretation

The solve step itself scales much better than the full wall time: `117.8 s -> 42.3 s`, or `2.78x` from `j=1` to `j=16`. Full wall time is only `1.66x` faster because initialization is still a large mostly serial block: about `68.6 s` even at `j=16`.

The major OpenMP win is still force integration:

- `forces.integrate_internal_forces`: `81.65 s -> 15.50 s`, `5.27x`
- strain/stress work: `47.61 s -> 8.40 s`, `5.67x`
- element force/scatter work: `32.05 s -> 5.07 s`, `6.32x`

Going from `j=8` to `j=16` gives only modest additional improvement:

- step elapsed improves `124.08 s -> 120.01 s`
- solve step improves `44.69 s -> 42.34 s`
- force integration improves `17.82 s -> 15.50 s`
- linear solve is effectively flat

The remaining bottlenecks at `j=16` are mostly serial or weakly scaling: initialization, matrix factorization/preconditioner setup/apply, simplex volumetric strains, and stiffness update.
