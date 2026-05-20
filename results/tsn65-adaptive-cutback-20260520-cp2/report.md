# CP2 TS-N65 Adaptive Cutback/Stagnation

Generated: 2026-05-20T18:17:19+02:00

## Metadata

```text
created_at=2026-05-20T18:17:19+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=acbb5f95e20a7d7da9eef53e1c250180b795907e
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
variants=G2-stagnation-only,G2-backtracking-actual,G2-backtracking-frozen-actual,G2-backtracking-actual-adaptive-K,G2-fixed-damping-05
```

## Adaptive Deck Override

```text
max_iterations 40
min_time_step 1.25e-5
max_time_step 1.25e-3
critical_step_decrease 0.5
step_decrease 0.5
step_increase 1.2
enlargeIt 6
shortenIt 15
limit_tolerance 0
nonlinear_line_search_cutback_on_fail 1
nonlinear_stagnation_cutback 1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25
```

## Result Table

| variant | exit | verdict | steps | accepted | rows | tail 6-8 | min dt | duration | warnings | NaNs | cutbacks | restarts | shortens | fallback | min alpha | max ls trials |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G2-stagnation-only | 1 | failed | 11 | 10 | 183 | 61 | 1.25e-05 | 00:08:54.383 | 1 | 0 | 5 | 5 | 2 | 0 | 1 | 0 |
| G2-backtracking-actual | 1 | failed | 9 | 8 | 168 | 80 | 1.25e-05 | 00:17:18.914 | 3 | 0 | 5 | 5 | 2 | 0 | 0.125 | 4 |

## Checkpoint Verdict

CP2 did not pass for the executed variants. The runs avoided NaNs and fallback acceptance, but reached `min_time_step` and failed before `total_time=1e-2`.

Furthest accepted physical time across executed variants: `0.0078125`. This is a useful negative checkpoint and points the next experiment toward CP3 indirect displacement control before implementing arc-length.

## Per-Variant Step Tables

### G2-stagnation-only

| step | time | dt | rows | converged | residual | displacement | energy | min alpha | ls trials | cutback reasons |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 0.00125 | 0.00125 | 6 | 1 | 7.70584e-06 | 0.000674361 | 7.13406e-05 | 1 | 0 |  |
| 2 | 0.0025 | 0.00125 | 6 | 1 | 0.000419858 | 0.00010217 | 0.000625529 | 1 | 0 |  |
| 3 | 0.00375 | 0.00125 | 10 | 1 | 0.000730889 | 0.00013466 | 0.000926 | 1 | 0 |  |
| 4 | 0.005 | 0.00125 | 13 | 1 | 0.000606437 | 0.000145131 | 0.000837195 | 1 | 0 |  |
| 5 | 0.00625 | 0.00125 | 17 | 1 | 0.000618948 | 0.000113648 | 0.00093227 | 1 | 0 |  |
| 6 | 0.006875 | 0.000625 | 16 | 1 | 0.000607355 | 9.52886e-05 | 0.00094584 | 1 | 0 |  |
| 7 | 0.0071875 | 0.0003125 | 15 | 1 | 0.000582407 | 0.000115925 | 0.000996541 | 1 | 0 |  |
| 8 | 0.0075 | 0.0003125 | 30 | 1 | 0.000524742 | 0.000112247 | 0.000963856 | 1 | 0 | nonlinear_stagnation |
| 9 | 0.0075 | 0.00015625 | 39 | 1 | 0.000527003 | 0.000109691 | 0.000972715 | 1 | 0 | nonlinear_stagnation, nonlinear_stagnation |
| 10 | 0.007421875 | 3.90625e-05 | 13 | 1 | 0.000513589 | 0.000117245 | 0.000992428 | 1 | 0 | nonlinear_stagnation |
| 11 | 0.007421875 | 1.953125e-05 | 18 | 0 | 0.00049634 | 0.000119682 | 0.00101014 | 1 | 0 | nonlinear_stagnation |

### G2-backtracking-actual

| step | time | dt | rows | converged | residual | displacement | energy | min alpha | ls trials | cutback reasons |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 0.00125 | 0.00125 | 6 | 1 | 7.70584e-06 | 0.000674361 | 7.13406e-05 | 1 | 1 |  |
| 2 | 0.0025 | 0.00125 | 6 | 1 | 0.000419858 | 0.00010217 | 0.000625529 | 1 | 1 |  |
| 3 | 0.00375 | 0.00125 | 10 | 1 | 0.000730889 | 0.00013466 | 0.000926 | 1 | 1 |  |
| 4 | 0.005 | 0.00125 | 13 | 1 | 0.000606437 | 0.000145131 | 0.000837195 | 1 | 1 |  |
| 5 | 0.00625 | 0.00125 | 17 | 1 | 0.000618948 | 0.000113648 | 0.00093227 | 1 | 1 |  |
| 6 | 0.006875 | 0.000625 | 14 | 1 | 0.000910757 | 7.20066e-05 | 0.000958197 | 0.5 | 2 |  |
| 7 | 0.0075 | 0.000625 | 28 | 1 | 0.000992464 | 5.30719e-05 | 0.000870299 | 0.25 | 3 |  |
| 8 | 0.0078125 | 0.0003125 | 38 | 1 | 0.000758386 | 5.80621e-05 | 0.000772382 | 0.125 | 4 | line_search_failed, line_search_failed |
| 9 | 0.00765625 | 7.8125e-05 | 36 | 0 | 0.00077557 | 0.000191126 | 0.00158254 | 1 | 1 | nonlinear_stagnation, nonlinear_stagnation, nonlinear_stagnation |


## Baseline Comparison

- Strict target sequence: `6,6,10,13,17,183,187,163`.
- Strict target rows: `585`.
- CP2 pass requires reaching `total_time=1e-2` without fallback acceptance or NaNs.

## Files

- `cp2_adaptive.tsv`: machine-readable aggregate table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/summary.json`: ignored parsed full summary.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
