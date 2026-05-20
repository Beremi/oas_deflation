# CP1 TS-N65 Globalization Sweep

Generated: 2026-05-20T17:21:28+02:00

Updated: 2026-05-20T17:24:00+02:00

## Metadata

```text
created_at=2026-05-20T15:28:27+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=f1cae734dedf5a65e1b8a917404679e9fc90eaf3
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
variants=G1-baseline,G1-backtracking-frozen-actual,G1-backtracking-actual,G1-bisection-frozen-actual,G1-fixed-damping-05,G1-adaptive-damping,G1-trust-stepnorm,G1-line-search-adaptive-K
```

## Result Table

| variant | exit | verdict | steps | rows | tail 6-8 | duration | warnings | NaNs | cutbacks | fallback | min alpha | max ls trials |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G1-baseline | 0 | neutral | 8 | 585 | 533 | 00:22:58.080 | 0 | 0 | 0 | 0 |  | 0 |
| G1-backtracking-frozen-actual | 124 | failed | 8 | 495 | 443 |  | 0 | 0 | 0 | 0 | 0.0625 | 6 |
| G1-backtracking-actual | manual stop | partial worse | 7 | 161 | 109 partial |  | 0 | 0 | 0 | 0 | 0.125 | 6 |

## Findings So Far

- `G1-baseline` reproduced the strict baseline exactly on the current build:
  `8/8`, `585` rows, no warnings, no NaNs, no fallback acceptance.
- `G1-backtracking-frozen-actual` reduced step 6 from `183` rows to `63`, but
  step 7 grew from `187` rows to `270` and took `00:50:36.666`; step 8 was
  still unconverged when the 90-minute runner timeout fired.
- `G1-backtracking-actual` also reduced step 6 to `63` rows with lower step-6
  cost (`00:08:09.736` versus `00:11:33.801` for frozen/actual), but it showed
  the same step-7 fixed-load pathology. It was manually stopped at step 7 row
  `45` after errors rose to residual `6.728244e-03`, displacement
  `2.237895e-03`, and energy `1.310600e-02`.
- Neither tested backtracking mode is promotable from CP1. Both suggest the next
  useful experiment is CP2 adaptive cutback/stagnation, not more strict fixed
  stepping with the same line-search policy.

## Remaining CP1 Variants

These variants were not run in this partial CP1 pass:

- `G1-bisection-frozen-actual`
- `G1-fixed-damping-05`
- `G1-adaptive-damping`
- `G1-trust-stepnorm`
- `G1-line-search-adaptive-K`

They remain queued, but the first two line-search variants already show that
strict fixed-step globalization can trade fewer rows in step 6 for a much worse
step 7. Run the remaining variants only if a direct comparison is still needed
before CP2.

## Baseline Gate

- Target sequence: `6,6,10,13,17,183,187,163`.
- Target total rows: `585`.
- Promote only variants that complete `8/8` and reduce total rows or the step 6-8 tail.

## Files

- `cp1_sweep.tsv`: machine-readable aggregate table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/summary.json`: ignored parsed full summary.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
