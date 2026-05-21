# CP1 TS-N65 Globalization Sweep

Generated: 2026-05-21T21:39:38+02:00

## Metadata

```text
created_at=2026-05-21T20:46:25+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=c0ca1e8b6479cec5a94605de1116f136ccea9762
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
variants=G1-baseline,G1-firstit-elastic,G1-firstit-secant,G1-firstit-archived-csl
```

## Result Table

| variant | exit | verdict | steps | rows | tail 6-8 | duration | warnings | NaNs | cutbacks | fallback | min alpha | max ls trials |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G1-baseline | 0 | neutral | 8 | 585 | 533 | 00:23:03.491 | 0 | 0 | 0 | 0 |  | 0 |
| G1-firstit-elastic | -15 | failed | 6 | 188 | 130 |  | 0 | 0 | 0 | 0 |  | 0 |
| G1-firstit-secant | -15 | failed | 7 | 367 | 315 |  | 0 | 0 | 0 | 0 |  | 0 |
| G1-firstit-archived-csl | -15 | failed | 6 | 132 | 80 |  | 0 | 0 | 0 | 0 |  | 0 |

## Interpretation

These variants screen whether the first Newton correction direction alone can
improve the late TS-N65 tail.

- `first_iteration_stiff_matrix_type elastic` changed the early cheap steps
  from baseline `6,6,10,13,17` to `6,6,11,15,20` and was stopped in step 6.
- `first_iteration_stiff_matrix_type secant` matched baseline through step 6,
  then entered the same step-7 degradation pattern seen in earlier failed
  globalization runs. It was stopped at step 7 with no tail-speed signal.
- `first_iteration_stiff_matrix_type archived_csl_damage_tangent` matched
  baseline through step 5, then worsened inside step 6 and was stopped.

Verdict: the first-correction matrix-type branch is not a useful speed path for
the strict eight-quarter TS-N65 baseline.

## Baseline Gate

- Target sequence: `6,6,10,13,17,183,187,163`.
- Target total rows: `585`.
- Promote only variants that complete `8/8` and reduce total rows or the step 6-8 tail.

## Files

- `cp1_sweep.tsv`: machine-readable aggregate table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/summary.json`: ignored parsed full summary.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
