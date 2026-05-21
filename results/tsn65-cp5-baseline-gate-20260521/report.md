# CP1 TS-N65 Globalization Sweep

Generated: 2026-05-21T07:51:01+02:00

## Metadata

```text
created_at=2026-05-21T07:27:50+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=718f2ec94ef46a02905b19b463560a0248713088
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
variants=G1-baseline
```

## Result Table

| variant | exit | verdict | steps | rows | tail 6-8 | duration | warnings | NaNs | cutbacks | fallback | min alpha | max ls trials |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G1-baseline | 0 | neutral | 8 | 585 | 533 | 00:23:07.903 | 0 | 0 | 0 | 0 |  | 0 |

## Baseline Gate

- Target sequence: `6,6,10,13,17,183,187,163`.
- Target total rows: `585`.
- Promote only variants that complete `8/8` and reduce total rows or the step 6-8 tail.

## Files

- `cp5_baseline_gate.tsv`: machine-readable aggregate table.
- `cp1_sweep.tsv`: original aggregate filename emitted by the reused CP1 sweep runner.
- `baseline_step_summary.tsv`: compact strict baseline step table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/summary.json`: ignored parsed full summary.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
