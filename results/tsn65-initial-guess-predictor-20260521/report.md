# CP1 TS-N65 Globalization Sweep

Generated: 2026-05-21T22:26:03+02:00

## Metadata

```text
created_at=2026-05-21T21:39:55+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=c0ca1e8b6479cec5a94605de1116f136ccea9762
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
variants=G1-baseline,P1-laststep-a025,P1-laststep-a050,P1-laststep-a075,P1-laststep-a100
```

## Result Table

| variant | exit | verdict | steps | rows | tail 6-8 | duration | warnings | NaNs | cutbacks | fallback | min alpha | max ls trials | IG accept | IG reject | IG min ratio |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| G1-baseline | 0 | neutral | 8 | 585 | 533 | 00:22:42.111 | 0 | 0 | 0 | 0 |  | 0 | 0 | 0 |  |
| P1-laststep-a025 | -15 | failed | 6 | 132 | 80 |  | 0 | 0 | 0 | 0 |  | 0 | 5 | 0 | 0.994465 |
| P1-laststep-a050 | -15 | failed | 6 | 102 | 50 |  | 0 | 0 | 0 | 0 |  | 0 | 5 | 0 | 0.983589 |
| P1-laststep-a075 | -15 | failed | 6 | 92 | 40 |  | 0 | 0 | 0 | 0 |  | 0 | 5 | 0 | 0.952484 |
| P1-laststep-a100 | -15 | failed | 6 | 125 | 73 |  | 0 | 0 | 0 | 0 |  | 0 | 5 | 0 | 0.000114001 |

## Interpretation

The strict baseline gate was rerun from the predictor-patched executable with
all new controls disabled. It reproduced the target exactly: `8/8` steps,
sequence `6,6,10,13,17,183,187,163`, total `585` nonlinear rows, no warnings,
no NaNs, no cutbacks, and no fallback acceptance.

All guarded `last_step` predictor variants were accepted on steps 2-6. They
improved the predictor guard merit at the beginning of the step, but none gave
a usable step-6 tail improvement. Each was manually stopped once its step-6
trajectory matched or slightly worsened the strict baseline.

Step-6 comparison:

| variant | row 0 residual | row 0 energy | later comparison |
| --- | ---: | ---: | --- |
| baseline | `2.548287e-02` | `1.001264e-01` | row 69 `3.831386e-03` / `7.130556e-03` |
| alpha `0.25` | `2.467055e-02` | `8.334072e-02` | row 69 `3.845413e-03` / `7.152297e-03` |
| alpha `0.50` | `2.407001e-02` | `6.225239e-02` | row 39 slightly worse than baseline |
| alpha `0.75` | `2.372770e-02` | `2.848583e-02` | row 39 slightly worse than baseline |
| alpha `1.00` | `2.547188e-02` | `5.114121e-02` | row 69 `3.857636e-03` / `7.181479e-03` |

Verdict: the missing speed lever is not a simple last-step free-DOF initial
guess. The predictor can improve the first residual/energy row, but the hard
late-step nonlinear tail quickly reappears.

## Baseline Gate

- Target sequence: `6,6,10,13,17,183,187,163`.
- Target total rows: `585`.
- Promote only variants that complete `8/8` and reduce total rows or the step 6-8 tail.

## Files

- `cp1_sweep.tsv`: machine-readable aggregate table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/summary.json`: ignored parsed full summary.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
