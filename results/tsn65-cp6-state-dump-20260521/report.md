# TS-N65 Nonlinear Solver Sweep

Generated: 2026-05-21T23:49:59+02:00

## Metadata

```text
created_at=2026-05-21T23:27:00+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=5ba72d14ccfeb26240edea67d121ae1c292695a8
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
variants=D0-baseline-state-dump
```

## Result Table

| variant | exit | verdict | steps | rows | tail 6-8 | duration | warnings | NaNs | cutbacks | fallback | min alpha | max ls trials | IG accept | IG reject | IG min ratio |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| D0-baseline-state-dump | 0 | neutral | 8 | 585 | 533 | 00:22:56.082 | 0 | 0 | 0 | 0 |  | 0 | 0 | 0 |  |

## Baseline Gate

- Target sequence: `6,6,10,13,17,183,187,163`.
- Target total rows: `585`.
- Promote only variants that complete `8/8` and reduce total rows or the step 6-8 tail.

The CP6 state-dump instrumentation did not disturb the strict baseline:

```text
iteration_sequence=6,6,10,13,17,183,187,163
total_iters=585
fallback_acceptance_count=0
nan_lines=0
warnings=0
```

## State Dump Checkpoints

Step 6 accepted state:

```text
CSL_statuses=6171036
active_damage_growth_statuses=125947
damaged_statuses=148183
damage_mean=0.010117183416654532
damage_p99=0.48690121350286242
damage_p999=0.8559466170215374
damage_increment_p99=0.1270366795664456
damage_increment_p999=0.49600226479357001
material_status_hash_global=6622989350471731203
```

Step 7 accepted state:

```text
CSL_statuses=6171036
active_damage_growth_statuses=135027
damaged_statuses=184849
damage_mean=0.013442314591728025
damage_p99=0.58104234307705516
damage_p999=0.98331295283067233
damage_increment_p99=0.10831286900670489
damage_increment_p999=0.53224834736487392
material_status_hash_global=14127883338335681505
```

## Files

- `cp1_sweep.tsv`: machine-readable aggregate table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/summary.json`: ignored parsed full summary.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
- `runs/D0-baseline-state-dump/state/step_006_accepted_summary.json`: ignored step-6 state summary.
- `runs/D0-baseline-state-dump/state/step_007_accepted_summary.json`: ignored step-7 state summary.
