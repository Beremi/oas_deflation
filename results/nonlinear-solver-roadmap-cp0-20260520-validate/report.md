# CP0 Nonlinear Solver Roadmap Harness

Generated: 2026-05-20T15:25:31+02:00

## Result

CP0 is complete. The new parser reproduces the strict TS-N65 baseline from the
existing local full-deck replication log.

- Parsed steps: `8`
- Parsed nonlinear row sequence: `6,6,10,13,17,183,187,163`
- Parsed total nonlinear rows: `585`
- Target total nonlinear rows: `585`
- Sequence match: `yes`
- OAS total duration in parsed log: `00:22:47.554`
- End of calculation found: `yes`
- Warnings: `0`
- NaN lines: `0`
- Cutbacks: `0`
- Fallback accepts: `0`

## Code State

```text
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit_at_validation=71020a64685fafac36865278c404688b8a139903
baseline_run_commit=16c60d903808d69c7a06273d9f6d445f3603bae9
executable_for_baseline_run=/tmp/oas_tsn65_full_baseline_build/bin/OAS
```

## Commands

Parser validation:

```bash
scripts/summarize_oas_run.py \
  results/tsn65-full-baseline-replication-20260520-123252/runs/current-branch-full-baseline/solver.out \
  --baseline tsn65 --format text
```

Generated machine-readable files:

```bash
scripts/summarize_oas_run.py \
  results/tsn65-full-baseline-replication-20260520-123252/runs/current-branch-full-baseline/solver.out \
  --baseline tsn65 --format tsv \
  > results/nonlinear-solver-roadmap-cp0-20260520-validate/step_summary.tsv

scripts/summarize_oas_run.py \
  results/tsn65-full-baseline-replication-20260520-123252/runs/current-branch-full-baseline/solver.out \
  --baseline tsn65 --format json \
  > results/nonlinear-solver-roadmap-cp0-20260520-validate/summary.json
```

## Step Table

| step | time | dt | rows | target | match | duration | conv | residual | displacement | energy |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `0.00125` | `0.00125` | 6 | 6 | yes | `00:00:27.025` | 1 | `7.705844e-06` | `6.743608e-04` | `7.134058e-05` |
| 2 | `0.0025` | `0.00125` | 6 | 6 | yes | `00:00:18.418` | 1 | `4.198577e-04` | `1.021696e-04` | `6.255293e-04` |
| 3 | `0.00375` | `0.00125` | 10 | 10 | yes | `00:00:22.178` | 1 | `7.308889e-04` | `1.346596e-04` | `9.259996e-04` |
| 4 | `0.005` | `0.00125` | 13 | 13 | yes | `00:00:37.974` | 1 | `6.064366e-04` | `1.451312e-04` | `8.371949e-04` |
| 5 | `0.00625` | `0.00125` | 17 | 17 | yes | `00:00:41.066` | 1 | `6.189476e-04` | `1.136484e-04` | `9.322705e-04` |
| 6 | `0.0075` | `0.00125` | 183 | 183 | yes | `00:06:40.024` | 1 | `6.428869e-04` | `8.777311e-05` | `9.986343e-04` |
| 7 | `0.00875` | `0.00125` | 187 | 187 | yes | `00:06:45.088` | 1 | `6.721379e-04` | `1.022238e-04` | `9.840756e-04` |
| 8 | `0.01` | `0.00125` | 163 | 163 | yes | `00:05:50.751` | 1 | `5.838771e-04` | `1.075313e-04` | `9.817164e-04` |

## Files

- `docs/nonlinear_solver_roadmap.md`: living checkpoint roadmap.
- `scripts/summarize_oas_run.py`: reusable OAS log parser.
- `results/nonlinear-solver-roadmap-cp0-20260520-validate/step_summary.tsv`: parsed baseline step table.
- `results/nonlinear-solver-roadmap-cp0-20260520-validate/summary.json`: full parsed baseline summary.

## Verdict

CP0 passes. The parser and roadmap are ready to support CP1 sweep reporting.
