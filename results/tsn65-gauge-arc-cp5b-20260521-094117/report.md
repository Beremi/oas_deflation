# CP5b TS-N65 Gauge Arc-Length Experiments

Generated: 2026-05-21T10:29:34+02:00

## Baseline Target

- Strict eight-quarter iteration sequence: `6,6,10,13,17,183,187,163`
- Strict total nonlinear rows: `585`
- First-quarter gauge target: `2.613423e-06`
- Full gauge target: `1.519697e-05`
- First-quarter internal-control calibration scale: `1.20576597808`
- Full internal-control calibration scale from CP3: `1.23219208515`

## Result Table

| variant | exit | verdict | steps | rows | duration | warnings | NaNs | fallback | cutbacks | final time | lambda | radius | gauge | target | gauge rel err |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CP5b-gauge-first-quarter | 0 | wrong_gauge | 1 | 25 | 00:04:07.195 | 0 | 0 | 0 | 0 | 0.00125 | 0.001035815 | 0.001 | 2.613423e-06 | 2.613423e-06 | 0.1706516702424368 |
| CP5b-gauge-first-quarter-calibrated | 0 | worsens | 1 | 25 | 00:04:07.583 | 0 | 0 | 0 | 0 | 0.00125 | 0.00124895 | 0.001 | 3.151177e-06 | 3.15117653973e-06 | 0.0 |
| CP5b-gauge-first-quarter-calibrated-line-search | 0 | worsens | 1 | 25 | 00:04:25.795 | 0 | 0 | 0 | 0 | 0.00125 | 0.00124895 | 0.001 | 3.151177e-06 | 3.15117653973e-06 | 0.0 |
| CP5b-gauge-first-quarter-line-search | 0 | wrong_gauge | 1 | 25 | 00:04:25.720 | 0 | 0 | 0 | 0 | 0.00125 | 0.001035815 | 0.001 | 2.613423e-06 | 2.613423e-06 | 0.1706516702424368 |
| CP5b-gauge-full-calibrated | 130 | failed | 7 | 235 |  | 0 | 0 | 0 | 0 | 0.0075 | 0.007355714 | 0.001 | 1.638489e-05 | 1.87255861523e-05 | 0.24681696417114724 |

## Notes

- `arc_length_constraint gauge` uses the CP3 `v01` displacement gauge as the continuation constraint.
- The runner appends a local target function to `functions.inp` and writes the actual zero-based `ic_function` index into each run's `solver.inp`.
- Uncalibrated runs target the internal IDC gauge directly. Calibrated runs scale that internal target because the exported `LD.out` interpolation gauge differs from the internal IDC coordinate selection.
- A result is promotable only if it reaches the target gauge with no fallback acceptance or NaNs and fewer rows or less wall time than the matching strict baseline target.
- `CP5b-gauge-full-calibrated` was manually stopped after step 7 started drifting away from convergence. Last parsed row: step 7, iteration 89, residual `3.963647e-3`, displacement `2.671126e-3`, energy `4.556947e-3`.
- The full calibrated run uses a linear gauge target to the final value, so it is a continuation experiment to the same final gauge rather than an exact strict-baseline loading-history replay.

## Files

- `cp5b_gauge_arc_length.tsv`: compact aggregate table.
- `runs/<variant>/solver.out`: ignored full solver log.
- `runs/<variant>/arc_length.tsv`: ignored arc-length trace with gauge columns.
- `runs/<variant>/results/LD.out`: ignored load/displacement gauge output.
