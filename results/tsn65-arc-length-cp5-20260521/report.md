# CP5 TS-N65 Arc-Length Experiments

Generated: 2026-05-21T07:27:10+02:00

## Baseline Target

- Final pseudo-time / load parameter: `0.00125`
- Strict baseline iteration sequence: `6,6,10,13,17,183,187,163`
- Strict baseline total rows: `585`

## Result Table

| variant | exit | verdict | steps | rows | duration | warnings | NaNs | fallback | cutbacks | final lambda | min radius | min dlambda |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CP5-arc-proportional-smoke | 0 | unsupported_zero_reference | 1 | 0 | 00:01:05.173 | 1 | 0 | 0 | 0 |  |  |  |
| CP5-arc-fd-r1p25e-3 | 1 | failed | 1 | 3 |  | 0 | 0 | 0 | 0 | 6.433438e-06 | 0.00125 | 5.549069e-06 |

## Notes

- `CP5-arc-proportional-smoke` checks the CP4 proportional-load arc-length path directly on TS-N65.
- `CP5-arc-fd-*` uses `arc_length_reference finite_difference`, which derives the reference direction from residual changes caused by the prescribed displacement/load parameter.
- A variant is only promotable if it reaches the same final parameter as the strict `0.01` baseline with no fallback acceptance or NaNs and fewer rows or less wall time.

## Files

- `cp5_arc_length.tsv`: compact aggregate table.
- `arc_fd_r1p25e-3_trace.tsv`: compact arc-length iteration trace for the finite-difference radius screen.
- `runs/<variant>/solver.out`: ignored full solver log.
- `runs/<variant>/arc_length.tsv`: ignored arc-length iteration trace.
