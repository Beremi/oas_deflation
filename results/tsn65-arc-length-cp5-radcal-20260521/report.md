# CP5 TS-N65 Arc-Length Radius Calibration

Generated after manual stop.

## Result Table

| variant | status | steps | rows | final lambda | note |
| --- | --- | ---: | ---: | ---: | --- |
| CP5-arc-fd-r1p35e-2 | manual stop | 3 | 39 | 0.001196307 | not promotable: already 39 rows before first baseline quarter-step was stably reached; strict baseline step 1 is 6 rows |

## Interpretation

- Radius `1.35e-2` calibrates the first predictor to `lambda=0.00125`, but the corrector oscillates and shrinks the effective arc progress.
- The run reached only a partial early-state trace and was stopped before completion because it was already slower than the strict baseline for the same physical target.
- Compact trace: `arc_fd_r1p35e-2_trace.tsv`.
