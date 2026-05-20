# CP4 Arc-Length Prototype Report

Generated: 2026-05-21T01:23:12+02:00

## Code State

- Repository: `/home/beremi/repos/oas_deflation`
- Commit: `dd80cb8ea89b3f1ea73266ccdd2a9dd90e5f35cd`
- Git status at report time: `M  docs/nonlinear_solver_roadmap.md
A  results/arc-length-cp4-20260521-cp4/cp4_arc_length.tsv
A  results/arc-length-cp4-20260521-cp4/metadata.txt
A  results/arc-length-cp4-20260521-cp4/report.md
A  results/tsn65-cp4-baseline-gate-20260521-0058/cp1_sweep.tsv
A  results/tsn65-cp4-baseline-gate-20260521-0058/metadata.txt
A  results/tsn65-cp4-baseline-gate-20260521-0058/report.md
A  scripts/run_arc_length_checkpoint.py
M  src/solver/src/solver_implicit.cpp
M  src/solver/src/solver_implicit.h
?? results/tsn65-full-baseline-replication-20260520-123252/
?? results/tsn65-local-baseline-20260520-100925/`
- Executable: `/tmp/oas_tsn65_full_baseline_build/bin/OAS`

## Benchmark

Small benchmark source: `src/benchmark/Timoshenko_beam3D`, copied into each run directory.
The runner rewrites the load functions as proportional PWL loads so the prototype reference load is `f_ext(1)-f_ext(0)`.

## Results

| variant | verdict | steps | rows | final time/lambda | arc rows | exit | elapsed s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `CP4-legacy-load-default` | passed | 3 | 6 | 0.3 | 0 | 0 | 0.064 |
| `CP4-load-control-regression` | passed | 3 | 6 | 0.3 | 0 | 0 | 0.064 |
| `CP4-arc-length-prototype` | passed | 3 | 6 | 0.3 | 6 | 0 | 0.064 |
| `CP4-arc-length-line-search` | passed | 4 | 18 | 0.3 | 18 | 0 | 0.064 |

## Checkpoint Notes

- `CP4-legacy-load-default` is the no-keyword legacy load-control run.
- `CP4-load-control-regression` exercises the explicit `nonlinear_control load` path.
- `CP4-arc-length-prototype` exercises the spherical Schur-complement predictor/corrector path.
- `CP4-arc-length-line-search` runs the same prototype with snapshot rollback and arc-length trial line-search logging enabled.
- Large run folders remain under `results/*/runs/` and are not intended for commit.
