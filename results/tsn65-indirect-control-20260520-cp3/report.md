# CP3 TS-N65 Indirect Displacement Control

Generated: 2026-05-20T21:16:09+02:00

## Metadata

```text
created_at=2026-05-20T21:16:09+02:00
repo=/home/beremi/repos/oas_deflation
branch=nonlinear_solver_testing
commit=573dd62ec53cf2601a746aa79825ef7b13414aed
executable=/tmp/oas_tsn65_full_baseline_build/bin/OAS
deck=/home/beremi/repos/oas_deflation/results/tsn65-zip-runs-20260519-190000/input/TS-N_65
target_v01_final=1.519697e-05
```

## Function Index Note

OAS function indices are zero-based. The handout block listed `ic_function 2`, but in the TS-N65 deck function index `2` is the steel material envelope, not the displacement ramp. CP3 appends a local IDC target function and uses `ic_function 4` inside each ignored run directory.

The appended target is:

```text
PWLFunction 2 0 0.01 0 1.519697e-05 # CP3 IDC v01 target
```

The final-v01 calibrated strict run scales that target by `1.23219208515`, derived from `target_v01 / CP3-idc-minus-strict_final_v01`.

## Baseline Gauges

| point | time | displ[m] | load[N] | v01 |
| --- | ---: | ---: | ---: | ---: |
| first strict step | 0.00125 | -6e-05 | -37401.27 | 2.613423e-06 |
| strict target | 0.01 | -0.00048 | -216962.5 | 1.519697e-05 |

## Result Table

| variant | exit | verdict | steps | rows | duration | warnings | NaNs | fallback | cutbacks | min dt | weights | sign score | final time | idc time | displ[m] | load[N] | v01 | v01 rel err |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CP3-indirect-control-benchmark | 1 | failed | 146 | 21910 | 00:02:00.750 | 1 | 0 | 138 | 0 | 0.0005 |  |  | 0.0778 | 0.07509961 |  |  |  |  |
| CP3-idc-minus-calibration | 0 | passed | 1 | 25 | 00:03:03.077 | 0 | 0 | 0 | 0 | 0.00125 | -1 1 | 3 | 0.00125 |  | -3.61514e-05 | -22597.26 | 1.575315e-06 | 0.89634 |
| CP3-idc-plus-calibration | 1 | failed | 1 | 0 |  | 0 | 0 | 0 | 0 | 0.00125 | 1 -1 | 0 |  |  |  |  |  |  |
| CP3-idc-minus-strict | 0 | failed | 8 | 193 | 00:16:06.195 | 0 | 0 | 0 | 0 | 0.00125 | -1 1 | 3 | 0.01 |  | -0.0003193611 | -176538.5 | 1.233328e-05 | 0.188438 |
| CP3-idc-minus-final-v01-calibrated-strict | -15 | failed | 7 | 218 |  | 0 | 0 | 0 | 0 | 0.00125 | -1 1 | 3 | 0.0075 |  | -0.0002900148 | -163956.9 | 1.144418e-05 | 0.246943 |

## Failed Or Incomplete Last Rows

| variant | exit | last step | last iter | residual | displacement | energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CP3-indirect-control-benchmark | 1 | 146 | 149 | 1.009874e-01 | 6.389661e-06 | 5.432050e-04 |
| CP3-idc-plus-calibration | 1 | 1 | None | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 |
| CP3-idc-minus-final-v01-calibrated-strict | -15 | 7 | 69 | 2.352053e-03 | 4.419994e-04 | 4.186416e-03 |

## Checkpoint Verdict

Sign calibration selected `minus` weights based on the first-step direction score.
CP3 did not pass for the strict TS-N65 IDC run(s). Preserve the run directories and use the compact table to decide whether to try adaptive IDC or CP4 arc-length.

## Files

- `cp3_indirect_control.tsv`: machine-readable aggregate table.
- `runs/<variant>/solver.out`: ignored full logs.
- `runs/<variant>/results/LD.out`: ignored gauge output.
- `runs/<variant>/step_summary.tsv`: ignored parsed step table.
