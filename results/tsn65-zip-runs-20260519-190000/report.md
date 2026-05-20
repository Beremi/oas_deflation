# TS-N65 Zip Benchmark Runs

Date: 2026-05-19
Branch: `nonlinear_solver_testing`
Build: `/tmp/oas_deflation_build_material_audit`

## Scope

The newly added benchmark archives were inspected and extracted into:

```text
results/tsn65-zip-runs-20260519-190000/input/
```

Archive/input size summary:

| item | size |
|---|---:|
| `TS-N_65.zip` | 110M |
| `Dogbone.zip` | 1.5M |
| extracted `TS-N_65/` | 925M |
| extracted `Dogbone/` | 11M |
| generated TS-N65 run directories | 3.7G |

`TS-N_65` contains the full target deck:

```text
auxnodes.inp
bc.inp
constraint.inp
elems.inp
exporters.inp
functions.inp
gov_nodes.inp
master.inp
materials.inp
nodes.inp
regions.inp
solver.inp
```

`Dogbone` was inspected and extracted. It includes a smaller deck plus previous `results/` outputs. This checkpoint ran TS-N65 because that is the current nonlinear-solver target benchmark.

## Original Deck Attempt

Run directory:

```text
results/tsn65-zip-runs-20260519-190000/runs/original/
```

The deck loaded successfully:

```text
DoF: 567923
constraint DoF: 19254
directly prescribed DoF: 13
elements: 514253
aux nodes: 18513108
```

The original `solver.inp` requests:

```text
solver_type PardisoLDLT
```

The local executable aborted before stepping:

```text
Solver type PardisoLDLT is not implemented
```

## Strict Tangent Hard-State Element Diagnostic

Run directory:

```text
results/tsn65-zip-runs-20260519-190000/runs/strict_tangent_element_top_step6_it10/
```

This used the strict quarter-step TS-N65 settings with `DeflatedFGMRES`, `stiff_matrix_type tangent`, `stiffness_matrix_iter_update 10`, and stopped after the hard-state diagnostic at step 6, iteration 10.

Command shape:

```bash
OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores \
  timeout 2400 /tmp/oas_deflation_build_material_audit/bin/OAS -j 16 master.inp
```

The diagnostic run reached step 6, iteration 10 in `00:05:24.737`.

Global Newton-direction check at the hard state:

| step | iteration | eps | direction | relative error | cosine | `||Kp||` | `||fd||` |
|---:|---:|---:|---|---:|---:|---:|---:|
| 6 | 10 | 1e-6 | newton | 5.320964 | 0.707442 | 1.186577e7 | 1.983825e6 |

Top element-local rows all identify the same material class behind the former plain `LDPMTetra` finding:

```text
element_name: LDPMTetra
material_id: 0
material_name: CSL material
material_status_names: CSL mat. statusx12
matrix_type: tangent
```

Top-50 element-local summary:

| matrix type | rows | mean relative error | max relative error | mean cosine | mean mismatch norm | max mismatch norm |
|---|---:|---:|---:|---:|---:|---:|
| `tangent` | 50 | 3.75268 | 7.00578 | -0.308147 | 213602 | 358831 |

First top row:

```text
element_id=300453 relative_error=3.70595 cosine=-0.160985 kp_norm=330284 fd_norm=96825.5 mismatch_norm=358831
```

## Consistent Element-Matrix Diagnostic At Same Hard State

Run directory:

```text
results/tsn65-zip-runs-20260519-190000/runs/strict_tangent_element_top_consistentmatrix_step6_it10/
```

This kept the nonlinear solve itself on the production `tangent` matrix, but ran the element-top diagnostic with:

```text
nonlinear_tangent_check_matrix_type consistent
```

The run reached the same hard state in `00:05:27.345`.

The top element rows are again:

```text
element_name: LDPMTetra
material_id: 0
material_name: CSL material
material_status_names: CSL mat. statusx12
matrix_type: consistent
```

Top-50 element-local summary:

| matrix type | rows | mean relative error | max relative error | mean cosine | mean mismatch norm | max mismatch norm |
|---|---:|---:|---:|---:|---:|---:|
| `consistent` | 50 | 0.133558 | 0.466171 | 0.986366 | 939.346 | 10748.7 |

First top row:

```text
element_id=74397 relative_error=0.223764 cosine=0.996811 kp_norm=37840.8 fd_norm=48036.1 mismatch_norm=10748.7
```

Interpretation:

- The TS-N65 hard-state top mismatch elements are `LDPMTetra` elements carrying CSL material/status objects.
- Replacing the element-local matrix used by the diagnostic with the new state-restored numerical `consistent` branch greatly reduces the top element mismatch and restores near-parallel finite-difference directionality.
- This confirms the local CSL active-damage tangent audit at the actual benchmark hard state.

## Global Consistent Matrix Smoke

Run directory:

```text
results/tsn65-zip-runs-20260519-190000/runs/strict_consistent_step1_smoke/
```

This attempted a one-quarter-step run with:

```text
stiff_matrix_type consistent
total_time 1.25e-3
```

The model loaded and assembled, but the first linear solve was not practical with the current `DeflatedFGMRES`/HYPRE setup:

```text
DeflatedFGMRES warning: performed 500 iterations and reached true relative residual 0.843808, required true tolerance is 0.1, active_basis_size=0
0   6.890472e-02            ---   2.543559e-03
```

The process was stopped after the first nonlinear row because it was already using about 27 GB RSS and had not completed the first step.

## Conclusions

1. The benchmark zips are valid local decks; TS-N65 is the large target deck, while Dogbone is a small deck with existing archived outputs.
2. The original TS-N65 solver settings cannot be run as-is in this build because `PardisoLDLT` is unavailable.
3. The strict TS-N65 hard-state rerun now attributes the previous top `LDPMTetra` mismatch rows to `CSL material` with `CSL mat. statusx12`.
4. The new numerical `consistent` branch fixes the element-local hard-state mismatch strongly enough to make the CSL active damage tangent the current primary culprit.
5. A production global `stiff_matrix_type consistent` run is not yet viable with the current linear solver/preconditioner path, despite the element-local diagnostic improvement.

## Next Steps

1. Keep the default `tangent` production path unchanged unless a later strict phase benchmark proves an improvement.
2. The follow-up analytical CSL active-damage tangent was derived, tested, and then archived as `archived_csl_damage_tangent` because it improved diagnostics but did not close TS-N65.
3. Treat `consistent` and `archived_csl_damage_tangent` as diagnostic/reference branches for the next path-following, tangent-limiting, or constitutive-regularization experiment.
