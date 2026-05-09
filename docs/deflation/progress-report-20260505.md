# OAS Deflation and Iterative Solver Progress Report

Date: 2026-05-05

This report summarizes the local OAS deflation workspace from the bootstrap commit through the latest adaptive-tolerance work. It is based on the committed implementation notes, solver/profiling code, and ignored local experiment reports under `results/`.

## Executive Summary

We started from a direct-solver OAS baseline where `PardisoLDLT` was the only credible large-case solver for TS-N_65, and where the immediate risk was not knowing where time was spent or how repeated nonlinear linear systems evolved.

The work moved in four stages:

1. Build reliable linear-solve and runtime profiling around OAS.
2. Add and tune AMG-style iterative solvers on the small Dogbone case.
3. Gate those solvers on the large TS-N_65 case.
4. Add native `DeflatedFGMRES` with Newton-increment recycling and test whether deflation improves the repeated-system workload.

The main result is that `HypreBoomerAMGCG` is currently the most practical iterative path on TS-N_65. It completed the first two TS-N_65 load steps and can trade accuracy for speed by relaxing the linear tolerance. `DeflatedFGMRES+hypre` also works and deflation is active; at strict `1e-6` tolerance it beat hypre-CG `1e-6` on the two-step run, but it did not beat the practical loose hypre-CG runs, and loose DFGMRES tolerances can converge globally while producing a different local displacement/damage path. That makes loose DFGMRES unsafe for production fracture results.

## Starting Point

The initial deflation note identified the relevant OAS hooks:

- `src/solver/src/linalg.h` and `src/solver/src/linalg.cpp` hold linear solver wrappers.
- `src/solver/src/solver_implicit.cpp` owns the implicit Newton loops and calls `linalgsolver->solve(ddr, f)`.
- `computeKeff()` and `factorizeLinearSystem()` create and update the effective stiffness matrix.

The planned deflation design was to store a small basis of accepted Newton increments, compute `AZ` and `Z^T A Z` on factorization, and use that basis to project repeated Krylov solves. The early estimate was that 20 vectors and their images for TS-N_65 would cost roughly 183 MiB, which is small compared with the reported 20 GB direct-solver memory footprint.

Initial benchmark roles:

| case | role | size/notes |
| --- | --- | --- |
| Dogbone | fast development and correctness gate | about 6.6k reduced DOF, 100 load steps |
| TS-N_65 | hard target | about 568k reduced rows in profiles, about 49.4M nnz, reported multi-day Pardiso runtime |

## Implementation Added

The repository now contains the following local solver infrastructure:

| area | implementation |
| --- | --- |
| Profiling | `linear_solver_profile`, matrix value/hash tracking, per-solve residual/iteration/status/timing, runtime phase profiling |
| Analysis scripts | `scripts/analyze-linear-profile.py`, `scripts/render-vtu-state.py`, profile runners and campaign scripts |
| AMGCL | `AmgclCG` / `AmgclCGElastic`, elastic lifting, near-nullspace modes, hybrid backend, diagonal scaling, true physical residual checks |
| hypre | `HypreBoomerAMGCG`, tuned BoomerAMG-CG options, persistent setup reuse, true residual reporting |
| Deflation | `DeflatedFGMRES`, accepted Newton-increment collection, A-orthogonal sliding basis, basis rebuild on matrix changes, orthogonality diagnostics |
| Adaptive tolerance | runtime tolerance updates and `requested_tolerance` recorded per linear solve |

The local change set since `origin/main` is about 7.6k inserted lines across solver code, profiling code, and campaign scripts.

## What Was Tested

### Direct Baselines

Dogbone direct runs established the timing and nonlinear-iteration baseline:

| solver | steps | nonlinear iterations | linear solves | factorize/setup time | solve time | wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `EigenLDLT` | 100 | 4685 | 4685 | 22.114 s | 11.310 s | 285.078 s |
| `PardisoLDLT` | 100 | 4685 | 4685 | 2.067 s | 7.058 s | 252.136 s |

TS-N_65 Pardiso two-step reference:

| solver | steps | nonlinear iterations | linear solves | setup/factor time | solve time | wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `PardisoLDLT` | 2 | 945 | 945 | 53.468 s | 731.026 s | 4589 s |

### AMGCL on Dogbone

AMGCL moved through several versions:

- Initial scalar/near-nullspace settings often hit iteration caps and were not robust.
- Elastic full-lift AMGCL with rigid-body near-nullspace fixed Dogbone convergence.
- The best Dogbone AMGCL settings solved strictly at tolerances around `5e-6` to `1e-5`.

Dogbone tuning highlights:

| solver | tolerance | status | median Krylov iters | max Krylov iters | linear time | wall time | agreement |
| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |
| `AmgclCGElastic` | `5e-6` | passed | 7 | 7 | 17.599 s | 273.618 s | strict same |
| `AmgclCGElastic` | `1e-5` | passed | 7 | 7 | 17.903 s | 279.141 s | same nonlinear trajectory |
| `AmgclCGElastic` | `1e-4` and looser | completed | 5 to 6 | 6 | lower | about 277 s | changed nonlinear trajectory |

On Dogbone, AMGCL became robust but did not beat Pardiso in wall time because Dogbone is small and much of the run is nonlinear/assembly work.

### hypre on Dogbone

hypre BoomerAMG-CG passed Dogbone gates and strict full-run checks.

| solver | tolerance | steps | nonlinear iterations | median Krylov iters | max Krylov iters | linear time | wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `HypreBoomerAMGCG` | `1e-6` | 100 | 4685 | 15 | 18 | 94.159 s | 355.738 s |
| `HypreBoomerAMGCG` | `1e-5` | 100 | 4685 | 13 | 16 | 85.756 s | 343.610 s |
| `HypreBoomerAMGCG` | `3e-5` | 100 | 4685 | 11 | 15 | 80.895 s | 339.210 s |

On Dogbone, hypre was robust but slower than both Pardiso and AMGCL. That changed on TS-N_65.

### TS-N_65 AMG Gate

The first TS-N_65 gate was first step, five nonlinear solves, strict `1e-6` true residual:

| solver/setup | result | details |
| --- | --- | --- |
| AMGCL block ILU0 | failed | NaN on first solve |
| AMGCL scalar lifted | failed | 500 iterations each solve, true residual around `5e-3` median |
| AMGCL block SPAI0 | failed | 500 iterations each solve, true residual around `3e-1` median |
| hypre BoomerAMG-CG | passed | 56, 117, 92, 68, 68 iterations; max true residual below `1e-6` |

This made hypre the first large-case iterative solver eligible for a longer TS-N_65 comparison.

### AMGCL TS-N_65 Rescue

The large-case AMGCL failure was traced mainly to scaling. The lifted elastic system mixed translational and rotational DOFs in raw stiffness units. Adding symmetric diagonal equilibration plus physical true-residual validation fixed the first TS-N_65 linear solve.

Best first-system comparison:

| solver | success | iterations | true residual | solve time | setup time | total linear time |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| diagonal-scaled AMGCL-CG | yes | 122 | `9.504e-7` | 17.995 s | 5.444 s | 23.439 s |
| hypre BoomerAMG-CG | yes | 60 | `9.877e-7` | 8.731 s | 8.253 s | 16.984 s |

AMGCL is no longer broken on the first large system, but hypre is still faster and more mature for TS-N_65 progression.

### TS-N_65 Two-Step hypre Comparison

The two-step campaign compared Pardiso and hypre at three tolerances:

| solver | tolerance | status | wall time | solve time | median/max iters | relative displacement L2 vs Pardiso at step 2 |
| --- | ---: | --- | ---: | ---: | --- | ---: |
| `PardisoLDLT` | direct | completed | 4589 s | 731 s | direct | reference |
| `HypreBoomerAMGCG` | `1e-1` | completed | 4927 s | 935 s | 7 / 16 | `9.487e-4` |
| `HypreBoomerAMGCG` | `1e-3` | completed | 6842 s | 2921 s | 21 / 36 | `5.625e-6` |
| `HypreBoomerAMGCG` | `1e-6` | completed | 13580 s | 9657 s | 71 / 139 | `1.026e-8` |

Interpretation:

- hypre-CG at `1e-1` is closest to Pardiso wall time and may be the fastest already verified TS-N_65 iterative path.
- `1e-3` is a conservative accuracy/speed compromise.
- `1e-6` is effectively reference-quality but much slower than Pardiso for the first two steps.

The attempted ten-step hypre validation at `1e-1` did not produce a completed ten-step validation. The preserved run had 2140 linear solves, median/max Krylov 7/16, about 10,890 s wall time, and available step-1/step-2 VTK comparison only.

### Dogbone DeflatedFGMRES

The first native `DeflatedFGMRES` implementation was tested with AMGCL and hypre preconditioners.

AMGCL as a raw V-cycle was not a usable FGMRES preconditioner. Nested AMGCL-CG worked but was slower than using AMGCL-CG directly.

hypre as a right-preconditioner was the useful route. The first Dogbone sweep exposed that `dfgmres_deflation_eps=1e-3` was too high and admitted only one vector. After correcting newest-first A-orthogonal rebuilds and using `eps=1e-15`, the basis grew correctly and deflation reduced Krylov work.

Dogbone strict `1e-6` DFGMRES+hypre results:

| run | max basis | outer iterations | median/max iters | linear solve time | wall time |
| --- | ---: | ---: | --- | ---: | ---: |
| N0 | 0 | 70008 | 15 / 18 | 48.814 s | 305.456 s |
| N5 | 5 | 35956 | 7 / 18 | 30.020 s | 294.660 s |
| N10 | 10 | 33774 | 7 / 18 | 28.095 s | 292.276 s |
| N20 | 20 | 31425 | 6 / 17 | 27.692 s | 304.933 s |
| N30 | 30 | 29745 | 6 / 17 | 27.707 s | 317.619 s |

Deflation clearly reduced linear iterations and solve time on Dogbone. The best wall time was N10, while N20/N30 paid more deflation overhead and nonlinear/other variation.

### TS-N_65 DeflatedFGMRES

At strict `1e-6` true residual, `DeflatedFGMRES+hypre` completed two TS-N_65 steps with N5 and N10:

| run | status | steps | nonlinear iterations | linear solves | outer iterations | median/max iters | solve time | wall time |
| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| N0 | partial | 1 | 18 | 382 | 19580 | 51 / 63 | 3113 s | 4701 s |
| N5 | completed | 2 | 945 | 945 | 33669 | 36 / 80 | 5105 s | 8839 s |
| N10 | completed | 2 | 945 | 945 | 34162 | 36 / 134 | 5323 s | 9218 s |

Deflation was active and the runs completed. At strict `1e-6`, DFGMRES N5 was faster than hypre-CG `1e-6`: hypre-CG used 9657 s of solve time and 13580 s wall time, while DFGMRES N5 used 5105 s of solve time and 8839 s wall time. It was still slower than Pardiso and slower than the practical loose hypre-CG runs, and the later loose-tolerance DFGMRES checks showed sensitivity.

At loose `1e-1`, DFGMRES N5 completed two steps much faster:

| run | tolerance | steps | nonlinear iterations | median/max iters | solve time | wall time |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| DFGMRES+hypre N5 | `1e-1` | 2 | 960 | 1 / 18 | 397.584 s | 4410 s |

But the VTK comparison showed this is not accurate enough for TS-N_65:

| metric | value |
| --- | ---: |
| step-2 relative L2 displacement difference vs Pardiso | `6.193e-3` |
| max nodal difference / Pardiso max displacement | `1.09e-1` |
| max pointwise relative difference on nonzero-reference points | `2.50e-1` |

The nonlinear solver accepted the step because global errors were below `1e-3`, but the local displacement/damage path differed. This is the key correctness warning for loose inexact solves on this fracture problem.

### Adaptive Tolerance

Adaptive tolerance was added and recorded through `requested_tolerance`. The first policy tried loose `1e-1` followed by tight `1e-6` when the nonlinear residual approached the trigger.

Result for TS-N_65 DFGMRES+hypre N5:

| item | value |
| --- | ---: |
| completed load steps | 1 |
| linear solves before stop | 60 |
| loose solves | 7 |
| tight solves | 53 |
| outer iterations | 2495 |
| stopped wall time | 854 s |

The run stalled in step 2 after switching to `1e-6`: residual error rose from `7.415510e-03` to `8.008438e-03`, above the nonlinear tolerance `1e-3`. The conclusion is that "loose first, strict later" is not enough for path-dependent fracture if the early loose increments already move the material state onto a different branch.

### Dogbone Accuracy Visualization

A final Dogbone VTK comparison between Pardiso and hypre-CG `1e-1` showed that Dogbone is much less sensitive than TS-N_65:

| metric | value |
| --- | ---: |
| final full-field relative L2 displacement difference | `3.625e-7` |
| max nodal displacement difference | `5.188e-10` |
| crack-opening relative L2 difference | `1.045e-6` |
| max displacement-difference over loading animation | `2.687e-9` |

This supports using Dogbone as a fast mechanical gate, but it also shows Dogbone is not sufficient to validate large-case crack-path accuracy.

## How The Work Moved

| date | movement |
| --- | --- |
| 2026-05-02 | Bootstrapped the deflation workspace and documented the OAS solver structure, benchmark policy, and proposed deflation design. |
| 2026-05-02 | Added Dogbone profiling, direct baselines, and first profiling reports. |
| 2026-05-03 | Added AMGCL and hypre iterative solvers, Dogbone gates, Dogbone full tuning, and TS-N_65 first-step gates. hypre became the first large-case iterative path to pass. |
| 2026-05-03 | Diagnosed and partly rescued AMGCL on TS-N_65 by adding diagonal scaling and true residual validation. |
| 2026-05-04 | Added native `DeflatedFGMRES`, first with AMGCL and then with hypre right preconditioning. Corrected basis rebuilds and showed Dogbone iteration reductions from deflation. |
| 2026-05-04 to 2026-05-05 | Ran TS-N_65 DFGMRES two-step tests, VTK comparisons, and inexact-solve investigations. Found that loose DFGMRES tolerance is not production-safe. |
| 2026-05-05 | Added adaptive tolerance tracking and tested a loose-to-tight DFGMRES policy. It worked mechanically but did not recover the safe TS-N_65 path. |

## Current Position

The solver stack is now instrumented enough to make decisions from data. The current best practical recommendations are:

| goal | current best path |
| --- | --- |
| Fast small-case development | Dogbone with profiling enabled |
| Dogbone robust iterative solver | `AmgclCGElastic` at about `5e-6` to `1e-5` |
| TS-N_65 fastest verified iterative path | `HypreBoomerAMGCG`, `hypre_tolerance=1e-1`, but validate beyond two steps before production use |
| TS-N_65 conservative iterative path | `HypreBoomerAMGCG`, `hypre_tolerance=1e-3` |
| TS-N_65 reference-quality iterative path | `DeflatedFGMRES+hypre` N5/N10 at `1e-6` is the faster strict two-step result; `HypreBoomerAMGCG` at `1e-6` remains the simpler strict baseline |
| Deflation research path | Continue `DeflatedFGMRES+hypre` with strict tolerance and stronger acceptance checks |
| Avoid for production | fixed or early loose `DeflatedFGMRES+hypre` at `1e-1` on TS-N_65 |

## Open Risks and Next Work

- The TS-N_65 two-step hypre evidence is strong enough for direction, but not enough for production claims. A stable longer-run validation is still needed.
- Loose linear tolerances need a problem-aware forcing policy. The first adaptive DFGMRES policy did not protect the material path.
- Accuracy checks should include crack/damage/internal variables and load-displacement gauges, not only displacement VTU fields.
- AMGCL is worth keeping, but hypre is ahead for large-case performance and robustness.
- Deflation should be evaluated with stricter acceptance rules, possibly a final tight correction before accepting a step, because global nonlinear convergence alone did not guarantee local state agreement.

## Key Evidence Artifacts

- `results/amgcl-dogbone-tsn65-20260503-002607/summary.md`
- `results/amgcl-hypre-dogbone-full-20260503-055920/summary.md`
- `results/dogbone-amg-tuning-20260503-072917/summary.md`
- `results/tsn65-amg-firststep-gate-20260503-062440/summary.md`
- `results/amgcl-tsn65-first-solve-rescue-20260503-214332/summary.md`
- `results/tsn65-first-system-amgcl-hypre-sweep-20260503-224637/summary.md`
- `results/tsn65-two-step-comparison-20260503-093031/summary.md`
- `results/dogbone-dfgmres-hypre-deflation-fixed-20260504-135000/summary.md`
- `results/tsn65-dfgmres-hypre-deflation-20260504-144750/summary.md`
- `results/tsn65-step2-pardiso-vs-dfgmres-tol1em1-vtk-20260505-132547/comparison.md`
- `results/tsn65-step2-pardiso-vs-dfgmres-tol1em1-vtk-20260505-132547/inexact-linear-solve-investigation.md`
- `results/tsn65-dfgmres-hypre-adaptive-tolerance-20260505-165534/adaptive-tolerance-findings.md`
- `results/tsn65-hypre-10step-validation-20260505-181811/summary.md`
- `results/dogbone-pardiso-vs-hypre-tol1em1-vtk-20260505-153145/comparison.md`
- `results/dogbone-pardiso-vs-hypre-tol1em1-animation-20260505-160920/animation.md`

Note: `scripts/run-tsn65-hypre-10step-campaign.py` was present as an untracked local script while this report was written.
