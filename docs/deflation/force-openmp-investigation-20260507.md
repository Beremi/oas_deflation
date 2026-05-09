# Force Integration OpenMP Investigation

Date: 2026-05-07

Scope: `forces.total active` and `forces.integrate_internal_forces active` for the TS-N_65 DFGMRES/Hypre deflation path after the OpenMP corrections. Old long runs were not rerun.

## Baseline Evidence

Two-step corrected run:

- Current artifact: `results/tsn65-openmp-corrected-two-step-20260507/dfgmres-hypre-tol1em1-N20`
- Previous artifact: `results/tsn65-dfgmres-hypre-tol-basis-sweep-20260506-053205/dfgmres-hypre-tol1em1-N20`
- Config: `tol=1e-1`, `N=20`, `total_time=0.01`, 16 OAS threads.

| phase | previous | current | speedup |
| --- | ---: | ---: | ---: |
| wall | 3805.4 s | 1222.2 s | 3.11x |
| `forces.total active` | 3352.5 s | 804.5 s | 4.17x |
| `forces.integrate_internal_forces active` | 3267.5 s | 701.3 s | 4.66x |
| `forces.simplex_volumetric_strains active` | 80.8 s | 98.2 s | 0.82x |

The current two-step run completed, but took more nonlinear iterations than the older run:

| metric | previous | current |
| --- | ---: | ---: |
| nonlinear iterations | 807 | 964 |
| DFGMRES outer iterations | 1278 | 438 |
| max linear iterations/solve | 23 | 4 |
| final nonlinear residual | 9.982675e-4 | 9.928774e-4 |

So the force loop is substantially faster per evaluation, while the nonlinear trajectory did more evaluations in this specific run.

## Short Scaling Test

To avoid repeatedly running the full two-step case, I used a first-step TS-N_65 profile:

- `results/tsn65-force-scaling-firststep-current-20260507-j1/dfgmres-hypre-tol1em1-N20`
- `results/tsn65-force-scaling-firststep-current-20260507-j16/dfgmres-hypre-tol1em1-N20`

Command shape:

```bash
python3 scripts/run-tsn65-dfgmres-hypre-deflation.py \
  --oas-bin ../oas_deflation-build/release/bin/OAS \
  --threads <1-or-16> \
  --total-time 5.000000e-03 \
  --linear-tol 1e-1 \
  --true-tol 1e-1 \
  --nvecs 20 \
  --only dfgmres-hypre-tol1em1-N20 \
  --out-dir <result-dir>
```

Both runs followed the same numerical path:

| metric | `j=1` | `j=16` |
| --- | ---: | ---: |
| accepted steps | 1 | 1 |
| nonlinear iterations | 19 | 19 |
| linear solves | 19 | 19 |
| DFGMRES outer iterations | 37 | 37 |
| max linear iterations/solve | 4 | 4 |
| final true relres | 6.2263925e-2 | 6.2263925e-2 |
| final nonlinear residual | 7.875648e-4 | 7.875648e-4 |

Scaling:

| phase | `j=1` | `j=16` | speedup |
| --- | ---: | ---: | ---: |
| step elapsed | 201.4 s | 118.5 s | 1.70x |
| wall by runner | 205.0 s | 122.3 s | 1.68x |
| `forces.total active` | 82.82 s | 17.40 s | 4.76x |
| `forces.integrate_internal_forces active` | 80.51 s | 15.15 s | 5.31x |
| `forces.simplex_volumetric_strains active` | 2.20 s | 2.14 s | 1.02x |
| linear solve time | 15.84 s | 15.29 s | 1.04x |
| stiffness update | 18.57 s | 9.82 s | 1.89x |

Interpretation:

- The force integration kernel scales, but only to about 33% parallel efficiency at 16 threads (`5.31 / 16`).
- Whole-step speedup is much lower because Hypre apply is intentionally serialized for iteration stability, and setup/factorization/export-like parts remain serial or weakly threaded.
- Simplex volumetric strain is not the cause of poor `integrate_internal_forces` scaling in this sample.

TS-N_65 element makeup in this case:

| item | value |
| --- | ---: |
| elements | 514253 |
| element type | `LDPMTetra` only |
| nodes file rows | 92416 |
| free DoF reported by OAS | 567923 |

## Tested Small Code Hypotheses

I tested and rejected two small changes.

1. Use a const-reference DOF accessor and reuse one gather vector per OpenMP thread.
   - Result: `forces.integrate_internal_forces active` changed from 15.15 s to 15.69 s on the first-step `j=16` profile.
   - Decision: reverted. It removes some allocation/copy churn in source code, but did not improve this real kernel.

2. Change the strain/stress OpenMP loop from `schedule(dynamic)` to `schedule(static)`.
   - Result: `forces.integrate_internal_forces active` changed from 15.15 s to 16.24 s.
   - Decision: reverted. Dynamic scheduling appears to help, likely because material branches and damage states make element costs nonuniform even though all elements are `LDPMTetra`.

## Current Hot-Path Shape

`ElementContainer::integrateInternalForces()` currently does:

1. Serial `full_f.setZero()`.
2. Serial `resetEigenStrain(time)`.
3. Serial `materials->runPreparationForStressEvaluation(this)`.
4. For each `solution_order`, parallel loop over all elements:
   - copy element DOF IDs,
   - gather element DOF values,
   - `Element::evaluateStrains()`,
   - `Element::evaluateStresses()`.
5. Allocate and zero one dense global force vector per OpenMP thread.
6. Parallel loop over elements:
   - compute local internal force vector,
   - scatter into the thread-local dense global vector.
7. Serial deterministic reduction of thread-local force vectors into `full_f`.

For TS-N_65, `solution_order` is effectively not the bottleneck because all elements in `elems.inp` are `LDPMTetra`.

## Likely Bottlenecks

The measurements point away from simple OpenMP scheduling overhead and toward the structure of the kernel:

- Each force evaluation touches 514k elements and 12 integration points per LDPM tetra.
- Element data is object-heavy: element pointers, material status pointers, dynamic Eigen `Vector`/`Matrix` operations, virtual calls, and scattered DOF reads/writes.
- Thread-local dense force vectors avoid races but multiply the global force-vector footprint by thread count. At 16 threads and ~568k free rows this is not huge in memory capacity, but every force call still zeroes, writes, reads, and reduces many large vectors.
- The scatter pattern is poor for cache locality. Each thread writes to a large dense vector at element DOF locations; adjacent elements do not necessarily touch adjacent DOFs.
- The active stress path has branchy material logic, so static scheduling under-loads some threads.

## External References

The current design is close to the standard shared-memory FEM pattern, but the references explain why scaling can stop early:

- The WorkStream FEM pattern separates parallel local contribution computation from the global reduction/scatter, and explicitly calls out the need to avoid races and keep reduction order repeatable. It also recommends moving expensive scratch allocation/precomputation out of the per-cell operation when possible. Source: Turcksin et al., "WorkStream -- A Design Pattern for Multicore-Enabled Finite Element Computations" (`https://www.math.colostate.edu/~bangerth/publications/2013-pattern.pdf`).
- deal.II 9.7 reports that sequential copy operations can cap assembly scalability to fewer than 10 threads, while cell coloring removes that bottleneck in their step-32 assembly benchmark. Source: deal.II 9.7 report (`https://www.osti.gov/servlets/purl/3006478`).
- Bošanský and Patzák evaluate OpenMP FEM assembly strategies in OOFEM, including synchronization, block locking, and lock-free coloring, and emphasize benchmarking because the best strategy is workload-dependent. Source: Acta Polytechnica paper (`https://doi.org/10.14311/AP.2020.60.0025`).
- MFEM performance notes emphasize that finite element kernels often have low arithmetic intensity and can be limited by memory traffic/access patterns. Source: MFEM GPU support notes (`https://mfem.org/gpu-support/`).

## Recommendations

### 1. Add subphase profiling before further solver edits

Before changing architecture, split `forces.integrate_internal_forces` into internal subphases:

- eigenstrain/material preparation,
- strain/stress evaluation,
- local force computation,
- thread-local scatter,
- thread-local force-vector zero/reduction.

Right now the profiler only sees the combined `integrate_internal_forces` block from `Solver::computeInternalExternalForces()`. Without subphase timing, we cannot tell whether the 5.3x limit comes from material stress evaluation, local force computation, scatter, or reduction.

Implementation difficulty: low-medium. The cleanest version is to let `ElementContainer` optionally receive a profiling callback or a lightweight profiler pointer from `Solver`.

### 2. Prototype blocked WorkStream-style force assembly

Instead of one dense global vector per thread:

1. Process elements in blocks, e.g. 2048 or 8192 elements.
2. Compute local element force vectors in parallel into block-local storage.
3. Scatter those local vectors into `full_f` serially in element order.

Pros:

- No per-thread dense global force vectors.
- Deterministic serial scatter order.
- Better bounded memory footprint.
- Matches the WorkStream pattern and is not too invasive.

Cons:

- If scatter is a small fraction of time, this may not improve wall time.
- If local force vector allocation remains expensive, use fixed-size or reusable per-thread local storage next.

Implementation difficulty: medium. This is the next practical experiment because it is reversible and does not require graph coloring.

### 3. Prototype colored direct scatter if blocked assembly shows scatter/reduction is costly

Build element colors where no two elements in a color share a global DOF, then run:

```cpp
for color in colors:
    #pragma omp parallel for
    for element in color:
        compute local force
        scatter directly to full_f
```

Pros:

- Avoids atomics.
- Avoids thread-local full vectors.
- Literature suggests coloring is the best path when the copy/scatter stage limits scalability.

Cons:

- More code and memory for graph coloring.
- Summation order changes by color, so compare with tolerances.
- Needs careful invalidation if elements or DOF numbering change.

Implementation difficulty: medium-hard.

### 4. Do not add atomics as the first force scatter alternative

Atomic scatter is easy to code, but it would create heavy contention on shared nodal DOFs and change summation order. It is useful only as a quick diagnostic, not as the preferred final implementation.

### 5. LDPM-specific fused kernels are probably the high-payoff architecture change

For TS-N_65, all elements are `LDPMTetra`, and the generic path pays for:

- virtual `Element`/`MaterialStatus` calls,
- dynamic Eigen vector/matrix temporaries,
- per-element local force vector return values,
- pointer-heavy data access.

A specialized LDPM force kernel could gather DOFs, evaluate strain/stress, compute `B^T stress`, and scatter local force in one pass with fixed-size stack arrays. That is likely the biggest single-node performance opportunity, but it is also a larger architectural step and should follow subphase profiling.

Implementation difficulty: hard.

## Kept Tooling

Added:

```bash
python3 scripts/compare-runtime-profiles.py <baseline-run-dir> <candidate-run-dir>
```

Example:

```bash
python3 scripts/compare-runtime-profiles.py \
  results/tsn65-force-scaling-firststep-current-20260507-j1/dfgmres-hypre-tol1em1-N20 \
  results/tsn65-force-scaling-firststep-current-20260507-j16/dfgmres-hypre-tol1em1-N20 \
  --baseline-label j1 \
  --candidate-label j16
```

This produces Markdown tables for:

- run-level metrics,
- runtime phase totals and speedups,
- accepted-step elapsed times,
- linear solve counts/iterations/times.

This is the short-loop comparison harness to use before rerunning the full two-step or campaign profiles.

## Follow-Up: `forces.total` Active Review And Patch

Follow-up run scope: current `DeflatedFGMRES` + hypre configuration with `hypre_relax_type=16`, `hypre_cheby_order=3`, `hypre_num_sweeps=3`, `dfgmres_elastic_reorder=2`, `N=20`, `tol=1e-1`, `j=16`.

Current two-step artifact:

- `results/tsn65-dfgmres-hypre-cheby3-sweeps3-twostep-20260507-j16/dfgmres-hypre-tol1em1-N20`

The current two-step `forces.total active` profile is:

| phase | count | total s | mean s | share of wall |
| --- | ---: | ---: | ---: | ---: |
| `forces.total active` | 964 | 808.570 | 0.8388 | 67.4% |
| `forces.integrate_internal_forces active` | 964 | 704.731 | 0.7310 | 58.7% |
| `forces.integrate.strain_stress active` | 500,960,916 | 384.920 | 7.684e-7 | 32.1% |
| `forces.integrate.element_force_scatter active` | 500,960,916 | 227.513 | 4.542e-7 | 19.0% |
| `forces.simplex_volumetric_strains active` | 964 | 98.916 | 0.1026 | 8.2% |
| `forces.integrate.reset_eigenstrain active` | 964 | 86.528 | 0.0898 | 7.2% |
| `forces.integrate.thread_force_reduce active` | 15,424 | 2.229 | 1.445e-4 | 0.2% |

Interpretation:

- The main force bottleneck is real element work: strain/stress evaluation plus local internal-force/scatter.
- Thread-local dense-vector reduction is not currently a meaningful bottleneck.
- `reset_eigenstrain` was an avoidable cost for TS-N_65: the case input has no eigenstrain loads, but the code still traversed every element and material status every force evaluation.

### Implemented Patch

Added `BCContainer::hasEigenStrainLoads()` and changed `ElementContainer::resetEigenStrain()` to:

1. return immediately when there are no eigenstrain loads;
2. otherwise clear element eigenstrains in parallel;
3. keep `applyEigenStrainLoads(time)` serial after the clear phase.

This keeps the existing semantics for eigenstrain cases while avoiding a full element/status traversal in cases like TS-N_65.

### First-Step Test

Baseline artifact:

- `results/tsn65-dfgmres-hypre-cheby3-sweeps3-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

Candidate artifact:

- `results/tsn65-force-eigenstrain-skip-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

| metric | baseline | candidate |
| --- | ---: | ---: |
| accepted steps | 1 | 1 |
| nonlinear iterations | 19 | 19 |
| linear solves | 19 | 19 |
| DFGMRES outer iterations | 37 | 37 |
| max linear iterations | 5 | 5 |
| final true relres | 0.0688078061 | 0.0688078061 |
| final nonlinear residual | 8.13799e-4 | 8.13799e-4 |

Timing:

| phase | baseline s | candidate s | speedup |
| --- | ---: | ---: | ---: |
| `forces.total active` | 17.395 | 15.888 | 1.095x |
| `forces.integrate_internal_forces active` | 15.148 | 13.656 | 1.109x |
| `forces.integrate.reset_eigenstrain active` | 1.839 | 0.000006 | effectively removed |
| `forces.integrate.strain_stress active` | 8.341 | 8.352 | 0.999x |
| `forces.integrate.element_force_scatter active` | 4.843 | 5.176 | 0.936x |
| `forces.simplex_volumetric_strains active` | 2.143 | 2.125 | 1.009x |
| `linear.solve_total` | 17.179 | 17.290 | 0.994x |

The force kernel improvement is clear, but whole first-step wall time did not improve because init, matrix, and linear timings varied by more than the saved force time in this short run. On the two-step run, the same removed `reset_eigenstrain` work corresponds to about 86.5 seconds of previously measured active-force time.

### Tested And Rejected

I also tested changing the strain/stress loop from `schedule(dynamic)` to `schedule(dynamic, 64)`.

Artifact:

- `results/tsn65-force-eigenstrain-skip-dyn64-firststep-20260507-j16/dfgmres-hypre-tol1em1-N20`

Result:

| metric | skip only | skip + dynamic chunk 64 |
| --- | ---: | ---: |
| `forces.total active` | 15.888 s | 15.983 s |
| `forces.integrate_internal_forces active` | 13.656 s | 13.652 s |
| `forces.integrate.strain_stress active` | 8.352 s | 8.384 s |
| `forces.simplex_volumetric_strains active` | 2.125 s | 2.220 s |
| nonlinear / linear iteration counts | unchanged | unchanged |

Decision: reverted. Chunking was neutral for `integrate_internal_forces` and slightly worse for total force time in this sample.

### Remaining Improvement Suggestions

1. **Avoid repeated eigenstrain traversal where possible.** Done for no-eigenstrain cases. This is low risk and directly removes about 7% of current two-step wall time from the old profile.

2. **Do not focus on thread-local force-vector reduction first.** It costs only 2.23 s out of 1200 s in the current two-step run.

3. **Next practical target: reduce `element_force_scatter`.** The current generic path calls `giveInternalForces()` and `giveDoFs()` in a second full pass over all elements. A prototype should evaluate whether a fused stress/internal-force pass is safe for all `solution_order` cases. This could remove one element traversal, one `giveDoFs()` copy per element, and some temporary `Vector` traffic.

4. **Highest-payoff architectural target: LDPM-specific fixed-size kernels.** TS-N_65 is dominated by `LDPMTetra`. The generic virtual/dynamic-Eigen path does many tiny matrix-vector operations and temporary allocations. A specialized LDPM mechanical force kernel using fixed-size stack storage would attack both dominant rows: `strain_stress` and `element_force_scatter`.

5. **Coloring/direct scatter remains a later experiment.** Current thread-local reduction is cheap, so coloring is only worth it if a fused/local-force prototype shows the scatter pass itself, not local force computation, is the limiter.

### Verification

After keeping only the eigenstrain fast path:

```bash
PKG_CONFIG_PATH=/home/beremi/local/librsb-1.3.0.2-native-openmp/lib/pkgconfig \
cmake --build ../oas_deflation-build/release --target OAS -j16

OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE \
ctest --test-dir ../oas_deflation-build/release -L tests --output-on-failure
```

CTest result: 7/7 passed.
