# OpenMP Bottleneck Review

Date: 2026-05-06

Scope: OAS-side serial bottlenecks relevant to the deflation and iterative-solver work, especially TS-N_65. This is a code review, not a measured OpenMP implementation report.

## Executive Finding

The current executable is linked with OpenMP and `main.cpp` sets the requested thread count, but the OAS assembly and nonlinear update code has effectively no OpenMP work-sharing loops. Threading today mostly comes from external solvers such as MKL/Pardiso and hypre.

After the linear solver was improved, the expected serial bottleneck is no longer only the linear solve. The large TS-N_65 profiles show a large serial floor in force evaluation:

| profile/run | wall time | linear total | `forces.total active` | `forces.integrate_internal_forces active` | `forces.simplex_volumetric_strains active` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `dfgmres-hypre-tol1em2-N10`, 2026-05-06 sweep | 4187 s | 659 s | 3385 s | 3298 s | 82.3 s |
| `dfgmres-hypre-tol1em1-N20`, 2026-05-06 sweep | 3805 s | 279 s | 3353 s | 3267 s | 80.8 s |
| `dfgmres-hypre-N5 tol1e-1`, 2026-05-05 profile | 4410 s | 409 s | 3860 s | 3753 s | 101.5 s |

For the best loose DFGMRES runs, about 85-90% of wall time is now in OAS force evaluation. That is the main OpenMP opportunity.

## Current OpenMP Surface

Relevant observations:

- `CMakeLists.txt` calls `find_package(OpenMP REQUIRED)` and links `OpenMP::OpenMP_CXX`.
- `src/solver/src/main.cpp:115-125` disables dynamic OpenMP teams and sets `omp_set_num_threads()` from `-j` or `OMP_NUM_THREADS`.
- `src/solver/src` has no meaningful `#pragma omp parallel for` loops in the solver, element, node, or material containers.
- MKL/Pardiso explicitly consumes `MKL_NUM_THREADS`/`OMP_NUM_THREADS` in `src/solver/src/linalg.cpp`.
- hypre is built/configured with OpenMP support and the matrix path uses `HYPRE_IJMatrixSetOMPFlag`, but this only helps hypre internals, not OAS assembly.

## Bottleneck Assessment

| priority | bottleneck | evidence | OpenMP payoff | implementation difficulty | recommendation |
| ---: | --- | --- | --- | --- | --- |
| 1 | Element stress/force integration | `ElementContainer::integrateInternalForces()` dominates TS-N_65 runtime | Very high | Medium-hard | First target |
| 2 | DFGMRES sparse matrix-vector products | DFGMRES profiles show 59-403 s in matvec depending on tolerance/basis | Medium | Medium | Useful after force loop, or for strict solves |
| 3 | Simplex volumetric strain update | 80-110 s over two TS-N_65 steps | Medium-low | Medium | Good second OAS target |
| 4 | Stiffness matrix assembly | about 27 s over 3 updates in two-step TS-N_65 runs | Low-medium | Hard | Do after force loop |
| 5 | Step-end material status update | about 0.5 s per accepted TS-N_65 step | Low | Easy-medium | Cheap cleanup target |
| 6 | Dense vector scatter/reduction kernels | 1-5 s totals over two steps for reductions/reactions | Low | Easy-medium | Opportunistic |
| 7 | Constraint matrix transform | about 3.6-4.0 s over 3 transforms | Low | Hard | Avoid manual OpenMP initially |
| 8 | DFGMRES deflation/orthogonalization vector ops | up to tens of seconds in strict runs; smaller in loose runs | Low-medium | Medium | Only after higher-payoff work |

## 1. Element Stress/Force Integration

Code:

- `src/solver/src/solver.cpp:322-347`: `Solver::computeInternalExternalForces()`
- `src/solver/src/element_container.cpp:658-690`: `ElementContainer::integrateInternalForces()`
- `src/solver/src/element.cpp:272-323`: per-element strain, stress, and internal-force work

Current behavior:

1. Zero the global internal force vector.
2. Reset eigenstrain and run material preparation.
3. Loop by `solution_order`.
4. For every element in that order, gather element DOFs, evaluate strains, and evaluate stresses.
5. Loop over all elements again, compute element internal force, and scatter into the shared global `full_f`.

Expected OpenMP approach:

- Keep the outer `solution_order` loop serial, because `DiscreteTrsprtCoupledElem` uses `solution_order = 1` and can depend on already evaluated LDPM tetras.
- Parallelize the element loop inside each `solution_order`.
- Split force assembly into thread-local accumulation buffers, then reduce them into `full_f`.

Sketch:

```cpp
std::vector<Vector> thread_forces(num_threads, Vector::Zero(full_f.size()));

for (unsigned so = 0; so <= max_sol_order; ++so) {
    #pragma omp parallel for schedule(dynamic)
    for (long idx = 0; idx < elems.size(); ++idx) {
        Element *e = elems[idx];
        if (e->giveSolutionOrder() != so) continue;
        // gather element dofs, evaluate strains, evaluate stresses
    }
}

#pragma omp parallel for schedule(dynamic)
for (long idx = 0; idx < elems.size(); ++idx) {
    int tid = omp_get_thread_num();
    // compute element internal forces and add to thread_forces[tid]
}

// reduce thread_forces into full_f
```

Main hazards:

- The global scatter `full_f[elDoFs[i]] += ...` is a race. Use thread-local dense vectors, graph coloring, or atomics. Thread-local dense vectors are the simplest first implementation. For roughly 600k DOFs and 16 threads, one `double` vector per thread is about 77 MiB, which is acceptable compared with direct-solver memory.
- `materials->runPreparationForStressEvaluation(this)` is a global material hook. It should stay outside the parallel region until material-specific implementations are audited.
- `resetEigenStrain(time)` currently loops over all elements and then applies eigenstrain loads. Keep it serial initially.
- Element stress updates mutate temporary material status owned by the element. That is usually thread-safe per element, but coupled elements that read neighbor/tetra state need the `solution_order` barrier.
- Some material implementations may use shared material objects for constants and per-element status objects for state. This is a good pattern for OpenMP, but it must be checked for hidden static/shared scratch buffers.

Difficulty: medium-hard.

Payoff: very high. If only active force integration in `dfgmres-hypre-tol1em2-N10` got an ideal 8x speedup, wall time would drop from about 4187 s to roughly 1300 s before secondary effects. Real speedup will be lower because of memory bandwidth, force-vector reduction, and element imbalance, but this is still the main serial floor.

## 2. DFGMRES Sparse Matrix-Vector Products

Code:

- `src/solver/src/linalg.cpp:1470-1472`
- `src/solver/src/linalg.cpp:1504-1506`
- `src/solver/src/linalg.cpp:1562-1565`

Current behavior:

- `DeflatedFGMRES` uses `impl->activeMatrix * vector`, i.e. Eigen sparse matrix-vector multiplication.
- The class already builds row-major CRS arrays in `makeRowMajorCrs()` for AMGCL/hypre setup, but those arrays are not used for DFGMRES matvecs.

Expected OpenMP approach:

- Add a local CRS SpMV routine over rows:

```cpp
#pragma omp parallel for schedule(static)
for (ptrdiff_t row = 0; row < rows; ++row) {
    double sum = 0.0;
    for (ptrdiff_t p = ptr[row]; p < ptr[row + 1]; ++p) {
        sum += val[p] * x[col[p]];
    }
    y[row] = sum;
}
```

Main hazards:

- The matrix may be lifted/scaled, so the CRS arrays must correspond to `activeMatrix`, not `reducedMatrix`.
- Need strict residual reproducibility checks because parallel summation changes roundoff.
- Memory bandwidth will cap scaling, but row-parallel SpMV should still be better than serial Eigen sparse multiply if Eigen is not parallelizing this path.

Difficulty: medium.

Payoff: medium. In the 2026-05-06 sweep, DFGMRES matvec time ranged from about 59 s to 403 s over two steps. This matters most for strict or low-basis runs. It will not address the 3.3-4.3 ks force-evaluation floor.

## 3. Simplex Volumetric Strain Update

Code:

- `src/solver/src/solver.cpp:327`
- `src/solver/src/node_container.cpp:151-164`
- `src/solver/src/simplex.cpp:104-140`

Current behavior:

- First pass computes valid simplex volumetric strain.
- Second pass lets invalid simplices steal/average volumetric strain from neighbors.

Expected OpenMP approach:

- Parallelize the first pass if each node owns a unique simplex pointer, or better, build a unique simplex list and parallelize over that.
- Add a barrier between the valid-simplex update and neighbor-stealing pass.
- Parallelize the stealing pass only after verifying each iteration writes to a unique simplex.

Main hazards:

- The current loop iterates over nodes, not explicitly over unique simplices. If multiple nodes can reference the same `Simplex`, naive OpenMP would race on `volstrain` and `updated`.
- Invalid-simplex stealing reads neighbors and writes the current simplex. That is safe only if no two iterations write the same simplex.

Difficulty: medium.

Payoff: medium-low. This costs roughly 0.10 s per active force evaluation on TS-N_65, or 80-110 s over two-step runs. It is worth doing after the main element integration loop, but it is not the dominant cost.

## 4. Stiffness Matrix Assembly

Code:

- `src/solver/src/solver_implicit.cpp:791-806`
- `src/solver/src/element_container.cpp:505-598`
- `src/solver/src/element.cpp:254-269`

Current behavior:

- `K = K * 0` clears the sparse matrix.
- The code loops over elements, builds the local stiffness matrix, maps element DOFs to solver DOFs, and updates `K.coeffRef(...)`.
- This is serial and `coeffRef` writes are not thread-safe.

Expected OpenMP approach options:

1. Precompute an element-to-global-nonzero map after `prepareStructuralMatrix()` and update `K.valuePtr()` directly.
2. Use graph coloring so concurrently processed elements do not write the same global entries.
3. Use thread-local triplet/value buffers and combine after the parallel region.
4. Use atomics on `K.valuePtr()` entries only if a stable entry index map exists.

Main hazards:

- Sparse matrix random writes through `coeffRef()` are unsuitable for naive OpenMP.
- Duplicate element contributions to the same global entry must be summed deterministically enough for nonlinear reproducibility.
- The matrix pattern is stable in the TS-N_65 profiles, so a one-time element-entry map is probably the right design.

Difficulty: hard.

Payoff: low-medium for current TS-N_65 iterative runs. Matrix update is about 27 s over three updates, so it is far below force integration. It becomes more important only if matrix updates are much more frequent or force integration is already parallelized.

## 5. Step-End Material Status Update

Code:

- `src/solver/src/solver.cpp:200-203`
- `src/solver/src/element_container.cpp:415-425`
- `src/solver/src/material_ldpm.cpp:495-502`

Current behavior:

- After an accepted step, OAS loops over all elements and commits temporary material state into persistent state.

Expected OpenMP approach:

- `#pragma omp parallel for schedule(static/dynamic)` over elements.

Main hazards:

- Most element material statuses appear element-owned and independent.
- Audit any material that updates shared/global material data during `update()`.
- Keep step-level solver state copies (`r`, `f_int_old`, `f_ext_old`) serial.

Difficulty: easy-medium.

Payoff: low. About 0.5 s per accepted TS-N_65 step in current profiles. It is worth doing because the change is small, but it will not move two-step wall time much.

## 6. Dense Vector Kernels and Reductions

Code:

- `src/solver/src/solver.cpp:274-290`: full DOF expansion and increment application
- `src/solver/src/node_container.cpp:279-312`: full/reduced array conversion
- `src/solver/src/node_container.cpp:315-341`: external-force/reaction update
- `src/solver/src/solver_implicit.cpp:1221-1277`: error accumulation
- `src/solver/src/solver.cpp:358-383`: energy accumulation

Current behavior:

- Several full-DOF loops are serial.
- The expensive parts are small compared with force integration. For example, in `dfgmres-hypre-tol1em2-N10`, reduced-force extraction is 1.52 s total, field expansion is 1.66 s, error accumulation is 1.26 s, and external reactions are 4.0 s.

Expected OpenMP approach:

- Use `parallel for` for pure copy/scatter loops with independent writes.
- Use OpenMP reductions for physical-field norm accumulation. Since `numPhysicalFields` is tiny, either explicit scalar reductions or per-thread small arrays are cleaner than atomics.
- Keep constraint master-force propagation serial until constraints are audited, because it mutates shared master entries.

Difficulty: easy-medium.

Payoff: low. These are good cleanup changes after the dominant element work.

## 7. Constraint Matrix Transform

Code:

- `src/solver/src/solver_implicit.cpp:489-494`
- `src/solver/src/constraint.cpp:594-603`

Current behavior:

- Uses Eigen sparse products: `Knew = X.transpose() * K * X`.
- Costs roughly 1.2-1.4 s per transform on TS-N_65.

Expected OpenMP approach:

- Do not start with manual OpenMP here.
- First check whether Eigen sparse products are already using OpenMP in this build.
- Better optimization may be avoiding repeated transforms or caching the reduced pattern, not hand-parallelizing sparse triple products.

Difficulty: hard.

Payoff: low for current profiles.

## 8. DFGMRES Deflation and Orthogonalization Operations

Code:

- `src/solver/src/linalg.cpp:1080-1089`: deflation projection
- `src/solver/src/linalg.cpp:1461-1468`: initial coarse correction
- `src/solver/src/linalg.cpp:1508-1520`: modified Gram-Schmidt

Current behavior:

- Uses serial-looking Eigen dot/axpy operations over large vectors.
- In the latest TS-N_65 sweep, orthogonalization and deflation were normally smaller than preconditioner apply and matvec, but can still reach tens of seconds.

Expected OpenMP approach:

- Parallelize dot products with reductions and axpy updates over vector rows, or rely on Eigen if configured to parallelize dense operations.
- For each Krylov vector, Gram-Schmidt itself remains sequential over previous basis vectors, but the vector dot/axpy operations inside each step are parallelizable.
- Deflation projection over a small basis can parallelize vector operations, but basis loop order remains serial unless rewritten as blocked dense GEMV.

Difficulty: medium.

Payoff: low-medium. Do this only after force integration and SpMV, unless strict DFGMRES runs become the primary production path.

## Suggested Implementation Order

1. Add a guarded OpenMP path for `ElementContainer::integrateInternalForces()`.
   - Parallel stress/strain evaluation by `solution_order`.
   - Use thread-local dense `full_f` buffers for force scatter.
   - Validate Dogbone against Pardiso/hypre final gauges and VTK.

2. Add OpenMP to `NodeContainer::updateSimplexVolumetricStrains()` using a unique simplex list.
   - Validate that duplicate simplex pointers cannot race.

3. Add a row-parallel CRS SpMV for `DeflatedFGMRES`.
   - Use existing CRS arrays from factorization.
   - Compare residual histories and VTK outputs against the current DFGMRES path.

4. Parallelize step-end material status updates.
   - Audit material `update()` functions for shared writes.

5. Parallelize low-risk dense vector loops and reductions.
   - Keep constraint propagation serial until a dedicated constraint audit is done.

6. Revisit stiffness matrix assembly only after the force path is parallel.
   - Implement a stable element-to-entry map instead of parallelizing `coeffRef()`.

## Practical Expectations

The first useful OpenMP milestone is not "make everything parallel"; it is "make active force integration scale without changing the nonlinear path." That one change attacks the current serial floor directly.

Expected difficulty summary:

| task | difficulty | reason |
| --- | --- | --- |
| Parallel element stress evaluation | Medium | per-element state is mostly independent, but solution-order barriers and material hooks matter |
| Parallel force scatter | Medium-hard | global force vector write races require thread-local buffers, coloring, or atomics |
| Parallel simplex update | Medium | ownership/duplicate-simplex semantics must be verified |
| Parallel DFGMRES SpMV | Medium | straightforward CSR row parallelism, but numerical reproducibility must be checked |
| Parallel material status commit | Easy-medium | likely independent per element, low payoff |
| Parallel stiffness assembly | Hard | sparse global writes and deterministic accumulation are nontrivial |
| Parallel constraint transform | Hard | sparse triple product; low payoff |

The dominant risk is correctness, not syntax. TS-N_65 is path-dependent: small changes in arithmetic order can change accepted local crack/damage states even if global nonlinear errors pass. Every OpenMP change should be gated by:

- exact or near-exact Dogbone gauge comparison,
- TS-N_65 first-system/first-step linear and nonlinear checks,
- TS-N_65 VTK displacement and, ideally, crack/damage/internal-variable comparisons,
- fixed thread-count reproducibility checks across repeated runs.
