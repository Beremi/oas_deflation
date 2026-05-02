# OAS Linear Solver Deflation Report

This is the compact project report for the OAS deflation workspace. The full source note is `docs/deflation/oas_linear_solver_review_parallel_amended.tex`.

## Current OAS Solver Structure

OAS selects the linear solver through `solver_type` in `solver.inp`. The main implementation points are:

- `src/solver/src/linalg.h` and `src/solver/src/linalg.cpp`: linear solver abstractions and backend wrappers.
- `src/solver/src/solver_implicit.cpp`: implicit linear and nonlinear solve loops.
- `SteadyStateNonLinearSolver::solve()`: the direct-control branch calls `linalgsolver->solve(ddr, f)`.
- `computeKeff()` and `factorizeLinearSystem()`: build the effective matrix, run one-time symbolic analysis, and refactor when the matrix is updated.

The existing `EigenConj` solver is conjugate gradients with Eigen's diagonal/Jacobi preconditioner. It is useful as an algorithmic baseline, not as a credible large fracture-case solver.

## Baselines

The Dogbone case has roughly 6.6k DOF and is the debugging target. It should be used for build validation, instrumentation, matrix export, residual checks, and first PCG experiments.

The TS-N_65 case has roughly 600k DOF and is the hard target. The email reports that only Pardiso is currently usable, with multi-day runtime and about 20 GB RAM. Plain CG is not expected to be competitive; any memory reduction depends on a strong preconditioner.

Use two separate baselines:

- Production baseline: `PardisoLDLT`, explicit `MKL_NUM_THREADS` and `OMP_NUM_THREADS`, same stiffness update policy as target runs.
- Algorithmic baseline: current `EigenConj`, instrumented for iterations and residuals.

## Recommended Solver Path

1. Add solver instrumentation first: matrix size, nonzeros, pattern hash, analyze/factorize/solve counts, CG iterations, final residual, failure status, and timing.
2. Add an incomplete-Cholesky PCG variant for Dogbone-level validation.
3. Prototype shared-memory AMG as a PCG preconditioner for the large case. AMGCL is the lowest-friction first integration; hypre/PETSc are stronger long-term options but more intrusive.
4. Add deflated or recycled PCG once the base preconditioner is stable.

For vector elasticity, AMG should receive DOF-block and near-nullspace information where possible. Black-box scalar AMG may run but can converge poorly on mechanics problems.

## Deflation Design

The OAS `LinAlgSolver` lifecycle already has useful hooks:

- `analyzePattern(A)`: cache graph metadata and allocate solver structures.
- `factorize(A)`: build or rebuild preconditioner state.
- `solve(x, b)`: run the current solve and record diagnostics.

A deflated PCG solver can store a small basis `Z`, compute `AZ` and `E = Z^T A Z` during factorization, and use a two-level balanced preconditioner in the PCG loop. For TS-N_65, storing 20 vectors and their images costs roughly 183 MiB, which is small compared with the reported direct-solver memory use.

The important engineering issue is vector lifecycle. Accepted increments or accepted Newton corrections are better history vectors than blindly storing every returned solve vector. Rejected nonlinear attempts and indirect-control right-hand sides should be labeled or filtered.

## Risks

- `analyzePattern()` is currently called once per solver object. Add a graph hash and rerun analysis if sparsity changes.
- Current CG warm-starts from the last solution vector, which can be misleading for semantically different right-hand sides.
- Iterative solver failure should return status and allow fallback, not terminate the run abruptly.
- Deflation can hurt when the basis is stale, rank deficient, or dominated by nearly collinear loading history.
