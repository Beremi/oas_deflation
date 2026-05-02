# OAS Deflation Workspace

This repository is an at-root fork of [Open Academic Solver (OAS)](https://gitlab.com/kelidas/OAS) for experiments with preconditioned and deflated Krylov linear solvers. The upstream OAS source is kept directly in this repository so edits to OAS files show up in ordinary `git diff` output.

The current upstream baseline is `upstream/master` at `3ce5a809c5cdc19bd057e96efac2e47bc167bd68`, matching the reviewed OAS snapshot from 2026-04-27. The generated OAS documentation is at <https://kelidas.gitlab.io/OAS-site/>.

## Repository Model

- `upstream` points to the public OAS repository.
- `main` tracks `upstream/master`.
- Deflation project docs live in `docs/deflation/`.
- Local benchmark inputs and run outputs live under `data/` and are ignored.
- Local agent working context lives in `agents.md` and `.agent_context/`; both are ignored.

Use this sync flow before starting solver work:

```sh
git fetch upstream
git rebase upstream/master
```

Use merge instead of rebase only for branches that have already been shared.

## Build

The root `Makefile` wraps the OAS out-of-source CMake build. The default build directory is:

```text
../oas_deflation-build/release
```

Configure and build:

```sh
make configure
make build
```

The default configuration enables VTK, Intel MKL Pardiso, and CHOLMOD. A smaller diagnostic configuration is available when optional solver dependencies are the suspected problem:

```sh
make configure-fast
```

If VTK discovery fails because the system VTK package is missing a transitive CMake package such as CLI11, keep the solver libraries enabled and disable only VTK:

```sh
make configure USE_VTK=OFF
make build USE_VTK=OFF
```

The Makefile sets `OAS_CXX_FLAGS=-DEIGEN_MKL_NO_DIRECT_CALL` by default. This keeps MKL/Pardiso available while avoiding Eigen/MKL direct-call compilation failures with GCC 15.

## Dogbone Benchmark

The Dogbone archive is local-only and ignored in `data/archives/Dogbone.zip`. It is extracted to `data/cases/Dogbone/`.

Run the benchmark with the solver already active in `solver.inp`:

```sh
make dogbone
```

Override the active solver in the ignored Dogbone input deck and run:

```sh
make dogbone-solver SOLVER=EigenLDLT
make dogbone-solver SOLVER=PardisoLDLT
```

Set runtime threads with:

```sh
make dogbone THREADS=8
```

For MKL/Pardiso runs, the Makefile sets both `MKL_NUM_THREADS` and `OMP_NUM_THREADS`.

Profile the current linear-solve behavior with local-only outputs:

```sh
make dogbone-profile USE_VTK=OFF THREADS=4 SOLVER=EigenLDLT
```

This temporarily enables `linear_solver_profile 1` in the ignored Dogbone `solver.inp`, writes raw TSV files under `data/cases/Dogbone/results/`, restores the input file, and generates a Markdown/PNG exchange report under `results/dogbone-<solver>-<timestamp>/`. The `results/` directory is excluded through `.git/info/exclude` so quick experiment reports stay out of project history.

## First Solver Goals

The immediate engineering target is not to replace Pardiso with plain CG. The reviewed OAS source already has `EigenConj`, but it is a weak Jacobi-preconditioned CG baseline. The staged solver plan is:

1. Instrument current linear solver calls, matrix updates, iteration counts, residuals, and matrix pattern stability.
2. Add an incomplete-Cholesky PCG variant for Dogbone-level validation.
3. Prototype AMG preconditioning for the larger TS-N_65 case.
4. Add history-based deflation/recycling on top of a stable base preconditioner.

See `docs/deflation/implementation-report.md` for the compact technical report and `docs/deflation/oas_linear_solver_review_parallel_amended.tex` for the full source note.
