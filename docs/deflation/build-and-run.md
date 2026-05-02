# Build and Run Notes

## Local Prerequisites

This workspace has been prepared for a Make-based CMake build. Expected local tools and libraries:

- CMake 3.18 or newer
- C++20 compiler
- GNU Make
- ccache
- VTK CMake package at `/usr/lib/cmake/vtk`
- Intel oneAPI MKL at `/opt/intel/oneapi/mkl/latest`
- CHOLMOD/SuiteSparse headers and libraries

Ninja is not assumed to be installed.

## Standard Build

```sh
make configure
make build
```

The default build directory is outside the source tree:

```text
../oas_deflation-build/release
```

The expected solver executable is:

```text
../oas_deflation-build/release/bin/OAS
```

The default configuration keeps VTK enabled. If the local VTK package is incomplete, for example CMake reports a missing `CLI11Config.cmake`, keep Pardiso and CHOLMOD enabled but bypass VTK for solver development:

```sh
make configure USE_VTK=OFF
make build USE_VTK=OFF
```

The Makefile also defines `EIGEN_MKL_NO_DIRECT_CALL` by default. This keeps Eigen's MKL/Pardiso support enabled while avoiding MKL direct-call macro compilation failures seen with GCC 15.

## Fast Diagnostic Configure

If configuration fails in optional solver or VTK detection, try:

```sh
make configure-fast
```

This creates:

```text
../oas_deflation-build/fast
```

It disables Pardiso, CHOLMOD, VTK, and docs. It is only for isolating build problems, not for production solver comparison.

## Dogbone

The benchmark data is ignored and local:

```text
data/archives/Dogbone.zip
data/cases/Dogbone/
```

Run:

```sh
make dogbone
```

Run with a solver override:

```sh
make dogbone-solver SOLVER=EigenLDLT
make dogbone-solver SOLVER=PardisoLDLT
```

Thread count:

```sh
make dogbone THREADS=8
```

For Pardiso, the Makefile sets `MKL_NUM_THREADS` and `OMP_NUM_THREADS`. The `-j` OAS argument is also passed for OpenMP-controlled paths.

## Linear-Solve Profiling

Enable the optional OAS profiler for a Dogbone run and generate a local exchange report:

```sh
make dogbone-profile USE_VTK=OFF THREADS=4 SOLVER=EigenLDLT
```

The target temporarily adds these keys to the ignored Dogbone `solver.inp` and restores the original file after the run:

```text
linear_solver_profile 1
linear_solver_profile_matrix_delta 1
linear_solver_profile_file linear_profile
```

Raw profile data stays with the OAS case output:

```text
data/cases/Dogbone/results/linear_profile_events.tsv
data/cases/Dogbone/results/linear_profile_iterations.tsv
```

The human-readable Markdown report and PNG plots are written under ignored local `results/`:

```text
results/dogbone-eigenldlt-<timestamp>/linear-profile.md
```

Regenerate a report from an existing profile without re-running OAS:

```sh
make linear-profile-report PROFILE_DIR=data/cases/Dogbone/results OUT_DIR=results/dogbone-eigenldlt-manual SOLVER=EigenLDLT
```

## Upstream Sync

```sh
make sync-upstream
```

This runs:

```sh
git fetch upstream
git rebase upstream/master
```

Use merge instead of rebase only when working on a branch already shared with others.
