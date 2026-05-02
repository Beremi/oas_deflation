# Benchmark Notes

## Dogbone

Dogbone is the fast development case. The archive contains the input deck plus reference results produced by an OAS build from `6b734518a4-dirty` on 2026-04-27.

Local path:

```text
data/cases/Dogbone/
```

The provided `solver.inp` uses:

```text
SteadyStateNonLinearSolver
time_step 0.01
total_time 1
solver_type EigenLDLT
stiffness_matrix_iter_update 20
stiffness_matrix_step_update -1
```

The archived reference `solver.out` reaches step 100 at `time = 1.0`. Use that as a shape check, not as an exact timing baseline, because local compiler, solver backend, thread count, and OAS revision can differ.

Recommended first runs:

```sh
make dogbone-solver SOLVER=EigenLDLT
make dogbone-solver SOLVER=PardisoLDLT
make dogbone-solver SOLVER=EigenConj
```

## TS-N_65

TS-N_65 is the hard benchmark and should not be extracted by default. The archive is local and ignored:

```text
data/archives/TS-N_65.zip
```

The input deck reports:

```text
SteadyStateNonLinearSolver
time_step 5.000000e-03
total_time 1
solver_type PardisoLDLT
stiffness_matrix_iter_update -5
stiffness_matrix_step_update 1
tolerance 1e-3
```

The email states approximately 600k DOF, multi-day runtime, and roughly 20 GB RAM with Pardiso. Do not run this benchmark casually from automation.

## Data Policy

Benchmark archives, extracted inputs, and generated results stay under `data/`, which is ignored. Commit scripts and small documentation that describe how to obtain or run the data; do not commit the large binary archives.
