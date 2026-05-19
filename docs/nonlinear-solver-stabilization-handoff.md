# Nonlinear Solver Stabilization Handoff

## Purpose

This note is a self-contained handoff for improving the nonlinear stabilization
of the steady-state solver used by the TS-N65 LDPM benchmark. It is written for
an expert who can advise on nonlinear finite element / discrete element
globalization strategies, especially for softening damage and reinforcement
plasticity.

The immediate problem is not that the linear solver cannot solve individual
linear systems. The observed failure mode is that some load steps require
hundreds to thousands of nonlinear iterations, and some are accepted only through
the fallback `limit_tolerance` path after failing the requested convergence
tolerances.

## Problem Being Solved

The TS-N65 case uses `SteadyStateNonLinearSolver` in 3D mechanics. The global
unknown is the reduced displacement vector `u` over free mechanical degrees of
freedom. At each pseudo-time/load step, the code solves static equilibrium:

```text
R(u, lambda) = f_ext(lambda) - f_int(u, history) = 0
```

The material model is path dependent. The internal force depends on the current
trial displacement and on material history variables such as damage, crack
opening, plastic strain, and cumulative plastic strain.

The current nonlinear iteration is a Newton-like fixed-point loop:

```text
assemble/update K
form residual r = f_ext - f_int
solve K * du = r
apply full du to the trial field
recompute internal forces and convergence errors
repeat until residual, displacement increment, and energy errors pass
```

In the steady-state direct-control path, the increment is always applied at full
length. There is no currently implemented globalization step between solving
`K * du = r` and committing `du`.

Relevant code:

- `src/solver/src/solver_implicit.cpp`
- `SteadyStateNonLinearSolver::solve()`
- `SteadyStateLinearSolver::updateSystemMatrices()`
- `SteadyStateNonLinearSolver::evaluateErrors()`

## Convergence Metrics

The solver evaluates three normalized errors:

```text
resErr = normalized residual norm
disErr = normalized displacement increment norm
eneErr = normalized residual work / energy norm
```

Convergence requires all active criteria to be below their tolerances:

```text
disErr <= maxDisErr
resErr <= maxResErr
eneErr <= maxEneErr
```

The default `tolerance` keyword sets all three tolerances. Separate controls are:

```text
tolerance_residuals
tolerance_energies
tolerance_increments
```

Important: `limit_tolerance` is not a multiplier on `tolerance`. It is parsed as
an absolute fallback threshold:

```text
limitEneErr = limitResErr = limitDisErr = input value
```

So `limit_tolerance 10` means a step can be accepted after exhausting
`max_iterations` as long as all normalized errors are below `10`, even if the
requested tolerance is `1e-3`.

## Current Local Capabilities

The following stabilization-related controls are implemented locally:

```text
max_iterations / maxIt
min_iterations / minIt
time_step
min_time_step
max_time_step
step_increase
step_decrease
critical_step_decrease
enlargeIt
shortenIt
limit_tolerance
stiff_matrix_type
first_iteration_stiff_matrix_type
stiffness_matrix_iter_update
stiffness_matrix_cumul_iter_update
stiffness_matrix_step_update
indirect_control / indirect_displacement_control
indirect_control_square_sum
```

Adaptive stepping exists, but it only acts after a step has failed or after a
converged step used many/few iterations. It does not reduce the Newton increment
inside the current iteration.

`indirect_control` changes the load multiplier to hit a prescribed displacement
or force-related control measure. It is useful for some unstable equilibrium
paths, but it is not a general Newton line search.

`stiff_matrix_type` can be `elastic`, `secant`, `tangent`, or `consistent`, but
material support varies. For `CSLMaterial`, `tangent` is not a true consistent
tangent; the code currently returns unloading/degraded stiffness for that branch.

## Current Missing Local Capabilities

The following are not implemented in the steady-state nonlinear solver:

```text
Newton line search
backtracking residual decrease check
constant Newton damping factor
adaptive Newton damping factor
trust-region method
Levenberg-Marquardt style regularization
arc-length / Riks method
early stagnation detection
automatic cutback before max_iterations is exhausted
trial alpha evaluation with guaranteed material-state rollback
```

There are `damping_matrix_*` options in transient solvers, but these refer to
damping/capacity matrices in time integration, not to damping of Newton
increments in `SteadyStateNonLinearSolver`.

There are `hypre_relax_*` and `dfgmres_restart` options, but those are linear
solver/preconditioner controls, not nonlinear stabilization.

## Observed TS-N65 Pathology

The active TS-N65 run used:

```text
time_step 1.25e-3
total_time 1e-1
max_iterations 1000
stiffness_matrix_iter_update 10
stiffness_matrix_step_update 1
limit_tolerance 10
tolerance 1e-3
solver_type DeflatedFGMRES
```

No `min_time_step` or `max_time_step` was provided, so the code used a fixed
step size with `dtmin == dtmax == dt`. That disables the useful cutback behavior.

Observed summary from `results/solver.out`:

```text
steps analyzed: 78
average nonlinear iterations: 353.9
steps above 100 nonlinear iterations: 68
steps hitting 1000 nonlinear iterations: 7
accepted through fallback with converged = 0: 7
```

Example failing/slow steps:

```text
step 31: 1000 iterations, resErr 1.65e-3, eneErr 2.22e-3, converged 0
step 32: 1000 iterations, resErr 1.27e-3, eneErr 1.66e-3, converged 0
steps 74-78: 1000 iterations each, converged 0
```

These steps hover close to the requested `1e-3` tolerance but do not robustly
cross all criteria. The current settings allow a full 1000 iterations before any
response, then accept through `limit_tolerance`.

## What Can Be Done Immediately With Existing Knobs

A safer baseline for testing nonlinear stabilization behavior is:

```text
max_iterations 40
min_time_step 1.25e-5
max_time_step 1.25e-3
critical_step_decrease 0.5
step_decrease 0.5
step_increase 1.2
enlargeIt 6
shortenIt 15
limit_tolerance 0
```

Matrix-update tests worth running:

```text
stiffness_matrix_iter_update 1
stiffness_matrix_iter_update 3
stiffness_matrix_iter_update 10
```

This does not solve globalization, but it prevents a single bad increment from
burning 1000 iterations and prevents unconverged acceptance through a very loose
fallback tolerance.

## Implementation Questions For Expert Guidance

The main implementation decision is which nonlinear globalization method is most
appropriate for this problem class:

```text
1. Backtracking line search on the Newton increment.
2. Fixed/adaptive damping of du.
3. Trust region or regularized Newton.
4. Arc-length / Riks / indirect displacement control for snapback.
5. Viscous or rate regularization for softening damage.
6. Material-level substepping versus global load-step cutback.
```

The simplest local code change is probably a residual-based line search:

```text
du_full = linear solve
for alpha in [1, 0.5, 0.25, ...]:
    trial u = u_old_iter + alpha * du_full
    evaluate residual/error with temporary material state
    accept alpha if merit decreases sufficiently
```

The hard part is not the scalar search itself. The hard part is making trial
evaluations safe with path-dependent material state. The code must be audited so
that rejected alpha trials do not accidentally advance damage/plasticity history
or leave temporary variables in a state inconsistent with the accepted alpha.

Questions needing expert design:

```text
What merit function should be used: residual norm, energy norm, mixed norm, or
incremental potential?

Should the line search allow material damage/plasticity to evolve during trial
residual evaluations, or should it freeze history variables?

What rollback boundary is required between trial alpha values?

Should a failed line search trigger global load-step cutback immediately?

How should line search interact with `limit_tolerance`?

Is load control fundamentally inappropriate in the softening portions of TS-N65,
requiring displacement/arc-length control instead?
```

## Global Limitations Independent Of This Code

For softening damage and localization, there may be load steps where the
load-controlled equilibrium path is unstable or has snapback. In those regions,
line search can improve robustness but cannot create a physically stable
load-controlled equilibrium if the mathematical problem has lost stability.

For brittle or quasi-brittle damage, the tangent may become non-positive or
strongly nonsymmetric/effectively inconsistent with the residual. A better
linear solve alone does not guarantee nonlinear convergence.

For mesh-sensitive localization, nonlinear solver stabilization is not a
replacement for physical or numerical regularization. The expert should consider
whether the material model and mesh require fracture-energy regularization,
nonlocal/gradient damage, viscosity, or controlled crack-band behavior.

For severe post-peak response, robust path following may require displacement
control, indirect control, or arc-length methods. A local Newton line search is
still useful, but it is not a complete continuation strategy.

## External References Useful For Orientation

PETSc SNES manual:

```text
https://petsc.org/main/manual/snes/
```

OpenSees Newton line search:

```text
https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/algorithm/NewtonLineSearch.html
```

Abaqus static nonlinear incrementation:

```text
https://abaqus.uclouvain.be/English/SIMACAEANLRefMap/simaanl-c-static.htm
```

## Recommended Next Patch Scope

A focused first implementation should add:

```text
nonlinear_line_search 0/1
nonlinear_line_search_min_alpha
nonlinear_line_search_reduction
nonlinear_line_search_max_trials
nonlinear_line_search_merit residual|energy|mixed
nonlinear_stagnation_cutback 0/1
nonlinear_stagnation_iterations
nonlinear_stagnation_ratio
```

The first patch should avoid changing material laws. It should only add a safe
trial-evaluation mechanism, rollback semantics, and an early cutback path. Once
that exists, TS-N65 can be used to compare fixed-step, adaptive-step,
line-search, and indirect-control behavior.
