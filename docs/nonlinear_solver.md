# OAS nonlinear solver enhancement handout

**Repository:** `Beremi/oas_deflation`, branch `original`  
**Primary solver target:** `SteadyStateNonLinearSolver` in `src/solver/src/solver_implicit.cpp` / `src/solver/src/solver_implicit.h`  
**Primary benchmark pressure case:** TS-N65 LDPM steady-state mechanics case  
**Goal:** reduce extremely slow nonlinear convergence, prevent unconverged fallback acceptance, and add robust globalization for path-dependent softening/damage/plasticity.

---

## 1. Executive recommendation

The most practical implementation path is not to jump directly to arc-length or rewrite material laws. The solver already has load-step cutback, stiffness-update options, frozen force evaluation, and material-status reset/update hooks. The missing layer is a **nonlinear globalization layer** between solving `K * du = r` and accepting the displacement increment.

Implement in this order:

1. **Safety baseline and instrumentation**: disable loose fallback acceptance, enable adaptive cutback, log per-iteration metrics, and run a controlled matrix-update sweep.
2. **Early stagnation detection with immediate load-step cutback**: stop burning hundreds or thousands of iterations once the iteration has clearly stalled.
3. **Constant/adaptive Newton damping**: cheapest code change; useful as a control and fallback.
4. **Backtracking Newton line search using frozen material-state trial residuals**: highest-value solver enhancement; reuse existing frozen force evaluation where possible.
5. **Material-state snapshot/rollback for trial evaluations**: required for robust line search if trial stress evaluations mutate temporary state even when `frozen=true`.
6. **Tangent/matrix-update improvements**: update stiffness more often and audit material tangent consistency.
7. **Trust-region / Levenberg-Marquardt-style regularization**: second-line globalization for indefinite or highly inconsistent tangents.
8. **Indirect displacement control and then arc-length/Riks continuation**: needed if load control is mathematically inappropriate in snapback regions.
9. **Material substepping, viscosity, and localization regularization**: broader modeling changes for damage/softening pathologies not solvable by Newton globalization alone.
10. **External nonlinear solver integration, for example PETSc SNES**: powerful but intrusive; use after the local globalization API is clean.

The first mergeable patch should add **line search + stagnation cutback + safe trial evaluation hooks**, without changing constitutive laws.

---

## 2. Current solver diagnosis

### 2.1 Current nonlinear loop

The steady-state nonlinear solve is Newton-like:

```text
assemble/update K
form residual r = f_ext - f_int
solve K * du = r
apply full du to the trial displacement field
recompute internal forces and convergence errors
repeat until residual, displacement increment, and energy errors pass
```

The current direct-control path applies the full Newton increment. There is no implemented line search, trust region, adaptive damping, or trial-step merit check between the linear solve and displacement update.

### 2.2 Existing convergence metrics

The solver tracks three normalized errors:

```text
resErr = residual error
disErr = displacement increment error
eneErr = residual-work / energy error
```

A nonlinear iteration is accepted only when all active criteria pass:

```text
disErr <= maxDisErr
resErr <= maxResErr
eneErr <= maxEneErr
```

The important practical issue is that `limit_tolerance` is parsed as an **absolute fallback threshold**, not as a multiplier on `tolerance`. Therefore:

```text
tolerance 1e-3
limit_tolerance 10
```

means the solver may accept a step after exhausting `max_iterations` if normalized errors are below `10`, even though the requested tolerance is `1e-3`.

### 2.3 Observed TS-N65 pathology

The supplied TS-N65 run used approximately:

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

No `min_time_step` or `max_time_step` was specified, so the code effectively used a fixed step size. That disables useful adaptive cutback behavior.

Observed behavior:

```text
steps analyzed: 78
average nonlinear iterations: 353.9
steps above 100 nonlinear iterations: 68
steps hitting 1000 nonlinear iterations: 7
accepted through fallback with converged = 0: 7
```

This indicates a nonlinear globalization problem, not primarily a linear-solver failure. Many steps hover near the requested tolerance but do not robustly cross it.

---

## 3. Relevant code map

### 3.1 Primary files

```text
src/solver/src/solver_implicit.cpp
src/solver/src/solver_implicit.h
src/solver/src/solver.cpp
src/solver/src/solver.h
src/solver/src/element_container.cpp
src/solver/src/element_container.h
```

### 3.2 Important methods and hooks

#### `SteadyStateNonLinearSolver::solve()`

Main nonlinear loop. The key insertion point is immediately after:

```cpp
linalgsolver->solve(ddr, f);
```

and before the current full-step update:

```cpp
updateFieldVariables();
computeForcesAtIntegrationTime(true);
evaluateErrors();
```

At that point, `ddr` contains the reduced Newton increment. The globalization layer should decide whether to apply:

```text
1.0 * ddr
alpha * ddr
some trust-region-limited ddr
or no increment, followed by load-step cutback
```

#### `Solver::updateFieldVariables()`

Maps reduced `ddr` to full `full_ddr`, handles dependent DOFs, and increments `trial_r`:

```cpp
nodes->giveFullDoFArray(ddr, full_ddr);
nodes->updateFullDoFsByDependenciesOnConjugates(full_ddr, trial_r, f_ext);
trial_r.add(full_ddr);
```

Do **not** naively call this repeatedly inside a trial line-search loop unless state is restored between trials. It mutates `trial_r` by adding the current increment.

#### `computeInternalExternalForces(...)`

Computes internal and external force vectors and the residual. It supports a `frozen` flag:

```cpp
computeInternalExternalForces(trial_r, load, frozen, time, timeStep, f_int, f_ext);
```

For trial line search, prefer frozen material-variable evaluation first:

```cpp
computeForcesAtIntegrationTime(true);
```

or the lower-level form with `frozen=true`.

#### `ElementContainer::updateMaterialStatuses()` and `resetMaterialStatuses()`

The code exposes material-status commit/reset hooks. These are essential for safe trial evaluations:

```cpp
elems->updateMaterialStatuses();
elems->resetMaterialStatuses();
```

The line-search implementation must confirm exactly what “reset” means for every material status class used by LDPM/CSL/reinforcement models.

---

## 4. Prioritized option table

| Priority | Option                                               |              Promise | Implementation cost | Use first when                                                             |
| -------: | ---------------------------------------------------- | -------------------: | ------------------: | -------------------------------------------------------------------------- |
|        0 | Instrumentation + safer existing input baseline      |            Very high |            Very low | Always                                                                     |
|        1 | Early stagnation detection and automatic cutback     |            Very high |                 Low | Iterations plateau or oscillate                                            |
|        2 | Fixed/adaptive Newton damping                        |                 High |                 Low | Need immediate stabilizer and line-search baseline                         |
|        3 | Backtracking line search with frozen trial residuals |              Highest |              Medium | Newton direction usually good but full step overshoots                     |
|        4 | Trial material-state snapshot/rollback               |            Very high |         Medium-high | Trial evaluations mutate status variables                                  |
|        5 | Stiffness update/tangent consistency improvements    |                 High |      Medium to high | Residual decreases slowly due to stale/inconsistent K                      |
|        6 | Trust-region / regularized Newton                    |          Medium-high |              Medium | Tangent is indefinite, poor, or nonlinear residual grows                   |
|        7 | Nonlinear acceleration / quasi-Newton                |               Medium |              Medium | Matrix formation is expensive and residual behavior is smooth enough       |
|        8 | Better indirect displacement control                 |          Medium-high |              Medium | Load-control path is unstable but a displacement control variable is known |
|        9 | Arc-length / Riks continuation                       |    High for snapback |                High | Need post-peak/snapback path following                                     |
|       10 | Material-level substepping                           |          Medium-high |                High | Local constitutive updates overshoot or local plasticity/damage jumps      |
|       11 | Viscous/rate regularization                          |          Medium-high |       High/modeling | Damage localization causes ill-posed quasi-static solve                    |
|       12 | Nonlocal/gradient/crack-band regularization          |        High/modeling |           Very high | Mesh-sensitive localization dominates                                      |
|       13 | Linear solver/preconditioner tuning                  | Low for this failure |          Low-medium | Linear solves fail or are inaccurate                                       |
|       14 | PETSc SNES integration                               |   High but long-term |           Very high | Want mature nonlinear methods and can afford API refactor                  |

---

# Option 0 — Establish safety baseline and instrumentation

## Why this is first

Before adding algorithms, prevent misleading results. The current TS-N65 setup allows fixed-step iteration up to `max_iterations=1000` and then loose fallback acceptance through `limit_tolerance 10`. That makes nonlinear improvements hard to evaluate.

## Immediate input-file baseline

Use this baseline for testing all nonlinear enhancements:

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

Run a stiffness-update sweep:

```text
stiffness_matrix_iter_update 1
stiffness_matrix_iter_update 3
stiffness_matrix_iter_update 10
```

Recommended first TS-N65 matrix test matrix:

| Case | `max_iterations` | `min_time_step` | `limit_tolerance` | `stiffness_matrix_iter_update` |
| ---- | ---------------: | --------------: | ----------------: | -----------------------------: |
| A    |               40 |       `1.25e-5` |                 0 |                              1 |
| B    |               40 |       `1.25e-5` |                 0 |                              3 |
| C    |               40 |       `1.25e-5` |                 0 |                             10 |
| D    |               80 |       `1.25e-5` |                 0 |                              1 |

## Add instrumentation

### New log fields

Add these per nonlinear iteration:

```text
step
iteration
time
dt
resErr
disErr
eneErr
merit
alpha
line_search_trials
stiffness_rebuilt
cutback_reason
converged
fully_converged
linear_solver_iterations if available
linear_solver_residual if available
```

### Suggested output line

```text
NL step=%d it=%d t=%+.8e dt=%+.3e res=%+.3e dis=%+.3e ene=%+.3e merit=%+.3e alpha=%+.3e ls=%d K=%d reason=%s
```

### Implementation notes

Add fields to `SteadyStateNonLinearSolver`:

```cpp
double lastLineSearchAlpha = 1.0;
unsigned lastLineSearchTrials = 0;
bool lastStiffnessRebuilt = false;
std::string lastCutbackReason = "";
```

Add a helper:

```cpp
void SteadyStateNonLinearSolver::printNonlinearIterationLog(
    unsigned step,
    unsigned it,
    double merit,
    bool stiffnessRebuilt,
    const char* reason);
```

## Checklist

- [ ] Disable `limit_tolerance` or set it to `0` in benchmark inputs.
- [ ] Ensure `dtmin < dtmax` so adaptive cutback can actually operate.
- [ ] Reduce `max_iterations` from 1000 to a diagnostic value such as 40 or 80.
- [ ] Log `alpha`, line-search trials, and cutback reason even before line search is implemented.
- [ ] Save `solver.out` for each case under a unique run folder.
- [ ] Parse and summarize: mean iterations, max iterations, number of cutbacks, number of fallback accepts, final displacement/load curve.

## Acceptance criteria

- No step is accepted through a loose fallback tolerance.
- A failed step cuts back rather than spending 1000 iterations.
- You can compare nonlinear strategies using the same metrics.

---

# Option 1 — Early stagnation detection and automatic cutback

## Promise

This is the highest-return low-risk patch. It does not attempt to make a bad Newton step good; it stops the solver when the current load increment is clearly too large.

## New input keywords

```text
nonlinear_stagnation_cutback 0/1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25
nonlinear_cutback_after_nan 1
nonlinear_cutback_after_line_search_fail 1
```

Suggested defaults:

```text
nonlinear_stagnation_cutback 1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25
nonlinear_cutback_after_nan 1
nonlinear_cutback_after_line_search_fail 1
```

## Merit function for stagnation

Use a mixed normalized merit:

```cpp
phi = max(
    activeResidualCriterion ? resErr / maxResErr : 0.0,
    activeEnergyCriterion   ? eneErr / maxEneErr : 0.0
);
```

Do not use `disErr` as the primary stagnation merit because damping or line search can reduce it trivially by shrinking `du`. Keep `disErr` as a convergence criterion, but not as the main globalization merit.

Optionally include displacement increment as a capped secondary term:

```cpp
phi = max(phi, 0.1 * disErr / maxDisErr);
```

## Stagnation rules

Let `phi[k]` be the merit after iteration `k`.

Trigger cutback if any of these hold after at least `nonlinear_stagnation_iterations` iterations:

```text
1. phi[k] > nonlinear_growth_cutback * phi[k-1]
2. min(phi over last N iterations) > nonlinear_stagnation_ratio * phi[k-N]
3. phi oscillates without reducing below best_phi * nonlinear_stagnation_ratio
4. NaN or Inf appears in any error
```

The simplest robust rule:

```cpp
bool stagnated = false;
if (it >= stagnationIterations) {
    double oldPhi = meritHistory[meritHistory.size() - stagnationIterations];
    double recentBest = min_last_N(meritHistory, stagnationIterations);
    stagnated = recentBest > stagnationRatio * oldPhi;
}
```

Interpretation for `stagnationRatio=0.95`: if the best recent merit has not improved by at least 5% over `N` iterations, cut back.

## Implementation guide

### Step 1 — Add fields in `solver_implicit.h`

```cpp
bool nonlinearStagnationCutback = true;
unsigned nonlinearStagnationIterations = 8;
double nonlinearStagnationRatio = 0.95;
double nonlinearGrowthCutback = 1.25;
std::vector<double> nonlinearMeritHistory;
std::string nonlinearCutbackReason;
```

### Step 2 — Parse fields in `readFromFile`

Follow the existing parser style used for `max_iterations`, `min_time_step`, `stiffness_matrix_iter_update`, and `limit_tolerance`.

Example pattern:

```cpp
else if ( param.compare("nonlinear_stagnation_cutback") == 0 ) {
    ifs >> nonlinearStagnationCutback;
}
else if ( param.compare("nonlinear_stagnation_iterations") == 0 ) {
    ifs >> nonlinearStagnationIterations;
}
else if ( param.compare("nonlinear_stagnation_ratio") == 0 ) {
    ifs >> nonlinearStagnationRatio;
}
else if ( param.compare("nonlinear_growth_cutback") == 0 ) {
    ifs >> nonlinearGrowthCutback;
}
```

### Step 3 — Add helper methods

```cpp
double SteadyStateNonLinearSolver::currentNonlinearMerit() const {
    double phi = 0.0;
    if (maxResErr > 0.0) phi = std::max(phi, resErr / maxResErr);
    if (maxEneErr > 0.0) phi = std::max(phi, eneErr / maxEneErr);
    return phi;
}

bool SteadyStateNonLinearSolver::shouldCutbackForStagnation(double phi) {
    nonlinearMeritHistory.push_back(phi);

    if (!nonlinearStagnationCutback) return false;
    if (nonlinearMeritHistory.size() <= nonlinearStagnationIterations) return false;

    const size_t n = nonlinearMeritHistory.size();
    const double previous = nonlinearMeritHistory[n - nonlinearStagnationIterations - 1];
    double recentBest = nonlinearMeritHistory[n - nonlinearStagnationIterations];
    for (size_t i = n - nonlinearStagnationIterations; i < n; ++i) {
        recentBest = std::min(recentBest, nonlinearMeritHistory[i]);
    }

    if (recentBest > nonlinearStagnationRatio * previous) {
        nonlinearCutbackReason = "stagnation";
        return true;
    }
    if (n >= 2 && nonlinearMeritHistory[n-1] > nonlinearGrowthCutback * nonlinearMeritHistory[n-2]) {
        nonlinearCutbackReason = "merit_growth";
        return true;
    }
    return false;
}
```

### Step 4 — Integrate into `solve()`

After `evaluateErrors()`:

```cpp
double phi = currentNonlinearMerit();
if (shouldCutbackForStagnation(phi)) {
    restart_now = true;
    break;
}
```

Then let the existing cutback branch handle the time-step reduction. If the current code only cuts back after the while-loop exits with `converged=false`, set `converged=false`, set a `restart_now`/`force_cutback` flag, and break out of the nonlinear iteration.

### Step 5 — Reset history on every new load-step attempt

At the start of each step or restart:

```cpp
nonlinearMeritHistory.clear();
nonlinearCutbackReason.clear();
```

## Checklist

- [ ] Add new fields and parse input options.
- [ ] Use residual/energy merit, not displacement increment alone.
- [ ] Reset merit history after each load-step restart.
- [ ] Emit cutback reason in logs.
- [ ] Verify cutback happens before `max_iterations` on intentionally oversized increments.
- [ ] Confirm no material state is committed before the cutback branch restores the previous step state.

## Acceptance criteria

- TS-N65 no longer spends hundreds of iterations on a single stalled increment.
- Failed steps cut back before `max_iterations`.
- Results are independent of arbitrary `max_iterations=1000` settings.

---

# Option 2 — Constant and adaptive Newton damping

## Promise

Damping is the cheapest globalization method. It is less powerful than line search but easy to implement and useful as a fallback.

Instead of applying:

```text
u_{k+1} = u_k + du
```

apply:

```text
u_{k+1} = u_k + alpha * du
```

with either a fixed or adaptive `alpha`.

## New input keywords

```text
nonlinear_damping_type off|fixed|adaptive
nonlinear_damping_factor 0.5
nonlinear_damping_min 0.03125
nonlinear_damping_max 1.0
nonlinear_damping_increase 1.25
nonlinear_damping_decrease 0.5
```

Suggested defaults:

```text
nonlinear_damping_type off
nonlinear_damping_factor 1.0
nonlinear_damping_min 0.03125
nonlinear_damping_max 1.0
nonlinear_damping_increase 1.25
nonlinear_damping_decrease 0.5
```

For experiments:

```text
nonlinear_damping_type fixed
nonlinear_damping_factor 0.5
```

or:

```text
nonlinear_damping_type adaptive
nonlinear_damping_factor 0.5
```

## Implementation guide

### Step 1 — Store full Newton increment before damping

After the linear solve:

```cpp
linalgsolver->solve(ddr, f);
```

copy `ddr`:

```cpp
Vector ddrNewton = ddr;
```

### Step 2 — Apply damping before `updateFieldVariables()`

```cpp
double alpha = currentDampingFactor;
ddr = ddrNewton;
ddr.times(alpha);          // or equivalent vector scaling method
lastLineSearchAlpha = alpha;
updateFieldVariables();
computeForcesAtIntegrationTime(true);
evaluateErrors();
```

Use the repository's actual vector API for scaling. If no `times` method exists, add a small helper to scale vectors consistently.

### Step 3 — Adaptive damping update

After evaluating errors:

```cpp
double phi = currentNonlinearMerit();

if (phi < previousPhi) {
    currentDampingFactor = std::min(nonlinearDampingMax,
                                    nonlinearDampingIncrease * currentDampingFactor);
} else {
    currentDampingFactor = std::max(nonlinearDampingMin,
                                    nonlinearDampingDecrease * currentDampingFactor);
}
previousPhi = phi;
```

This is a crude but useful control. Prefer line search for production because it explicitly tests the trial step before accepting it.

## Checklist

- [ ] Add parser options.
- [ ] Initialize `currentDampingFactor` at the start of each load step.
- [ ] Log the applied alpha.
- [ ] Keep convergence criteria unchanged.
- [ ] Run fixed-alpha tests: `1.0`, `0.75`, `0.5`, `0.25`, `0.125`.
- [ ] Compare total nonlinear iterations and load-step cutbacks.

## Risks

- Damping can make convergence slow if alpha is too small.
- Displacement-increment error may look artificially better with small alpha.
- Damping does not guarantee residual decrease.
- Damping cannot follow unstable load-control snapback paths.

## Acceptance criteria

- The damped solver reduces divergence/NaN cases.
- Fixed damping does not hide unconverged steps through fallback tolerance.
- The best fixed damping factor gives a baseline for line-search alpha behavior.

---

# Option 3 — Backtracking Newton line search with frozen trial residuals

## Promise

This is the most promising local solver enhancement.

The Newton direction is often good, but the full step is too long. Backtracking reuses the expensive linear solve and evaluates several cheaper trial residuals:

```text
du_full = solution of K * du = r
for alpha = 1, 0.5, 0.25, ...:
    trial u = u_iter + alpha * du_full
    evaluate residual/error at trial u
    accept alpha if merit decreases sufficiently
```

This targets the exact missing layer in the current solver.

## New input keywords

```text
nonlinear_line_search off|backtracking
nonlinear_line_search_merit residual|energy|mixed
nonlinear_line_search_reduction 0.5
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_armijo 1e-4
nonlinear_line_search_accept_any_decrease 1
nonlinear_line_search_freeze_material 1
nonlinear_line_search_cutback_on_fail 1
nonlinear_line_search_allow_alpha_one 1
```

Suggested defaults for first implementation:

```text
nonlinear_line_search backtracking
nonlinear_line_search_merit mixed
nonlinear_line_search_reduction 0.5
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_armijo 1e-4
nonlinear_line_search_accept_any_decrease 1
nonlinear_line_search_freeze_material 1
nonlinear_line_search_cutback_on_fail 1
```

## Merit function

Recommended first merit:

```cpp
phi = max(
    residual_active ? resErr / maxResErr : 0.0,
    energy_active   ? eneErr / maxEneErr : 0.0
);
```

Do not use raw `disErr` as the line-search merit. Any smaller alpha reduces the increment norm and can fake improvement.

Accepted if either:

```text
phi(alpha) < phi(0)
```

or Armijo-like:

```text
phi(alpha) <= (1 - c * alpha) * phi(0)
```

where `c = nonlinear_line_search_armijo`, typically `1e-4`.

For the first implementation, accept any actual decrease. It is more robust for inconsistent tangents and noisy residuals:

```cpp
bool ok = phiTrial < phiBefore;
```

Then later add Armijo as a stricter option.

## Critical state rule

A rejected trial alpha must not change the accepted solver state.

For every trial alpha, restore:

```text
trial_r
full_ddr
ddr
f_int
f_ext
residuals
f
resErr/disErr/eneErr
any material trial/status variables touched during stress evaluation
```

The difficult part is material state. Start with frozen trial residuals and reset material statuses between alpha trials. Then audit whether `frozen=true` is truly side-effect safe.

## Implementation architecture

### Step 1 — Add option fields

In `SteadyStateNonLinearSolver`:

```cpp
enum class NonlinearLineSearchType {
    Off,
    Backtracking
};

enum class NonlinearMeritType {
    Residual,
    Energy,
    Mixed
};

NonlinearLineSearchType nonlinearLineSearch = NonlinearLineSearchType::Off;
NonlinearMeritType nonlinearLineSearchMerit = NonlinearMeritType::Mixed;
double nonlinearLineSearchReduction = 0.5;
double nonlinearLineSearchMinAlpha = 0.03125;
unsigned nonlinearLineSearchMaxTrials = 6;
double nonlinearLineSearchArmijo = 1e-4;
bool nonlinearLineSearchAcceptAnyDecrease = true;
bool nonlinearLineSearchFreezeMaterial = true;
bool nonlinearLineSearchCutbackOnFail = true;
```

### Step 2 — Add a solver-state snapshot

In `solver_implicit.h`:

```cpp
struct NonlinearTrialState {
    Vector trial_r;
    Vector ddr;
    Vector full_ddr;
    Vector f_int;
    Vector f_ext;
    Vector residuals;
    Vector f;
    double resErr = 0.0;
    double disErr = 0.0;
    double eneErr = 0.0;
    double W_int = 0.0;
    double W_ext = 0.0;
};
```

Add helpers:

```cpp
void captureNonlinearTrialState(NonlinearTrialState& s) const;
void restoreNonlinearTrialState(const NonlinearTrialState& s);
```

Implementation outline:

```cpp
void SteadyStateNonLinearSolver::captureNonlinearTrialState(NonlinearTrialState& s) const {
    s.trial_r = trial_r;
    s.ddr = ddr;
    s.full_ddr = full_ddr;
    s.f_int = f_int;
    s.f_ext = f_ext;
    s.residuals = residuals;
    s.f = f;
    s.resErr = resErr;
    s.disErr = disErr;
    s.eneErr = eneErr;
    s.W_int = W_int;
    s.W_ext = W_ext;
}

void SteadyStateNonLinearSolver::restoreNonlinearTrialState(const NonlinearTrialState& s) {
    trial_r = s.trial_r;
    ddr = s.ddr;
    full_ddr = s.full_ddr;
    f_int = s.f_int;
    f_ext = s.f_ext;
    residuals = s.residuals;
    f = s.f;
    resErr = s.resErr;
    disErr = s.disErr;
    eneErr = s.eneErr;
    W_int = s.W_int;
    W_ext = s.W_ext;
}
```

Use the actual vector copy semantics and naming in the repository.

### Step 3 — Add a controlled trial-increment method

Avoid calling `updateFieldVariables()` blindly in a loop. Create a method that constructs a trial field from a fixed base state:

```cpp
bool SteadyStateNonLinearSolver::evaluateTrialIncrement(
    const Vector& baseTrialR,
    const Vector& duReducedFullNewton,
    double alpha,
    TrialResult& out)
{
    ddr = duReducedFullNewton;
    ddr.times(alpha);

    nodes->giveFullDoFArray(ddr, full_ddr);

    // Use baseTrialR, not the already-mutated trial_r.
    trial_r = baseTrialR;
    nodes->updateFullDoFsByDependenciesOnConjugates(full_ddr, trial_r, f_ext);
    trial_r.add(full_ddr);

    if (nonlinearLineSearchFreezeMaterial) {
        computeForcesAtIntegrationTime(true);
    } else {
        computeForcesAtIntegrationTime(false);
    }

    evaluateErrors();

    out.alpha = alpha;
    out.resErr = resErr;
    out.disErr = disErr;
    out.eneErr = eneErr;
    out.merit = currentNonlinearMerit();
    out.converged = convergenceCriteriaSatisfied();
    return finite(out.merit);
}
```

Important: `nodes->updateFullDoFsByDependenciesOnConjugates` may depend on `f_ext`; if indirect control or dependent DOFs use external forces, verify that `f_ext` is consistent before trial evaluation.

### Step 4 — Add `performBacktrackingLineSearch`

```cpp
bool SteadyStateNonLinearSolver::performBacktrackingLineSearch(
    const Vector& duNewton,
    const Vector& baseTrialR,
    double phiBefore,
    TrialResult& accepted)
{
    NonlinearTrialState baseState;
    captureNonlinearTrialState(baseState);

    double alpha = 1.0;
    bool haveBest = false;
    TrialResult best;

    lastLineSearchTrials = 0;

    for (unsigned i = 0; i < nonlinearLineSearchMaxTrials; ++i) {
        if (alpha < nonlinearLineSearchMinAlpha) break;

        restoreNonlinearTrialState(baseState);
        elems->resetMaterialStatuses();  // conservative; audit cost and semantics

        TrialResult trial;
        bool ok = evaluateTrialIncrement(baseTrialR, duNewton, alpha, trial);
        ++lastLineSearchTrials;

        if (ok && (!haveBest || trial.merit < best.merit)) {
            best = trial;
            haveBest = true;
        }

        if (ok && trial.converged) {
            accepted = trial;
            lastLineSearchAlpha = alpha;
            return true;
        }

        bool sufficientDecrease = false;
        if (ok) {
            if (nonlinearLineSearchAcceptAnyDecrease) {
                sufficientDecrease = trial.merit < phiBefore;
            } else {
                sufficientDecrease = trial.merit <= (1.0 - nonlinearLineSearchArmijo * alpha) * phiBefore;
            }
        }

        if (sufficientDecrease) {
            accepted = trial;
            lastLineSearchAlpha = alpha;
            return true;
        }

        alpha *= nonlinearLineSearchReduction;
    }

    if (haveBest && best.merit < phiBefore) {
        restoreNonlinearTrialState(baseState);
        elems->resetMaterialStatuses();
        evaluateTrialIncrement(baseTrialR, duNewton, best.alpha, accepted);
        lastLineSearchAlpha = best.alpha;
        return true;
    }

    restoreNonlinearTrialState(baseState);
    elems->resetMaterialStatuses();
    nonlinearCutbackReason = "line_search_failed";
    return false;
}
```

### Step 5 — Integrate into `solve()`

Current logic is conceptually:

```cpp
linalgsolver->solve(ddr, f);
updateFieldVariables();
computeForcesAtIntegrationTime(true);
evaluateErrors();
```

Change to:

```cpp
linalgsolver->solve(ddr, f);

Vector duNewton = ddr;
Vector baseTrialR = trial_r;
double phiBefore = currentNonlinearMeritBeforeStep();

bool acceptedIncrement = true;

if (nonlinearLineSearch == NonlinearLineSearchType::Backtracking) {
    TrialResult accepted;
    acceptedIncrement = performBacktrackingLineSearch(duNewton, baseTrialR, phiBefore, accepted);
} else if (nonlinearDampingType != Off) {
    applyDampedIncrement(duNewton, baseTrialR);
    computeForcesAtIntegrationTime(true);
    evaluateErrors();
} else {
    ddr = duNewton;
    updateFieldVariables();
    computeForcesAtIntegrationTime(true);
    evaluateErrors();
}

if (!acceptedIncrement) {
    if (nonlinearLineSearchCutbackOnFail && dt > dtmin * 1.00001) {
        restart_now = true;
        converged = false;
        break;
    }
    // Fallback: apply minimum alpha only if explicitly allowed.
}
```

### Step 6 — Compute `phiBefore`

At the start of each nonlinear iteration, before solving `K * du = r`, ensure you have a merit for the current state. Options:

1. Use the previous iteration's `resErr/eneErr`.
2. On the first iteration of a load step, compute residual and call `evaluateErrors()` before entering the first solve.
3. Add a dedicated residual norm function that does not depend on `ddr`.

Recommended first implementation:

```cpp
if (it == 0) {
    computeForcesAtIntegrationTime(true);
    evaluateErrors();
}
double phiBefore = currentNonlinearMerit();
```

Be careful: `evaluateErrors()` may use `ddr/full_ddr` for displacement and energy terms. For merit before the Newton increment, use residual/energy only or write a dedicated function:

```cpp
double currentResidualEnergyMeritOnly() const;
```

## Checklist

- [ ] Add parser options and defaults.
- [ ] Add `NonlinearTrialState` snapshot/restore.
- [ ] Add trial-increment evaluation that uses a fixed `baseTrialR`.
- [ ] Do not call `updateFieldVariables()` repeatedly without restoring `trial_r`.
- [ ] Use residual/energy merit, not displacement increment alone.
- [ ] Reset material statuses between rejected trials, at least for the first implementation.
- [ ] Log selected alpha and number of alpha trials.
- [ ] Failed line search triggers load-step cutback when `dt > dtmin`.
- [ ] Keep `limit_tolerance` disabled while testing.
- [ ] Add unit test where alpha `1.0` fails but alpha `0.5` passes.

## Test cases

### Synthetic nonlinear spring

Residual:

```text
R(u) = P - k u + a u^3
```

Pick a load where full Newton overshoots. Verify:

```text
alpha=1 rejected
alpha<1 accepted
residual merit decreases
```

### Material-status mutation test

Create a mock material status with counters:

```text
trial_evaluations
committed_updates
resets
```

Run line search with rejected alphas. Verify:

```text
rejected alphas do not increment committed_updates
final committed state corresponds only to accepted alpha after step success
```

### TS-N65 regression

Metrics to collect:

```text
mean nonlinear iterations
median nonlinear iterations
max nonlinear iterations
number of cutbacks
number of accepted fallback steps
min alpha histogram
wall-clock time
load-displacement curve
```

## Risks

- If `frozen=true` still mutates trial material variables, rejected alpha trials can pollute state.
- If the tangent is very poor, residual may not decrease even for small alpha.
- If load control is unstable in snapback, line search can still fail correctly; then use displacement/arc-length control.

## Acceptance criteria

- No rejected alpha trial changes committed material history.
- For difficult steps, accepted alpha is often below 1 and residual/energy merit decreases.
- TS-N65 average nonlinear iterations drop materially compared with baseline.
- Fallback acceptance through `limit_tolerance` is not needed.

---

# Option 4 — Material-state snapshot/rollback for trial evaluations

## Promise

This is the robustness foundation for any serious line search, trust region, or arc-length implementation in a path-dependent material solver.

The current code has `updateMaterialStatuses()` and `resetMaterialStatuses()`, but line search needs a stronger guarantee:

```text
Any rejected trial evaluation must leave the solver exactly as it was before the trial.
```

For softening damage, plasticity, crack opening, and cumulative plastic strain, trial evaluations may update temporary variables even when not committed. You must audit and enforce rollback semantics.

## Two rollback levels

### Level A — Reset to last committed step state

Use existing:

```cpp
elems->resetMaterialStatuses();
```

Pros:

- Low implementation cost.
- Likely sufficient for a first frozen line search.

Cons:

- May not restore arbitrary within-step trial states.
- If the nonlinear iteration relies on evolving uncommitted state across iterations, resetting to last committed state can change behavior.

### Level B — Full snapshot/restore of material statuses

Add clone/restore methods for every material status object.

Pros:

- Correct for arbitrary trial evaluations.
- Supports line search, trust region, arc-length, and constitutive substepping.

Cons:

- Requires touching material/status hierarchy.
- Memory cost can be high.
- Must include all hidden history and temporary variables.

## Recommended path

Implement Level A for the first line-search patch, but design the API so Level B can replace it without changing solver logic.

## Proposed API

In `element_container.h`:

```cpp
class MaterialStatusSnapshot;

std::unique_ptr<MaterialStatusSnapshot> createMaterialStatusSnapshot() const;
void restoreMaterialStatusSnapshot(const MaterialStatusSnapshot& snapshot);
```

In material status base class:

```cpp
class MaterialStatus {
public:
    virtual std::unique_ptr<MaterialStatus> clone() const = 0;
    virtual void restoreFrom(const MaterialStatus& other) = 0;
};
```

For each concrete status:

```cpp
std::unique_ptr<MaterialStatus> CSLMaterialStatus::clone() const override {
    return std::make_unique<CSLMaterialStatus>(*this);
}

void CSLMaterialStatus::restoreFrom(const MaterialStatus& other) override {
    const auto& src = dynamic_cast<const CSLMaterialStatus&>(other);
    *this = src;
}
```

If dynamic cast overhead is a concern, store typed snapshots per material class.

## Snapshot container design

```cpp
class MaterialStatusSnapshot {
public:
    std::vector<std::vector<std::unique_ptr<MaterialStatus>>> statusesByElement;
};
```

The actual hierarchy may store statuses in elements, integration points, or material containers. Mirror the existing ownership structure.

## Trial evaluation pattern with snapshot

```cpp
auto solverState = captureNonlinearTrialState();
auto materialState = elems->createMaterialStatusSnapshot();

for (double alpha : alphaList) {
    restoreNonlinearTrialState(solverState);
    elems->restoreMaterialStatusSnapshot(*materialState);

    evaluateTrialIncrement(baseTrialR, duNewton, alpha, trial);

    if (accept(trial)) {
        // Leave solver/material state at accepted trial.
        return true;
    }
}

restoreNonlinearTrialState(solverState);
elems->restoreMaterialStatusSnapshot(*materialState);
return false;
```

## Audit checklist for material statuses

For every material used in TS-N65:

- [ ] List all committed history variables.
- [ ] List all trial/temporary variables.
- [ ] List all cached stiffness/tangent variables.
- [ ] List all crack/damage/plasticity variables.
- [ ] List random/eigenstrain/internal variable fields, if any.
- [ ] Confirm copy constructor copies every relevant field.
- [ ] Confirm `resetMaterialStatuses()` restores exactly last committed state.
- [ ] Confirm `updateMaterialStatuses()` commits only after a converged load step.
- [ ] Confirm `evaluateStresses(frozen=true)` does not commit history.
- [ ] Confirm rejected trial evaluations do not affect next alpha.

## Unit-test checklist

- [ ] Snapshot/restore round-trip leaves all scalar, vector, and tensor fields equal.
- [ ] Snapshot/restore works after multiple trial evaluations.
- [ ] Snapshot restore does not leak memory.
- [ ] Snapshot restore handles polymorphic material status classes.
- [ ] Parallel/MPI behavior is deterministic if the solver is run in parallel.

## Acceptance criteria

- Rejected line-search alphas are perfectly side-effect free.
- Trust-region and arc-length code can reuse the same trial-evaluation API.
- Material-state rollback tests pass for every material in TS-N65.

---

# Option 5 — Stiffness update and tangent consistency improvements

## Promise

If `K` is stale or inconsistent with the residual, Newton iterations degrade into fixed-point iterations. The note says TS-N65 used `stiffness_matrix_iter_update 10`. That means the solver can reuse the same matrix for many nonlinear iterations while the material is softening/damaging.

## Immediate experiments

Use existing controls:

```text
stiffness_matrix_iter_update 1
stiffness_matrix_iter_update 3
stiffness_matrix_iter_update 10
```

Also test:

```text
first_iteration_stiff_matrix_type elastic
first_iteration_stiff_matrix_type secant
stiff_matrix_type secant
stiff_matrix_type tangent
stiff_matrix_type consistent
```

Only trust `tangent`/`consistent` where the material implementation actually provides a consistent tangent. The supplied note indicates that for `CSLMaterial`, the `tangent` branch may return unloading/degraded stiffness rather than a true consistent tangent.

## Implementation improvements

### 5.1 Adaptive stiffness rebuild

Add criteria to rebuild `K` when the nonlinear solver stalls or the accepted alpha is too small.

New keywords:

```text
nonlinear_adaptive_matrix_update 0/1
nonlinear_rebuild_on_small_alpha 0.5
nonlinear_rebuild_on_stagnation 1
nonlinear_rebuild_on_residual_growth 1
```

Logic:

```cpp
bool forceRebuildK = false;

if (lastLineSearchAlpha < nonlinearRebuildOnSmallAlpha) {
    forceRebuildK = true;
}
if (meritGrew || stagnated) {
    forceRebuildK = true;
}

if (forceRebuildK) {
    updateSystemMatrices();
    factorizeLinearSystem();
}
```

This is useful when matrix assembly is costly but stale stiffness causes many iterations.

### 5.2 Tangent consistency audit

For each material tangent implementation:

- [ ] Compare analytical tangent with finite-difference derivative of internal force.
- [ ] Check tangent symmetry expectations.
- [ ] Check unloading/reloading branches.
- [ ] Check softening damage branches.
- [ ] Check reinforcement plasticity branches.
- [ ] Check whether tangent includes history-variable algorithmic derivatives.

Finite-difference test:

```text
Given converged state u and perturbation p:
K p ≈ [f_int(u + eps p) - f_int(u)] / eps
```

Run for:

```text
eps = 1e-4, 1e-5, 1e-6, 1e-7
```

Measure:

```text
relative_tangent_error = ||Kp - fd|| / max(||fd||, tiny)
```

### 5.3 Consistent tangent implementation

A true consistent tangent should linearize the same stress update used to compute residual. For path-dependent materials this includes the local return-mapping/damage update algorithm, not only elastic or secant stiffness.

Implementation checklist:

- [ ] Identify exact constitutive update called during residual evaluation.
- [ ] Derive algorithmic derivative of stress with respect to strain increment.
- [ ] Include damage evolution derivative where smooth.
- [ ] Add branch-safe derivatives for loading/unloading.
- [ ] Regularize discontinuous switches where necessary.
- [ ] Expose matrix type only when implementation is valid.
- [ ] Add finite-difference tangent tests.

## Risks

- A bad “consistent” tangent can be worse than secant stiffness.
- Tangent changes can alter physical response if they accidentally change stress update logic.
- More frequent matrix updates increase per-iteration cost; total runtime may still improve.

## Acceptance criteria

- `stiffness_matrix_iter_update 1` materially reduces iteration counts for difficult steps.
- Adaptive matrix rebuild gives near-`1` behavior only when needed.
- Finite-difference tangent error is acceptable in elastic/plastic/damage regimes.

---

# Option 6 — Trust-region and Levenberg-Marquardt-style regularization

## Promise

Trust-region methods help when the Newton direction is not reliable. In softening damage, the tangent can become indefinite or a poor local model of the residual. A trust region limits the step size based on model quality.

This is more complex than line search but less intrusive than arc-length.

## Two practical variants

### Variant A — Step-norm trust region

Compute Newton step `du`, then cap its norm:

```text
if ||du|| > Delta:
    du = du * Delta / ||du||
```

Accept/reject based on actual merit reduction. Update `Delta`:

```text
if merit reduction is good: increase Delta
if bad: decrease Delta
```

This is easy and resembles adaptive damping.

### Variant B — Levenberg-Marquardt-style regularized Newton

Solve:

```text
(K + mu D) du = r
```

where `D` is a positive diagonal scaling, for example:

```text
D = abs(diag(K)) + diagonal_floor
```

If the trial step improves merit, reduce `mu`. If not, increase `mu` and retry.

## New input keywords

```text
nonlinear_trust_region off|step_norm|lm
nonlinear_trust_radius_initial 1.0
nonlinear_trust_radius_min 1e-10
nonlinear_trust_radius_max 1e10
nonlinear_trust_shrink 0.5
nonlinear_trust_expand 2.0
nonlinear_trust_eta_accept 0.1
nonlinear_lm_mu_initial 1e-8
nonlinear_lm_mu_min 1e-14
nonlinear_lm_mu_max 1e8
nonlinear_lm_mu_shrink 0.1
nonlinear_lm_mu_expand 10.0
```

## Implementation guide for step-norm trust region

### Step 1 — Define norm

Use a displacement norm consistent with `disErr`:

```cpp
double normDu = computeReducedIncrementNorm(ddr);
```

Prefer a scaled norm:

```text
||du||_scaled = sqrt(sum_i (du_i / scale_i)^2)
```

where `scale_i` can be based on current displacement magnitude or DOF type.

### Step 2 — Cap step

```cpp
if (normDu > trustRadius) {
    double alpha = trustRadius / normDu;
    ddr.times(alpha);
    lastLineSearchAlpha = alpha;
}
```

Then use the same trial-evaluation API as line search.

### Step 3 — Update radius

```cpp
double reduction = phiBefore - phiAfter;

if (reduction > 0.0) {
    trustRadius = std::min(trustRadiusMax, trustExpand * trustRadius);
    accept;
} else {
    trustRadius = std::max(trustRadiusMin, trustShrink * trustRadius);
    reject or retry;
}
```

## Implementation guide for LM regularization

### Step 1 — Build a regularized matrix

After assembling `K`, construct:

```cpp
Kreg = K;
for each diagonal entry i:
    Kreg(i,i) += mu * max(abs(K(i,i)), diagonal_floor);
```

The actual sparse matrix class determines how to do this. If diagonal modification is not simple, add a matrix method:

```cpp
addScaledAbsoluteDiagonal(double mu, double floor);
```

### Step 2 — Solve and trial-evaluate

```cpp
for (attempt = 0; attempt < maxLmAttempts; ++attempt) {
    factorize(Kreg);
    solve(du, r);
    evaluateTrialIncrement(...);

    if (merit improves) {
        mu = max(muMin, muShrink * mu);
        accept;
        break;
    }

    mu = min(muMax, muExpand * mu);
}
```

### Step 3 — Interact with line search

Use either:

```text
LM regularization first, then line search on resulting direction
```

or:

```text
line search first, LM only when line search fails
```

Recommended first combination:

1. Standard Newton direction.
2. Backtracking line search.
3. If line search fails, cut back.
4. Add LM later as an optional fallback before cutback.

## Checklist

- [ ] Reuse the line-search trial-evaluation and rollback API.
- [ ] Do not overwrite the main `K` permanently unless intended.
- [ ] Log `trustRadius` or `mu`.
- [ ] Verify diagonal regularization does not break matrix format or solver assumptions.
- [ ] Test on indefinite tangent cases.

## Risks

- LM-style diagonal shifts can distort mechanics if overused.
- Sparse matrix modification may be expensive.
- Trust-region methods need careful scaling across DOF types.

## Acceptance criteria

- Trust-region mode recovers cases where line search repeatedly fails.
- `mu` or `trustRadius` adapts predictably.
- Accepted steps still satisfy original residual/energy/displacement tolerances.

---

# Option 7 — Nonlinear acceleration and quasi-Newton methods

## Promise

These methods can reduce matrix assembly cost or improve fixed-point convergence, but they are less directly targeted than line search and cutback.

Possible methods:

```text
Broyden update
BFGS/L-BFGS inverse Jacobian approximation
Anderson acceleration
nonlinear GMRES
```

## When to try

Try after line search if:

- Matrix assembly/factorization is expensive.
- Residual behavior is smooth enough.
- The tangent is not reliable but residual evaluations are cheap.
- You want fewer full matrix rebuilds.

Do not use these as the first fix for severe damage snapback or material-state rollback problems.

## Practical implementation: Anderson acceleration

For an iteration map:

```text
u_{k+1} = G(u_k) = u_k + alpha du_k
```

store recent displacement and residual differences, then compute a linear combination that minimizes residual norm.

New options:

```text
nonlinear_acceleration off|anderson
nonlinear_anderson_depth 5
nonlinear_anderson_start_iteration 3
nonlinear_anderson_regularization 1e-10
nonlinear_anderson_restart_on_cutback 1
```

Implementation outline:

```cpp
store u_k, r_k for last m iterations
solve small least-squares problem for coefficients
construct accelerated trial u
trial-evaluate with rollback
accept only if merit decreases
otherwise fall back to line-search Newton step
```

## Practical implementation: Broyden/L-BFGS

Maintain approximate inverse Jacobian action:

```text
du = H_k r
```

Update using:

```text
s_k = u_{k+1} - u_k
y_k = r_{k+1} - r_k
```

Use limited memory to avoid large dense matrices.

## Checklist

- [ ] Implement only after robust trial evaluation exists.
- [ ] Always accept/reject accelerated trial by merit decrease.
- [ ] Clear acceleration history after cutback, restart, material status reset, or stiffness rebuild.
- [ ] Bound accelerated step norm.
- [ ] Log when acceleration is used versus rejected.

## Risks

- Can be unstable for nonsmooth damage/plasticity switches.
- Can combine incompatible states if material history changes discontinuously.
- Requires careful handling of constraints/dependent DOFs.

## Acceptance criteria

- Reduces total matrix rebuilds or total iterations on smooth portions.
- Never worsens difficult softening steps because rejected acceleration falls back to line search.

---

# Option 8 — Improve indirect displacement control

## Promise

If load control is unstable in the post-peak regime, a better nonlinear step length may not be enough. The existing code has `indirect_control` / `indirect_displacement_control`, which changes the load multiplier to hit a prescribed measure. This can be used earlier and more systematically.

## When to use

Use indirect/displacement control when:

- Load-controlled steps fail even with line search and cutback.
- The response has snap-through or snapback.
- A meaningful control displacement/strain/crack-opening measure is known.
- You need to trace post-peak equilibrium rather than only reach stable load levels.

## Improvements to implement

### 8.1 Better control-variable selection

Add support for multiple control choices:

```text
single DOF displacement
relative displacement between two nodes
average displacement over node set
crack opening measure
reaction force target
square-sum measure already partly available
```

New keywords:

```text
indirect_control_type dof|relative_dof|node_set_average|square_sum|reaction
indirect_control_node <id>
indirect_control_dof x|y|z|rx|ry|rz
indirect_control_node_a <id>
indirect_control_node_b <id>
indirect_control_node_set <set_id>
indirect_control_target_increment <value>
```

### 8.2 Automatic switch from load control to indirect control

Add a trigger:

```text
nonlinear_auto_indirect_control 0/1
nonlinear_auto_indirect_after_cutbacks 3
nonlinear_auto_indirect_after_small_alpha 0.03125
```

Logic:

```cpp
if (consecutiveCutbacks >= threshold || repeatedLineSearchFailure) {
    switchToIndirectControl();
}
```

Only do this if the input defines a valid control variable.

### 8.3 Combine indirect control with line search

Indirect control chooses load multiplier correction. Line search still helps because the coupled update can overshoot.

Use:

```text
indirect/displacement control equation + line search on the combined increment
```

## Checklist

- [ ] Identify a physically meaningful displacement or crack-opening control for TS-N65.
- [ ] Add diagnostics: current load multiplier, control value, target control value.
- [ ] Keep line search active.
- [ ] Compare load-displacement curve with pure load control.
- [ ] Detect wrong control DOF sign/direction.

## Risks

- Bad control variable can hide global instability or trace the wrong path.
- Switching control modes mid-analysis can complicate reproducibility.
- Existing indirect control may not be general enough for all snapback cases.

## Acceptance criteria

- Cases failing under load control can progress under displacement/indirect control.
- The chosen control variable evolves monotonically enough for stable continuation.
- Load-displacement curve is physically interpretable.

---

# Option 9 — Arc-length / Riks continuation

## Promise

Arc-length is the correct continuation strategy when equilibrium has snap-through or snapback and neither load control nor a single displacement control is sufficient.

It augments the unknowns:

```text
u      = displacement vector
lambda = load factor
```

and solves:

```text
R(u, lambda) = lambda * f_ref - f_int(u) = 0
g(Delta u, Delta lambda) = 0
```

where `g` is an arc-length constraint, commonly:

```text
Delta u^T W Delta u + psi^2 Delta lambda^2 = Delta s^2
```

## New input keywords

```text
nonlinear_control load|indirect|arc_length
arc_length_radius_initial <value>
arc_length_radius_min <value>
arc_length_radius_max <value>
arc_length_psi <value>
arc_length_shrink 0.5
arc_length_expand 1.2
arc_length_target_iterations 8
arc_length_max_iterations 40
arc_length_sign_strategy previous_increment|positive_load|control_dof|monotone_lambda
arc_length_constraint spherical|gauge|cylindrical|updated_normal_plane
arc_length_auto_radius 0|1
arc_length_gauge_tolerance 1e-3
```

Current implementation note: `spherical` is the generic CP4 prototype.
`gauge` is a TS-N65/IDC continuation experiment that uses a displacement-only
`indirect_control` block as the corrector constraint and keeps defaults at
legacy load control unless explicitly enabled.

## Implementation plan

### Step 1 — Separate reference load from current load

You need a clean representation:

```text
f_ext(lambda) = f_ext_fixed + lambda * f_ext_ref
```

The current solver uses pseudo-time/load `load`. Arc-length requires treating load factor as an unknown, not merely a prescribed time value.

### Step 2 — Predictor step

Use previous converged tangent direction:

```text
K du_bar = f_ref
Delta lambda = sign * Delta s / sqrt(du_bar^T W du_bar + psi^2)
Delta u = Delta lambda * du_bar
```

### Step 3 — Corrector equations

At each corrector iteration solve a bordered system:

```text
[ K   -f_ref ] [du_corr     ] = [ R ]
[ c^T  c_lam ] [dlambda_corr]   [ -g]
```

Where `c` and `c_lam` are derivatives of the arc-length constraint.

### Step 4 — Reuse line search/trial evaluation

Arc-length corrector steps still need globalization:

```text
(u, lambda) trial = (u, lambda) + alpha * (du_corr, dlambda_corr)
```

Use the same merit and rollback rules.

### Step 5 — Adapt arc-length radius

```text
if iterations <= target: radius *= expand
if iterations > target or cutback: radius *= shrink
```

## Implementation checklist

- [ ] Refactor load application so `lambda` can be solved as an unknown.
- [ ] Add storage for previous converged increment direction.
- [ ] Add bordered linear system solver or Schur-complement solve.
- [ ] Add arc-length constraint evaluation and derivative.
- [ ] Add sign-selection logic to continue along the correct branch.
- [ ] Add load/displacement output that records `lambda` separately from pseudo-time.
- [ ] Reuse material rollback for all corrector trial evaluations.
- [ ] Start with small benchmark, not TS-N65.

## Schur-complement implementation

Avoid assembling a fully bordered matrix initially. Solve:

```text
K a = R
K b = f_ref
```

Then use constraint equation to solve for `dlambda`:

```text
du = a + dlambda * b
c^T du + c_lam dlambda = -g

dlambda = (-g - c^T a) / (c^T b + c_lam)
```

Then:

```text
du = a + dlambda * b
```

## Risks

- More intrusive than line search.
- Requires careful handling of load patterns and prescribed displacements.
- Branch switching can occur near bifurcations.
- Arc-length traces mathematical equilibria that may include unstable physical states.

## Acceptance criteria

- Reproduces standard snap-through/snapback benchmark curves.
- Continues TS-N65-like post-peak regimes where load control fails.
- Does not require `limit_tolerance` fallback.

---

# Option 10 — Material-level substepping

## Promise

Global load-step cutback reduces the whole problem step. Material substepping reduces only the local constitutive update increment at integration points. This helps when local damage/plasticity evolution overshoots even though the global Newton step is reasonable.

## When to use

Use if diagnostics show:

- Stress updates fail or produce NaN/Inf.
- Damage/plastic strain jumps abruptly in one global increment.
- Line search repeatedly accepts very small alpha due to local material instability.
- Local return mapping has convergence issues.

## New input keywords

```text
material_substepping 0/1
material_substepping_max 16
material_substepping_min_fraction 1e-4
material_substepping_error_tolerance 1e-3
material_substepping_damage_increment_limit 0.02
material_substepping_plastic_strain_limit 1e-4
```

## Implementation concept

Inside material stress update:

```text
strain_increment = epsilon_trial - epsilon_committed
split strain_increment into n local subincrements
for each subincrement:
    perform stress/damage/plasticity update
    if local criterion fails: increase n and retry
```

## Local criteria

Examples:

```text
abs(delta_damage) <= damage_increment_limit
norm(delta_plastic_strain) <= plastic_strain_limit
local residual converged
stress finite
energy dissipation nonnegative
```

## Implementation guide

### Step 1 — Add local trial/commit separation

The material update must distinguish:

```text
committed state at last converged global step
local substep state during stress evaluation
trial state returned to element
```

### Step 2 — Add adaptive local substep loop

```cpp
bool updateStressWithSubstepping(const Strain& epsTrial, MaterialStatus& status) {
    for (unsigned n = 1; n <= maxSubsteps; n *= 2) {
        auto local = status.committedClone();
        bool ok = true;

        for (unsigned i = 0; i < n; ++i) {
            Strain epsSub = interpolate(epsCommitted, epsTrial, double(i+1)/n);
            ok = updateOneSubstep(epsSub, local);
            if (!ok || violatesLocalLimits(local)) break;
        }

        if (ok) {
            status.setTrialFrom(local);
            return true;
        }
    }
    return false;
}
```

### Step 3 — Report local failures to global solver

If material substepping fails, return a status flag to the element/container/solver:

```text
material_update_failed
max_local_damage_increment_exceeded
local_return_mapping_failed
```

Then trigger global cutback.

## Checklist

- [ ] Implement on one material first.
- [ ] Track number of local substeps per integration point.
- [ ] Log worst local substep count per global iteration.
- [ ] Add local failure flags that global solver can read.
- [ ] Verify energy dissipation sign and stress continuity.

## Risks

- More expensive stress evaluation.
- Can mask underlying ill-posed damage localization.
- Requires constitutive-law expertise.

## Acceptance criteria

- Local NaNs/discontinuous jumps disappear or reduce.
- Global line search accepts larger alphas.
- Load-displacement curve remains physically consistent.

---

# Option 11 — Viscous or rate regularization

## Promise

Viscous regularization can stabilize softening damage by adding a rate-dependent term. It can make a quasi-static softening problem numerically tractable, especially when abrupt damage localization creates near-discontinuous response.

This is a modeling change, not a pure nonlinear-solver enhancement.

## New input keywords

```text
damage_viscosity 0.0
plasticity_viscosity 0.0
viscous_regularization_type none|damage|plasticity|both
viscous_energy_output 1
```

## Implementation idea

For damage variable `d`, instead of immediate inviscid update:

```text
d = d_eq(kappa)
```

use:

```text
eta * d_dot + d = d_eq(kappa)
```

Backward Euler:

```text
d_{n+1} = (eta / dt * d_n + d_eq) / (eta / dt + 1)
```

## Checklist

- [ ] Implement only with clear documentation that it changes the material model.
- [ ] Output viscous dissipation separately.
- [ ] Run viscosity sensitivity study.
- [ ] Ensure results converge as viscosity tends to zero, if physically expected.
- [ ] Do not use viscosity merely to force convergence without reporting it.

## Risks

- Alters peak load and post-peak response.
- Depends on pseudo-time step.
- Can hide mesh-dependence or localization issues.

## Acceptance criteria

- Convergence improves without unacceptable change in physical response.
- Sensitivity to viscosity is quantified.
- Viscous energy is small relative to fracture/damage energy when claiming quasi-static behavior.

---

# Option 12 — Nonlocal, gradient, or crack-band localization regularization

## Promise

If the root problem is mesh-sensitive localization, no nonlinear solver can fully fix it. The mathematical model may be ill-posed after softening. Regularization introduces a characteristic length or fracture-energy scaling.

Possible approaches:

```text
crack-band fracture-energy regularization
nonlocal damage averaging
gradient damage
phase-field fracture
LDPM-specific facet/element size energy scaling
```

## Implementation options

### Crack-band style regularization

Scale softening law by element/facet characteristic length so dissipated energy is mesh objective.

Checklist:

- [ ] Define characteristic length for each LDPM facet/element.
- [ ] Verify area/volume weighting.
- [ ] Integrate stress-crack opening relation to target fracture energy.
- [ ] Test mesh refinement.

### Nonlocal damage averaging

Compute nonlocal equivalent strain/history variable:

```text
kappa_bar(x) = integral w(|x-y|) kappa(y) dy / integral w(|x-y|) dy
```

Checklist:

- [ ] Build neighbor search.
- [ ] Choose length scale.
- [ ] Handle boundaries.
- [ ] Add consistent tangent or accept secant-like iteration with line search.

### Gradient damage

Introduce additional field variable and PDE-like regularization.

Checklist:

- [ ] Add new DOFs or internal solve.
- [ ] Assemble coupled stiffness blocks.
- [ ] Add boundary conditions for gradient field.
- [ ] Validate on standard localization benchmark.

## Risks

- Significant model development.
- Requires calibration.
- Not a local patch to nonlinear solver.

## Acceptance criteria

- Mesh-refinement response is stable.
- Damage zone width and dissipated energy are physically controlled.
- Nonlinear convergence improves as a consequence of better-posed equations.

---

# Option 13 — Linear solver/preconditioner tuning

## Promise

The supplied note says the immediate failure is not inability to solve individual linear systems. Therefore, linear solver tuning is lower priority for this problem. Still, inaccurate linear solves can degrade nonlinear convergence.

## What to check

- Linear residual tolerance.
- Number of linear iterations.
- Whether `DeflatedFGMRES` actually converges to requested tolerance.
- Whether preconditioner changes cause different nonlinear behavior.
- Whether stiffness matrix is singular/ill-conditioned in softening regimes.

## Practical additions

Add to logs:

```text
linear_solver_type
linear_iterations
linear_initial_residual
linear_final_residual
linear_converged
preconditioner_rebuild
```

## Checklist

- [ ] Verify linear solve accuracy on slow nonlinear steps.
- [ ] Tighten linear tolerance and compare nonlinear iteration count.
- [ ] Loosen linear tolerance and see whether nonlinear convergence worsens.
- [ ] Do not attribute nonlinear snapback to the linear solver without evidence.

## Acceptance criteria

- Linear solves are accurate enough that nonlinear globalization behavior is meaningful.
- Linear solver changes are not the primary mechanism of improvement unless measured.

---

# Option 14 — PETSc SNES or external nonlinear-solver integration

## Promise

Mature nonlinear libraries provide line search, trust region, quasi-Newton, nonlinear GMRES, convergence tests, and monitors. This is attractive long-term but intrusive.

## Why not first

The code still needs clean residual evaluation, trial state rollback, and material-status semantics. External solvers cannot fix unsafe residual callbacks.

## Integration plan

### Step 1 — Define residual callback

```cpp
void computeResidual(const Vector& u, Vector& R);
```

Requirements:

- Does not commit material state.
- Can evaluate at arbitrary trial `u`.
- Restores previous state after rejected trial.
- Applies boundary/dependent DOF logic consistently.

### Step 2 — Define Jacobian callback

```cpp
void computeJacobian(const Vector& u, SparseMatrix& J);
```

or use matrix-free approximation initially.

### Step 3 — Wrap OAS vectors/matrices

Either:

```text
copy OAS vectors into PETSc Vec/Mat every evaluation
```

or:

```text
refactor OAS vector/matrix storage to expose PETSc-compatible buffers
```

### Step 4 — Map convergence criteria

Map OAS criteria:

```text
resErr
disErr
eneErr
```

to SNES convergence tests and monitors.

## Checklist

- [ ] Implement local line search first to clarify callback semantics.
- [ ] Ensure material rollback is correct.
- [ ] Start with a small model and compare residuals exactly.
- [ ] Add SNES monitor output to existing solver logs.
- [ ] Preserve existing solver path as fallback.

## Acceptance criteria

- SNES residual callback is side-effect free.
- SNES line search/trust-region methods reproduce or improve local implementation.
- Integration does not require changing material physics.

---

# Recommended first patch: line search + stagnation cutback

## Patch scope

Add these features only:

```text
nonlinear_line_search off|backtracking
nonlinear_line_search_min_alpha
nonlinear_line_search_reduction
nonlinear_line_search_max_trials
nonlinear_line_search_merit residual|energy|mixed
nonlinear_line_search_armijo
nonlinear_line_search_accept_any_decrease
nonlinear_line_search_freeze_material
nonlinear_line_search_cutback_on_fail
nonlinear_stagnation_cutback
nonlinear_stagnation_iterations
nonlinear_stagnation_ratio
nonlinear_growth_cutback
```

Do not change material laws in the first patch.

## First-patch default behavior

```text
nonlinear_line_search off
nonlinear_stagnation_cutback 0
```

Keep defaults backward-compatible. Enable in TS-N65 input explicitly:

```text
nonlinear_line_search backtracking
nonlinear_line_search_merit mixed
nonlinear_line_search_reduction 0.5
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_accept_any_decrease 1
nonlinear_line_search_freeze_material 1
nonlinear_line_search_cutback_on_fail 1

nonlinear_stagnation_cutback 1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25

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

## First-patch pseudocode

```cpp
while (!converged && it < maxIt) {
    bool rebuildK = shouldUpdateStiffness(it, cumul_it, ...);
    if (rebuildK) {
        updateSystemMatrices();
        factorizeLinearSystem();
    }

    double phiBefore = residualEnergyMerit();

    linalgsolver->solve(ddr, f);

    Vector duNewton = ddr;
    Vector baseTrialR = trial_r;

    bool acceptedIncrement = true;

    if (nonlinearLineSearch == Backtracking) {
        TrialResult accepted;
        acceptedIncrement = performBacktrackingLineSearch(
            duNewton,
            baseTrialR,
            phiBefore,
            accepted
        );
    } else {
        updateFieldVariables();
        computeForcesAtIntegrationTime(true);
        evaluateErrors();
    }

    if (!acceptedIncrement) {
        converged = false;
        restart_now = true;
        nonlinearCutbackReason = "line_search_failed";
        break;
    }

    double phiAfter = residualEnergyMerit();
    if (shouldCutbackForStagnation(phiAfter)) {
        converged = false;
        restart_now = true;
        break;
    }

    converged = convergenceCriteriaSatisfied() && it >= minIt;
    ++it;
    ++cumul_it;
}
```

## First-patch code checklist

### Header changes

- [ ] Add line-search enum.
- [ ] Add merit enum.
- [ ] Add option fields.
- [ ] Add state snapshot struct.
- [ ] Add trial result struct.
- [ ] Add helper method declarations.

### Parser changes

- [ ] Parse all new keywords.
- [ ] Validate ranges: reduction in `(0,1)`, min alpha positive, max trials >= 1.
- [ ] Print warnings for invalid combinations.

### Solve-loop changes

- [ ] Store `duNewton` after linear solve.
- [ ] Store `baseTrialR` before any trial update.
- [ ] Call line-search helper before full-step update.
- [ ] On line-search failure, trigger existing cutback path.
- [ ] Reset line-search/stagnation state after step restart.

### Trial evaluation changes

- [ ] Apply trial increment relative to fixed `baseTrialR`.
- [ ] Evaluate forces with `frozen=true` initially.
- [ ] Evaluate errors.
- [ ] Restore solver vectors after rejection.
- [ ] Reset material statuses after rejection.

### Logging changes

- [ ] Log alpha and trial count.
- [ ] Log cutback reason.
- [ ] Log merit before/after.
- [ ] Log fallback acceptance separately from true convergence.

### Tests

- [ ] Unit test for alpha backtracking.
- [ ] Unit test for no state mutation on rejected alphas.
- [ ] Regression test with line search off: same result as before.
- [ ] Regression test with line search on: fewer stalls/cutbacks.

---

# Detailed solver-state handling for line search

## Why this matters

The numerical algorithm is easy. The state management is hard. The solver has global vectors, nodal trial fields, internal/external forces, residuals, energy variables, and material status variables. A trial alpha must be reversible.

## State to save before alpha loop

At minimum:

```text
trial_r
r if distinct from trial_r
full_ddr
ddr
f
f_int
f_ext
residuals
resErr
disErr
eneErr
W_int
W_ext
load if indirect control changes it
any indirect-control internal state
material statuses or at least reset boundary
```

## State to recompute for accepted alpha

After accepting alpha, leave solver state as if this code had run exactly once:

```cpp
ddr = alpha * duNewton;
updateFieldVariables();
computeForcesAtIntegrationTime(true);
evaluateErrors();
```

But because `updateFieldVariables()` mutates `trial_r`, the accepted state should be produced by a controlled helper, not by accumulating prior rejected trials.

## Trial helper invariant

For every alpha:

```text
trial_r_trial = baseTrialR + full(alpha * duNewton)
```

not:

```text
trial_r_trial = previousTrialR + full(alpha * duNewton)
```

## Material-state invariant

For every rejected alpha:

```text
committed material state after rejection = committed material state before alpha loop
trial material state after rejection = trial material state before alpha loop
```

If Level B snapshot is not implemented, at least enforce:

```cpp
elems->resetMaterialStatuses();
```

between alpha trials and before returning failure.

---

# Merit, convergence, and fallback policy

## Recommended merit policy

Use line search merit only to decide whether a trial step is better than the current iterate. Use original convergence tolerances to decide actual convergence.

```cpp
double meritResidual() { return resErr / maxResErr; }
double meritEnergy()   { return eneErr / maxEneErr; }
double meritMixed()    { return max(meritResidual(), meritEnergy()); }
```

If one criterion is disabled, exclude it.

## Do not let line search accept loose steps

Line search acceptance means:

```text
this alpha is an improvement and can become the next nonlinear iterate
```

It does not mean:

```text
the load step is converged
```

The load step still must pass:

```text
disErr <= maxDisErr
resErr <= maxResErr
eneErr <= maxEneErr
```

## Recommended fallback policy

For stabilization development:

```text
limit_tolerance 0
```

If fallback acceptance must exist for legacy workflows:

- make it opt-in,
- log it loudly,
- never count it as full convergence,
- do not use it in benchmark comparisons.

Possible new keyword:

```text
allow_unconverged_limit_acceptance 0/1
```

Recommended default:

```text
allow_unconverged_limit_acceptance 0
```

---

# TS-N65 experimental matrix

## Baseline group

```text
G0-A: current settings except limit_tolerance 0 and max_iterations 40
G0-B: adaptive dt enabled, stiffness update 10
G0-C: adaptive dt enabled, stiffness update 3
G0-D: adaptive dt enabled, stiffness update 1
```

## Damping group

```text
G1-A: fixed damping alpha=0.75
G1-B: fixed damping alpha=0.5
G1-C: fixed damping alpha=0.25
G1-D: adaptive damping
```

## Line-search group

```text
G2-A: backtracking mixed merit, max_trials=6, min_alpha=0.03125
G2-B: backtracking residual merit
G2-C: backtracking energy merit
G2-D: backtracking mixed merit + stiffness update 1
G2-E: backtracking mixed merit + stiffness update 3
```

## Stagnation group

```text
G3-A: line search + stagnation cutback N=8 ratio=0.95
G3-B: line search + stagnation cutback N=5 ratio=0.95
G3-C: line search + stagnation cutback N=8 ratio=0.90
```

## Indirect/path-following group

```text
G4-A: indirect control only
G4-B: indirect control + line search
G4-C: indirect control + line search + stagnation cutback
```

## Metrics to report

```text
total accepted steps
total rejected/cutback attempts
mean nonlinear iterations per accepted step
median nonlinear iterations
95th percentile nonlinear iterations
max nonlinear iterations
number of steps hitting max_iterations
number of fallback accepted steps
alpha histogram
line-search failure count
stagnation cutback count
wall-clock time
load-displacement curve
final damage/plasticity summaries
```

## Success threshold

A successful first patch should achieve most of these:

- No fallback accepted steps.
- No 1000-iteration steps.
- Average nonlinear iterations substantially below the observed 353.9.
- Difficult steps cut back or converge in tens of iterations.
- Load-displacement curve remains consistent with baseline where baseline is trustworthy.

---

# Implementation pitfalls and how to avoid them

## Pitfall 1 — Rejected alpha accumulates displacement

Bad pattern:

```cpp
for (alpha : alphas) {
    ddr = alpha * du;
    updateFieldVariables(); // accumulates into trial_r repeatedly
    ...
}
```

Correct pattern:

```cpp
baseTrialR = trial_r_before_alpha_loop;
for (alpha : alphas) {
    restoreState();
    trial_r = baseTrialR;
    ddr = alpha * du;
    buildFullIncrement();
    trial_r += full_ddr;
    evaluateResidual();
}
```

## Pitfall 2 — Using displacement increment error as merit

Small alpha automatically reduces displacement increment norm. A line search based on `disErr` alone will always prefer tiny steps.

Use residual/energy merit.

## Pitfall 3 — Assuming `frozen=true` means side-effect free

Audit material code. `frozen=true` may prevent committed history update but still alter trial caches.

Add rollback tests.

## Pitfall 4 — Accepting `limit_tolerance` results as success

Fallback acceptance is not convergence to the requested tolerance. During development, disable it.

## Pitfall 5 — Expecting line search to solve snapback

Line search improves local convergence. It does not change the fact that load-controlled equilibrium may be unstable after peak. Use displacement/arc-length control for that.

## Pitfall 6 — Rebuilding stiffness too rarely during damage

Softening material response can change stiffness rapidly. `stiffness_matrix_iter_update 10` may be too stale.

Use adaptive rebuild on stagnation or small alpha.

---

# Minimal review checklist for a merge request

## Algorithm

- [ ] The solver can run with all new features disabled and reproduce old behavior.
- [ ] Line search uses a clear merit function.
- [ ] Line search acceptance does not replace convergence criteria.
- [ ] Stagnation cutback uses merit history and logs reason.
- [ ] Cutback interacts correctly with existing `dtmin/dtmax` logic.

## State safety

- [ ] Rejected alpha restores solver vectors.
- [ ] Rejected alpha restores material status or resets to a documented safe boundary.
- [ ] Accepted alpha leaves state equivalent to one accepted update.
- [ ] Step failure restores previous converged step state.
- [ ] No hidden mutation of indirect-control state on rejected trials.

## Inputs and defaults

- [ ] Defaults preserve existing behavior.
- [ ] Invalid values are rejected or clamped with warning.
- [ ] New keywords are documented.
- [ ] Example TS-N65 input block is provided.

## Logging

- [ ] Alpha and trial count are logged.
- [ ] Cutback reason is logged.
- [ ] Fallback acceptance is clearly labeled.
- [ ] Linear solver status is available for diagnosis.

## Tests

- [ ] Unit test: line search chooses alpha < 1.
- [ ] Unit test: rejected trials do not commit material state.
- [ ] Regression: features disabled matches old solver.
- [ ] Regression: adaptive cutback triggers before max iterations.
- [ ] TS-N65 comparison table included.

---

# Suggested documentation block for users

```text
# Nonlinear globalization controls
#
# nonlinear_line_search:
#   off            Apply full Newton increment, legacy behavior.
#   backtracking   Try alpha=1, reduction*alpha, ... and accept first residual/energy merit decrease.
#
# nonlinear_line_search_merit:
#   residual       Use normalized residual error.
#   energy         Use normalized energy error.
#   mixed          Use max(residual/tol, energy/tol). Recommended.
#
# nonlinear_line_search_freeze_material:
#   1              Evaluate trial residuals with frozen internal variables where supported.
#                  Recommended for path-dependent material line-search trials.
#
# nonlinear_line_search_cutback_on_fail:
#   1              If no alpha gives merit decrease, reject current load step and reduce dt.
#
# nonlinear_stagnation_cutback:
#   1              Reject and cut back a load step if nonlinear merit stagnates for N iterations.
```

Example:

```text
nonlinear_line_search backtracking
nonlinear_line_search_merit mixed
nonlinear_line_search_reduction 0.5
nonlinear_line_search_min_alpha 0.03125
nonlinear_line_search_max_trials 6
nonlinear_line_search_accept_any_decrease 1
nonlinear_line_search_freeze_material 1
nonlinear_line_search_cutback_on_fail 1

nonlinear_stagnation_cutback 1
nonlinear_stagnation_iterations 8
nonlinear_stagnation_ratio 0.95
nonlinear_growth_cutback 1.25

max_iterations 40
min_time_step 1.25e-5
max_time_step 1.25e-3
limit_tolerance 0
```

---

# Final recommended roadmap

## Week 1-style patch

Implement:

```text
instrumentation
limit_tolerance diagnostic logging
stagnation cutback
fixed damping
backtracking line search with frozen residuals
```

Do not implement:

```text
arc-length
material-law changes
PETSc SNES integration
new consistent tangents
```

## Week 2-style patch

Implement:

```text
material-status snapshot/restore
automatic stiffness rebuild on small alpha/stagnation
finite-difference tangent tests
line-search test suite
```

## Later patches

Implement based on evidence:

```text
trust-region/LM if line search still fails locally
indirect/displacement control if load control is unstable
arc-length if post-peak snapback must be traced
material substepping/viscosity/localization regularization if constitutive/localization pathology dominates
```

---

# Bottom line

The first serious nonlinear solver enhancement should be a **state-safe backtracking line search plus early stagnation cutback**. It fits the current architecture, directly addresses the observed failure mode, and creates the trial-evaluation/rollback infrastructure needed for all later methods. Arc-length, viscosity, material substepping, and regularization remain important, but they should come after the local globalization layer is correct and instrumented.
