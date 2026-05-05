#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ctx="$root/.agent_context"
mkdir -p "$ctx"

timestamp="$(date -Is)"
head_sha="$(git -C "$root" rev-parse --short=12 HEAD 2>/dev/null || true)"
upstream_sha="$(git -C "$root" rev-parse --short=12 upstream/master 2>/dev/null || true)"
branch="$(git -C "$root" branch --show-current 2>/dev/null || true)"

{
  printf '# Current State\n\n'
  printf -- '- Updated: `%s`\n' "$timestamp"
  printf -- '- Branch: `%s`\n' "${branch:-unknown}"
  printf -- '- HEAD: `%s`\n' "${head_sha:-unknown}"
  printf -- '- Upstream master: `%s`\n' "${upstream_sha:-unknown}"
  printf -- '- Build dir: `%s`\n' "../oas_deflation-build/release"
  printf -- '- OAS executable: `%s`\n\n' "../oas_deflation-build/release/bin/OAS"
  printf '## Git Status\n\n```text\n'
  git -C "$root" status --short --branch || true
  printf '```\n'
} > "$ctx/current_state.md"

cat > "$ctx/build_run.md" <<'EOF'
# Build and Run

Use the root Makefile.

```sh
make configure
make build
make dogbone
make dogbone-solver SOLVER=EigenLDLT
make dogbone-solver SOLVER=PardisoLDLT
make dogbone-profile USE_VTK=OFF THREADS=4 SOLVER=EigenLDLT
make build USE_VTK=OFF USE_HYPRE=ON THREADS=8
make dogbone-amgcl-prefix USE_VTK=OFF THREADS=8 AMGCL_TOL=1e-6 AMGCL_MAXIT=500
make dogbone-hypre-prefix USE_VTK=OFF USE_HYPRE=ON THREADS=8 HYPRE_TOL=1e-6 HYPRE_MAXIT=500
make dogbone-amgcl-profile USE_VTK=OFF THREADS=8 AMGCL_TOL=1e-6 AMGCL_MAXIT=500
make tsn65-amgcl-prefix USE_VTK=OFF THREADS=8 AMGCL_TOL=1e-6 AMGCL_MAXIT=500
make tsn65-hypre-prefix USE_VTK=OFF USE_HYPRE=ON THREADS=8 HYPRE_TOL=1e-6 HYPRE_MAXIT=500
```

Default build directory:

```text
../oas_deflation-build/release
```

Dogbone input:

```text
data/cases/Dogbone/master.inp
```

Linear profiling:

```sh
make linear-profile-report PROFILE_DIR=data/cases/Dogbone/results OUT_DIR=results/dogbone-eigenldlt-manual SOLVER=EigenLDLT
scripts/sweep-amgcl.py --mode dogbone-prefix --threads 8 --tolerances 1e-6 --max-iterations 500 --eps-strong 0 --npre 1 --npost 1
```

Raw profile TSVs live in `data/cases/Dogbone/results/`. Markdown/PNG exchange reports live in ignored `results/`.
EOF

cat > "$ctx/repo_map.md" <<'EOF'
# Repo Map

- `src/solver/src/linalg.h`, `src/solver/src/linalg.cpp`: linear solver interface and solver wrappers.
- `src/solver/src/linalg_profile.h`, `src/solver/src/linalg_profile.cpp`: optional linear-solve profiling data collector.
- `src/solver/src/solver_implicit.cpp`: implicit solve loops, matrix update policy, factorization calls, and `linalgsolver->solve(...)` call sites.
- `docs/deflation/`: human-facing deflation project documentation.
- `scripts/analyze-linear-profile.py`: converts profiler TSV files and `solver.out` into Markdown/PNG reports.
- `scripts/run-dogbone-profile.sh`: local Dogbone profiling run wrapper that temporarily edits ignored `solver.inp`.
- `scripts/run-oas-profile.sh`: generic local profiling wrapper for Dogbone and TS-N_65; temporarily edits ignored `solver.inp`, restores it, and copies raw/profile artifacts to ignored `results/`.
- `scripts/sweep-amgcl.py`: local AMGCL parameter/tolerance sweep driver with Markdown and TSV summaries.
- `scripts/update-agent-context.sh`: regenerates ignored agent context.
- `data/`: ignored benchmark archives, extracted cases, and run outputs.
- `results/`: ignored local experiment exchange reports; check the newest `results/*/linear-profile.md`.
EOF

{
  printf '# Benchmark State\n\n'
  printf '## Local Archives\n\n'
  if [[ -d "$root/data/archives" ]]; then
    find "$root/data/archives" -maxdepth 1 -type f -printf '- `%f` (%s bytes)\n' | sort
  else
    printf -- '- No `data/archives` directory found.\n'
  fi
  printf '\n## Extracted Cases\n\n'
  if [[ -d "$root/data/cases" ]]; then
    find "$root/data/cases" -mindepth 1 -maxdepth 1 -type d -printf '- `%f`\n' | sort
  else
    printf -- '- No `data/cases` directory found.\n'
  fi
  printf '\n## Latest Linear Profile Reports\n\n'
  if [[ -d "$root/results" ]]; then
    mapfile -t reports < <(find "$root/results" -mindepth 2 -maxdepth 2 -name 'linear-profile.md' -printf '%T@ %p\n' | sort -nr | head -5 | cut -d' ' -f2-)
    if (( ${#reports[@]} )); then
      for report in "${reports[@]}"; do
        rel="${report#"$root"/}"
        title="$(sed -n '1s/^# //p' "$report" 2>/dev/null || true)"
        printf -- '- `%s`' "$rel"
        if [[ -n "$title" ]]; then
          printf ' - %s' "$title"
        fi
        printf '\n'
      done
    else
      printf -- '- No `results/*/linear-profile.md` reports found.\n'
    fi
  else
    printf -- '- No local `results/` directory found.\n'
  fi
  printf '\n## Latest Summary Reports\n\n'
  if [[ -d "$root/results" ]]; then
    mapfile -t summaries < <(find "$root/results" -mindepth 2 -maxdepth 2 -name 'summary.md' -printf '%T@ %p\n' | sort -nr | head -5 | cut -d' ' -f2-)
    if (( ${#summaries[@]} )); then
      for summary in "${summaries[@]}"; do
        rel="${summary#"$root"/}"
        title="$(sed -n '1s/^# //p' "$summary" 2>/dev/null || true)"
        printf -- '- `%s`' "$rel"
        if [[ -n "$title" ]]; then
          printf ' - %s' "$title"
        fi
        printf '\n'
      done
    else
      printf -- '- No `results/*/summary.md` reports found.\n'
    fi
  else
    printf -- '- No local `results/` directory found.\n'
  fi
} > "$ctx/benchmark_state.md"

cat > "$ctx/decision_log.md" <<'EOF'
# Decision Log

- Repository shape: OAS source at repository root, with `upstream` pointing to `https://gitlab.com/kelidas/OAS.git`.
- Data policy: Dogbone and TS-N_65 archives plus extracted cases are local-only under ignored `data/`.
- Build policy: use out-of-source CMake through the root Makefile; default generator is Unix Makefiles.
- Solver policy: first implementation work adds optional instrumentation, then PCG/IC, then AMG, then deflation/recycling.
- Profiling policy: `linear_solver_profile 1` writes raw TSVs into the OAS case `results/` directory; generated Markdown/PNG reports stay local-only under ignored `results/`.
- AMG preconditioner policy: Dogbone first-5 is the gate before full Dogbone/TS-N_65. Current passing gate report: `results/amgcl-hypre-dogbone-gate-20260503-054622/summary.md` with AMGCL elastic block CG at 6-7 iterations and hypre BoomerAMG-CG at 12-15 iterations for `1e-6` true relative residual.
- Dogbone full policy: current full Dogbone comparison is `results/amgcl-hypre-dogbone-full-20260503-055920/summary.md`; both AMGCL elastic block CG and hypre BoomerAMG-CG preserve the Pardiso accepted-step shape at `1e-6`.
- TS-N_65 policy: current first-step gate report is `results/tsn65-amg-firststep-gate-20260503-062440/summary.md`. hypre passes with iterations 56,117,92,68,68; AMGCL remains blocked because block ILU0 gives NaNs and scalar/block-SPAI variants hit 500 iterations.
- hypre policy: use nodal BoomerAMG with the lifted node-major elastic space. Keep explicit `dof_func` and interpolation-vector calls disabled until an ownership-safe wrapper is implemented.
- Agent policy: keep `agents.md` and `.agent_context/` ignored; refresh `.agent_context/` after meaningful build, run, or code changes.
EOF

cat > "$root/agents.md" <<'EOF'
# Agent Instructions

This repository is an at-root OAS fork for linear solver deflation work.

Read these files first:

1. `README.md`
2. `docs/deflation/implementation-report.md`
3. `.agent_context/current_state.md`
4. `.agent_context/build_run.md`
5. `.agent_context/repo_map.md`
6. `.agent_context/benchmark_state.md`
7. `.agent_context/decision_log.md`

Keep OAS C++ changes focused. The first solver implementation phase should instrument current solver behavior before adding new algorithms.

Before answering profiling questions, check the newest local report matching:

```sh
ls -td results/* 2>/dev/null | head
```

Then read `results/<latest>/linear-profile.md` if it exists.

For AMGCL/hypre work, also check the newest local summary matching:

```sh
ls -td results/*amg* results/*hypre* 2>/dev/null | head
```

Read `summary.md` first when present, then drill into the linked `linear-profile.md` reports and raw TSVs.

After meaningful changes to source, build configuration, benchmark state, or test results, run:

```sh
make agent-context
```

Then update `.agent_context/current_state.md` manually if the generated summary misses important context.
EOF

printf 'Updated agent context in %s\n' "$ctx"
