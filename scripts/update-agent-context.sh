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
```

Default build directory:

```text
../oas_deflation-build/release
```

Dogbone input:

```text
data/cases/Dogbone/master.inp
```
EOF

cat > "$ctx/repo_map.md" <<'EOF'
# Repo Map

- `src/solver/src/linalg.h`, `src/solver/src/linalg.cpp`: linear solver interface and solver wrappers.
- `src/solver/src/solver_implicit.cpp`: implicit solve loops, matrix update policy, factorization calls, and `linalgsolver->solve(...)` call sites.
- `docs/deflation/`: human-facing deflation project documentation.
- `scripts/update-agent-context.sh`: regenerates ignored agent context.
- `data/`: ignored benchmark archives, extracted cases, and run outputs.
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
} > "$ctx/benchmark_state.md"

cat > "$ctx/decision_log.md" <<'EOF'
# Decision Log

- Repository shape: OAS source at repository root, with `upstream` pointing to `https://gitlab.com/kelidas/OAS.git`.
- Data policy: Dogbone and TS-N_65 archives plus extracted cases are local-only under ignored `data/`.
- Build policy: use out-of-source CMake through the root Makefile; default generator is Unix Makefiles.
- Solver policy: bootstrap makes no OAS C++ behavior changes. First implementation work should add instrumentation, then PCG/IC, then AMG, then deflation/recycling.
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

After meaningful changes to source, build configuration, benchmark state, or test results, run:

```sh
make agent-context
```

Then update `.agent_context/current_state.md` manually if the generated summary misses important context.
EOF

printf 'Updated agent context in %s\n' "$ctx"
