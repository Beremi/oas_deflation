#!/usr/bin/env bash
set -euo pipefail

case_dir="${1:?case directory is required}"
master_file="${2:?master input file is required}"
oas_bin="${3:?OAS executable is required}"
mklroot="${4:?MKLROOT is required}"
threads="${5:?thread count is required}"
solver="${6:?solver type is required}"

solver_inp="$case_dir/solver.inp"
profile_dir="$case_dir/results"
backup="$(mktemp)"

restore_solver_input() {
  if [[ -f "$backup" ]]; then
    cp "$backup" "$solver_inp"
    rm -f "$backup"
  fi
}
trap restore_solver_input EXIT

if [[ ! -f "$solver_inp" ]]; then
  echo "Missing Dogbone solver input: $solver_inp" >&2
  exit 1
fi
if [[ ! -x "$oas_bin" ]]; then
  echo "Missing or non-executable OAS binary: $oas_bin" >&2
  exit 1
fi

cp "$solver_inp" "$backup"
mkdir -p "$profile_dir"
rm -f "$profile_dir/linear_profile_events.tsv" "$profile_dir/linear_profile_iterations.tsv"

python3 - "$solver_inp" "$solver" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
solver = sys.argv[2]
updates = {
    "solver_type": solver,
    "linear_solver_profile": "1",
    "linear_solver_profile_matrix_delta": "1",
    "linear_solver_profile_file": "linear_profile",
}

lines = path.read_text(encoding="utf-8").splitlines()
seen = set()
out = []

for line in lines:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        out.append(line)
        continue
    key = stripped.split()[0]
    if key in updates:
        out.append(f"{key}\t{updates[key]}")
        seen.add(key)
    else:
        out.append(line)

missing = [key for key in updates if key not in seen]
if missing:
    if out and out[-1].strip():
        out.append("")
    out.append("# Linear-solve profiling; added by scripts/run-dogbone-profile.sh for this local run.")
    for key in missing:
        out.append(f"{key}\t{updates[key]}")

path.write_text("\n".join(out) + "\n", encoding="utf-8")
PY

MKLROOT="$mklroot" MKL_NUM_THREADS="$threads" OMP_NUM_THREADS="$threads" \
  "$oas_bin" -j "$threads" "$master_file"
