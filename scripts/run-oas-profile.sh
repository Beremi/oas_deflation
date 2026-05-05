#!/usr/bin/env bash
set -euo pipefail

case_dir="${1:?case directory is required}"
master_file="${2:?master input file is required}"
oas_bin="${3:?OAS executable is required}"
threads="${4:?thread count is required}"
solver="${5:?solver type is required}"
title="${6:?report title is required}"
case_name="${7:?case name is required}"
out_dir="${8:?output directory is required}"
shift 8

solver_inp="$case_dir/solver.inp"
profile_dir="$case_dir/results"
backup="$(mktemp)"
timeout_duration="${OAS_TIMEOUT:-12h}"
replay_dir="linear_replay"
for item in "$@"; do
  if [[ "$item" == linear_solver_replay_dir=* ]]; then
    replay_dir="${item#linear_solver_replay_dir=}"
  fi
done

restore_solver_input() {
  if [[ -f "$backup" ]]; then
    cp "$backup" "$solver_inp"
    rm -f "$backup"
  fi
}
trap restore_solver_input EXIT

if [[ ! -f "$solver_inp" ]]; then
  echo "Missing solver input: $solver_inp" >&2
  exit 1
fi
if [[ ! -x "$oas_bin" ]]; then
  echo "Missing or non-executable OAS binary: $oas_bin" >&2
  exit 1
fi

mkdir -p "$profile_dir" "$out_dir"
cp "$solver_inp" "$backup"
cp "$backup" "$out_dir/solver.inp.original"

python3 - "$solver_inp" "$solver" "$@" <<'PY'
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

for item in sys.argv[3:]:
    if "=" not in item:
        raise SystemExit(f"Expected key=value override, got {item!r}")
    key, value = item.split("=", 1)
    updates[key] = value

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
    out.append("# Local profiling settings; added temporarily by scripts/run-oas-profile.sh.")
    for key in missing:
        out.append(f"{key}\t{updates[key]}")

path.write_text("\n".join(out) + "\n", encoding="utf-8")
PY

cp "$solver_inp" "$out_dir/solver.inp.effective"
rm -f \
  "$profile_dir/linear_profile_events.tsv" \
  "$profile_dir/linear_profile_iterations.tsv" \
  "$profile_dir/runtime_profile_events.tsv" \
  "$profile_dir/runtime_profile_summary.tsv" \
  "$profile_dir/solver.out"
rm -rf "$profile_dir/$replay_dir"

{
  printf 'Starting OAS profile run at %s\n' "$(date -Is)"
  printf 'Case: %s\nSolver: %s\nThreads: %s\nTimeout: %s\n' "$case_name" "$solver" "$threads" "$timeout_duration"
} | tee "$out_dir/oas.log"

set +e
MKLROOT="${MKLROOT:-/opt/intel/oneapi/mkl/latest}" MKL_NUM_THREADS="$threads" OMP_NUM_THREADS="$threads" \
  timeout --signal=INT --kill-after=10m "$timeout_duration" "$oas_bin" -j "$threads" "$master_file" >> "$out_dir/oas.log" 2>&1
run_status=$?
set -e

printf 'OAS exit status: %s\nFinished at %s\n' "$run_status" "$(date -Is)" | tee -a "$out_dir/oas.log"

cp -f "$profile_dir/linear_profile_events.tsv" "$out_dir/" 2>/dev/null || true
cp -f "$profile_dir/linear_profile_iterations.tsv" "$out_dir/" 2>/dev/null || true
cp -f "$profile_dir/runtime_profile_events.tsv" "$out_dir/" 2>/dev/null || true
cp -f "$profile_dir/runtime_profile_summary.tsv" "$out_dir/" 2>/dev/null || true
cp -f "$profile_dir/solver.out" "$out_dir/" 2>/dev/null || true
if [[ -d "$profile_dir/$replay_dir" ]]; then
  rm -rf "$out_dir/$replay_dir"
  cp -a "$profile_dir/$replay_dir" "$out_dir/"
fi

if [[ -f "$profile_dir/linear_profile_events.tsv" ]]; then
  python3 scripts/analyze-linear-profile.py \
    --profile-dir "$profile_dir" \
    --out-dir "$out_dir" \
    --title "$title" \
    --case "$case_name" \
    --threads "$threads"
else
  printf 'No linear_profile_events.tsv produced.\n' | tee -a "$out_dir/oas.log"
fi

exit "$run_status"
