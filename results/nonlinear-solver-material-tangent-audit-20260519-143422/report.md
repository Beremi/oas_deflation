# Nonlinear Solver Material Tangent Audit

Date: 2026-05-19
Branch: `nonlinear_solver_testing`
Build: `/tmp/oas_deflation_build_material_audit`

## Scope

This checkpoint continues the LDPM/CSL tangent investigation from `intro.md`.

Implemented:

- typed material-state clone/restore/hash support for LDPM-family statuses and coupled CSL statuses,
- a local `OAS_material_tangent_audit` executable,
- a `MaterialTangentAudit` CTest entry,
- element attribution output metadata for `material_id`, `material_name`, and compact material-status names,
- explicit CSL and LDPM `consistent` tangent branches using state-restored numerical material derivatives.
- an archived CSL analytical active-damage tangent branch, `archived_csl_damage_tangent`, kept out of the default `tangent` path because it did not close the TS-N65 phase benchmark in later testing.

No production nonlinear control is enabled by default by this checkpoint.

## Verification Commands

```bash
cmake -S . -B /tmp/oas_deflation_build_material_audit -DCMAKE_BUILD_TYPE=Release -DUSE_VTK=OFF -DUSE_HYPRE=ON
cmake --build /tmp/oas_deflation_build_material_audit --target OAS OAS_material_tangent_audit -j 16
/tmp/oas_deflation_build_material_audit/bin/OAS_material_tangent_audit --output results/nonlinear-solver-material-tangent-audit-20260519-143422/material_tangent_audit.tsv
ctest --test-dir /tmp/oas_deflation_build_material_audit -R MaterialTangentAudit --output-on-failure
```

Result:

```text
material tangent audit passed expected linear/frozen checks
MaterialTangentAudit: Passed, 0.04 sec
```

The only compile warning seen during the final rebuild is the pre-existing `r` loop-variable shadowing warning in `maybeRunNonlinearTangentCheck`.

## Audit Harness

The executable builds one-point material-level probes using `MaterialTestElement`.

For each material/case/branch it compares:

```text
D_material * perturbation
```

against:

```text
[stress(strain + eps * perturbation) - stress(strain)] / eps
```

The harness snapshots/restores material state around each probe. CTest fails only if expected linear or frozen paths drift. Active damage-growth mismatches are recorded as diagnostic findings because they are the current suspected defect.

Cases:

- `elastic_loading`
- `damage_growth`
- `damage_growth_frozen`
- `damaged_unloading`
- `damaged_unloading_frozen`

Branches:

- CSL: `elastic`, `secant`, `tangent`, `archived_csl_damage_tangent`, `consistent`; `consistent` is an actual-state numerical material tangent and is marked not applicable for frozen audit rows.
- LDPM: `elastic`, `secant`, `consistent`; `consistent` is an actual-state numerical material tangent and is marked not applicable for frozen audit rows; `tangent` is unsupported.

## Key Findings

Machine-readable output:

```text
results/nonlinear-solver-material-tangent-audit-20260519-143422/material_tangent_audit.tsv
```

Selected rows:

| material | case | evaluation | branch | relerr | cosine | note |
|---|---|---|---|---:|---:|---|
| CSLMaterial | elastic_loading | actual | tangent | 2.08e-13 | 1.0000 | closes |
| CSLMaterial | damage_growth | actual | secant | 3.4089 | -0.9963 | active damage mismatch |
| CSLMaterial | damage_growth | actual | tangent | 3.4089 | -0.9963 | legacy default tangent mismatch |
| CSLMaterial | damage_growth | actual | archived_csl_damage_tangent | 2.38e-5 | 1.0000 | archived analytical active-damage tangent closes |
| CSLMaterial | damage_growth | actual | consistent | 3.01e-5 | 1.0000 | numerical consistent branch closes |
| CSLMaterial | damage_growth_frozen | frozen | tangent | 6.12e-13 | 1.0000 | closes when frozen |
| CSLMaterial | damaged_unloading | actual | tangent | 7.77e-13 | 1.0000 | closes |
| LDPMMaterial | elastic_loading | actual | secant | 186.796 | 0.99996 | diagnostic mismatch |
| LDPMMaterial | elastic_loading | actual | consistent | 2.32e-6 | 1.0000 | numerical consistent branch closes |
| LDPMMaterial | damage_growth | actual | secant | 8.76e-12 | 1.0000 | closes in this local path |
| LDPMMaterial | damaged_unloading | actual | secant | 9.66e-12 | 1.0000 | closes |

The strongest local CSL finding is specific: the CSL secant branch and the default `tangent` branch agree with finite differences for elastic, frozen, and damaged-unloading probes, but not for active damage growth. The archived analytical branch adds the missing active damage derivative and closes the active damage-growth probe at the same accuracy level as the numerical `consistent` reference. It remains archived because the later TS-N65 8-quarter phase runs showed diagnostic improvement but no strict phase closure.

The LDPM one-point probe exposes a large mismatch in the selected `elastic_loading` path for the existing `secant` branch. The new LDPM `consistent` branch closes that same path, while the later damage/unloading paths already close with `secant`. The later TS-N65 rerun showed that the hard-state top `LDPMTetra` rows were backed by CSL statuses, so the LDPM local mismatch remains a separate diagnostic item rather than the current TS-N65 culprit.

## Element Attribution Metadata

`nonlinear_tangent_check_scope element_top` now writes:

```text
material_id
material_name
material_status_names
```

in addition to `element_id` and `element_name`.

The status summary de-duplicates repeated integration-point statuses and includes the mechanical substatus when a coupled status wraps one, for example:

```text
CoupledCSLMaterialStatus->CSLMaterialStatusx4
```

This closes the checklist item that the previous top-row packet could only identify `LDPMTetra`, not the exact material/status classes behind those rows. The previous hard-state run did not contain these columns, so it should be rerun before making a final attribution claim.

## Files Changed

Code:

```text
CMakeLists.txt
src/solver/src/material_csl.cpp
src/solver/src/material_csl.h
src/solver/src/material_ldpm.cpp
src/solver/src/material_ldpm.h
src/solver/src/material_tangent_audit_main.cpp
src/solver/src/solver_implicit.cpp
src/tests/CMakeLists.cmake
```

Results:

```text
results/nonlinear-solver-material-tangent-audit-20260519-143422/material_tangent_audit.tsv
results/nonlinear-solver-material-tangent-audit-20260519-143422/report.md
```

## Phase Closure

Strict TS-N65 phase closure was not rerun for this checkpoint.

Reason: this material-audit checkpoint originally preceded the TS-N65 zip rerun. The later benchmark rerun confirmed that the top hard-state `LDPMTetra` rows were backed by `CSL material` / `CSL mat. statusx12`, and the archived analytical CSL branch improved local diagnostics but still failed strict phase closure.

An attempted broader CTest smoke command was not useful: legacy mechanics tests hit their existing 10-second CTest timeout, and a direct `SpringMechElastic` harness invocation finished OAS quickly but failed because expected result/check directories were missing in that build/test invocation.

## Next Steps

1. Keep the archived analytical CSL branch available for reproduction as `archived_csl_damage_tangent`, but keep it out of the normal `tangent` path.
2. Treat numerical `consistent` and `archived_csl_damage_tangent` as diagnostics/references, not proven production matrix choices.
3. Move the next solver experiment toward path-following/indirect control, tangent limiting/blending, or constitutive regularization.
4. Rerun strict diagnostic-disabled TS-N65 8-quarter-step phase closure against the saved `8/8`, `585` iteration, `26:53` baseline before claiming solver improvement.
