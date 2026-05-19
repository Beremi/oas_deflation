# Nonlinear solver tangent attribution packet 1 report

Generated: 2026-05-19T12:47:40

## Scope

Implemented the first checkpoint diagnostic: element/material tangent attribution on top of the existing global nonlinear tangent finite-difference check. No production nonlinear convergence behavior is changed when `nonlinear_tangent_check 0`.

## Code changes

- Added `nonlinear_tangent_check_scope global|element_top`.
- Added `nonlinear_tangent_check_top_elements`.
- Added `nonlinear_tangent_check_element_output`.
- Added `nonlinear_tangent_check_matrix_type current|<matrix type>`.
- Added element-level local `Ke * p_e` versus local finite-difference internal-force contribution ranking.
- Optimized element attribution to restore only element-local material status inside each local FD test instead of restoring the full TS-N65 material snapshot for every element.
- Updated `docs/nonlinear_solver_next_steps_workplan.md` with the mandatory strict 8-quarter-step phase-closure rule and diagnostic stop-after note.

## Build and smoke

- Build used `/tmp/oas_deflation_build_hypre_novtk` with `USE_HYPRE=ON`, `USE_VTK=OFF`: passed.
- Default small smoke with diagnostics disabled: passed earlier in this packet.
- Small `element_top` smoke on `CantileverMechanicsNonlinear`: wrote both `tangent_check.tsv` and `tangent_check_elements.tsv`; local element errors were near machine precision.

## Global tangent check at TS-N65 hard state

Case: `diag-step6-it10-localrestore`, step 6 iteration 10, eps `1e-6`, matrix type `tangent`.

|direction|relerr|cosine|kp|fd|residual|
|---|---|---|---|---|---|
|newton|5.32096|0.707442|1.18658e+07|1.98382e+06|618.181|
|random_1|0.0329979|0.999465|2.92288e+09|2.91168e+09|618.181|

Interpretation: the Newton direction is the problematic one: relative error about `5.32` and cosine about `0.707`. The deterministic random direction is much better: relative error about `0.033` and cosine about `0.9995`.

## Top element/material contributors

Top Newton-direction mismatch rows:

|rank|elem|name|relerr|cosine|kp|fd|mismatch|
|---|---|---|---|---|---|---|---|
|1|300453|LDPMTetra|3.70595|-0.160985|330284|96825.5|358831|
|2|297651|LDPMTetra|1.71968|0.489191|327755|166265|285922|
|3|62838|LDPMTetra|2.3145|-0.0495626|244293|119849|277390|
|4|300222|LDPMTetra|4.74441|-0.228832|251143|56888.8|269904|
|5|300651|LDPMTetra|2.59754|-0.30437|218227|103317|268370|
|6|300102|LDPMTetra|4.15667|-0.836564|208600|63523.2|264045|
|7|510156|LDPMTetra|3.62304|-0.410115|218189|70468.6|255311|
|8|509344|LDPMTetra|5.93169|-0.67327|213638|40988.4|243130|
|9|69628|LDPMTetra|3.0373|-0.20246|209387|78344.4|237956|
|10|300092|LDPMTetra|3.99671|-0.889535|182077|59097.2|236194|

Top random-direction mismatch rows:

|rank|elem|name|relerr|cosine|kp|fd|mismatch|
|---|---|---|---|---|---|---|---|
|1|167177|LDPMTetra|2.37514|0.161229|1.24654e+06|536927|1.27527e+06|
|2|354562|LDPMTetra|2.34426|0.25022|1.23514e+06|517836|1.21394e+06|
|3|164317|LDPMTetra|1.6522|0.532929|1.38906e+06|711608|1.17572e+06|
|4|351733|LDPMTetra|1.38239|0.550144|1.39492e+06|844481|1.1674e+06|
|5|356636|LDPMTetra|1.93657|0.279638|1.1798e+06|601494|1.16484e+06|
|6|508678|LDPMTetra|1.21095|0.641937|1.43146e+06|906441|1.09766e+06|
|7|346009|LDPMTetra|1.9822|0.464403|1.23415e+06|551510|1.0932e+06|
|8|137873|LDPMTetra|2.00285|0.320231|1.12224e+06|538278|1.07809e+06|
|9|75029|LDPMTetra|3.6801|0.0043141|1.02845e+06|290035|1.06736e+06|
|10|343798|LDPMTetra|1.71868|0.399487|1.14157e+06|615981|1.05867e+06|

Interpretation: the top contributors are `LDPMTetra` elements. The mismatch is not concentrated in `NormalPlasticBeam` or boundary/control elements in this check. That points the next checkpoint toward the LDPM/CSL constitutive tangent branches.

## Diagnostic-enabled run status

The diagnostic-enabled run was stopped after reproducing the step-6 blow-up/stall. It generated the tangent attribution files successfully.

|step|dt|iters|duration|conv|res|disp|energy|
|---|---|---|---|---|---|---|---|
|1|0.00125|8|00:01:03.500|yes|2.35302e-07|7.85526e-05|2.44221e-05|
|2|1.250000e-03|9|00:00:41.448|yes|9.89544e-05|0.000260349|0.00053777|
|3|1.250000e-03|12|00:01:02.163|yes|0.000379676|0.000100366|0.000740182|
|4|1.250000e-03|14|00:01:05.889|yes|0.00048679|9.8543e-05|0.00080374|
|5|1.250000e-03|18|00:01:13.659|yes|0.000553742|8.63182e-05|0.000894886|
|6|1.250000e-03|132|stalled/killed|no|0.807411|1|1.43619|

Diagnostic run converged steps before stop: `5/6`; parsed nonlinear iterations: `193`.

## Diagnostic-disabled strict 8 quarter-step phase closure

Settings: `time_step = min_time_step = max_time_step = 1.25e-3`, `total_time = 1e-2`, `max_iterations = 300`, `stiff_matrix_type tangent`, rebuild every 10 nonlinear iterations, DFGMRES/HYPRE with tolerance `1e-1`, `N=20`, `j=16`, `nonlinear_tangent_check 0`.

|step|dt|iters|duration|conv|res|disp|energy|
|---|---|---|---|---|---|---|---|
|1|0.00125|8|00:00:54.810|yes|2.35302e-07|7.85526e-05|2.44221e-05|
|2|1.250000e-03|9|00:00:37.177|yes|9.89544e-05|0.000260349|0.00053777|
|3|1.250000e-03|12|00:00:53.544|yes|0.000379676|0.000100366|0.000740182|
|4|1.250000e-03|14|00:00:57.027|yes|0.00048679|9.8543e-05|0.00080374|
|5|1.250000e-03|18|00:01:06.505|yes|0.000553742|8.63182e-05|0.000894886|
|6|1.250000e-03|133|stalled/killed|no|1|0.00977435|1.41421|

Phase closure result: `5/8` strict steps completed before stop; parsed nonlinear iterations before stop: `194`.

Baseline target from saved notes: `8/8`, `585` nonlinear iterations, `26:53` total duration. This local run did **not** reproduce that target: it stalled in step 6 after iteration 132 with residual near `1`, energy near `1.414`, repeated DFGMRES warnings, and no further nonlinear output for more than 10 minutes. I killed it and marked the phase as non-converged/stalled.

## Conclusions

- The new element attribution diagnostic works and produces ranked local mismatch rows.
- The Newton-direction tangent mismatch at step 6 iteration 10 is dominated by `LDPMTetra` rows, so checkpoint 2 should audit LDPM/CSL material tangent branches first.
- The random direction is much closer to the assembled tangent than the Newton direction, which suggests the bad behavior is path/state/localization-direction dependent rather than a uniform matrix assembly sign error.
- The strict 8-quarter-step phase closure failed locally with diagnostics disabled, so the current local executable/settings do not reproduce the saved `585 iteration / 26:53` baseline. This must be resolved before claiming any solver improvement.
- Until continuation-after-diagnostic is proven non-intrusive, use `nonlinear_tangent_check_stop_after 1` for `element_top` hard-state diagnostics and run phase closure separately with diagnostics disabled.

## Recommended next checkpoint

Checkpoint 2: targeted LDPM/CSL tangent audit tests. Start with the branches used by the top `LDPMTetra` elements from the Newton-direction table. Compare material stress update derivatives against finite differences for loading/unloading/damage-growth states, with frozen-state and actual-state variants.
