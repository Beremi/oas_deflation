set (solver_files
    model.cpp
    boundary_condition.cpp
    constraint.cpp
    cross_section.cpp
    data_exporter.cpp
    element_beam.cpp
    element_container.cpp
    element.cpp
    element_discrete.cpp
    element_fiber.cpp
    element_ldpm.cpp
    element_continuous.cpp
    element_polyhedral.cpp
    geometry.cpp
    globals.cpp
    linalg.cpp
    material_container.cpp
    material_coulomb_friction.cpp
    material.cpp
    material_beam.cpp
    material_vectorial.cpp
    material_fiber.cpp
    material_misc.cpp
    material_plasticity.cpp
    material_slide_3_2.cpp
    material_rve.cpp
    material_csl.cpp
    material_tensorial.cpp
    material_thermomechanical.cpp
    material_ldpm.cpp
    material_fatigue.cpp
    node_container.cpp
    node.cpp
    solver.cpp
    solver_eigenvalue.cpp
    solver_explicit.cpp
    solver_implicit.cpp
    vtk_exporter.cpp
    exporter_model.cpp
    constraint.cpp
    pblock_constraints.cpp
    pblock_periodic_bc.cpp
    pblock_rebar.cpp
    pblock_container.cpp
    function.cpp
    indirect_control.cpp
    shape_functions.cpp
    integration.cpp
    simplex.cpp
    material_htc.cpp
    )

set (superelem_files
    element_superelem.cpp
    )

if (USE_ML_TORCH)
    list(APPEND solver_files material_neuralnetwork.cpp)
    add_library(ml_torch
        ${LIB_TYPE}
        ${superelem_files}
        )
endif ()

add_library(solver
    ${LIB_TYPE}
    ${solver_files}
    )
