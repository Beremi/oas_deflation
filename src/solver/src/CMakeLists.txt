set (solver
    boundary_condition.cpp
    data_exporter.cpp
    element_container.cpp
    element.cpp
    linear_algebra.cpp
    #main.cpp
    material_container.cpp
    material.cpp
    mars_material.cpp
    fatigue_material.cpp
    node_container.cpp
    node.cpp
    solver.cpp
    )

    add_library(solver
        ${LIB_TYPE}
        ${solver}
        )
