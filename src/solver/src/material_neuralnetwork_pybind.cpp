#include "material.h"
#include "element.h"
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <Eigen/Dense>

namespace py = pybind11;

int test(fs :: path python_script_path) {
    py :: scoped_interpreter guard{};  // start Python interpreter

    Eigen :: VectorXd x(13);
    x.setZero();
    Eigen :: VectorXf xf = x.cast< float >();

    try {
        // Add path to where model_def.py is
        py :: module sys = py :: module :: import("sys");
        sys.attr("path").attr("append")( python_script_path.parent_path().string() );

        // Import the Python module
        py :: object scope = py :: module :: import("__main__").attr("__dict__");
        py :: eval_file(python_script_path.string(), scope);

        // Grab model object from the file
        py :: object model = scope [ "model" ];

        // Create input tensor in Python
        py :: module torch = py :: module :: import("torch");

        // py::object displ = torch.attr("zeros")(13);  // shape [13]

        py :: array_t< float >x_array( { xf.size() }, xf.data() );
        py :: object displ = torch.attr("from_numpy")(x_array).attr("clone")();

        std :: cout << "Input: " << py :: str(displ).cast< std :: string >() << std :: endl;

        // Call forward pass (force computation)
        py :: object result = model(displ);
        py :: tuple outputs = result.cast< py :: tuple >();

        // First tensor (tangent)
        py :: object tangent = outputs [ 0 ];
        py :: object force   = outputs [ 1 ];

        // Convert tangent to numpy
        py :: object tangent_numpy = tangent.attr("detach")().attr("cpu")().attr("numpy")();
        py :: array_t< float >tangent_array = tangent_numpy.cast< py :: array_t< float > >();
        py :: buffer_info buf_tangent = tangent_array.request();

        // Convert force to numpy
        py :: object force_numpy = force.attr("detach")().attr("cpu")().attr("numpy")();
        py :: array_t< float >force_array = force_numpy.cast< py :: array_t< float > >();
        py :: buffer_info buf_force = force_array.request();

        // Map into Eigen (no copy)
        Eigen :: Map< Eigen :: VectorXf >force_eigen(
            static_cast< float * >( buf_force.ptr ),
            buf_force.shape [ 0 ]
            );

        Eigen :: Map< Eigen :: VectorXf >tangent_eigen(
            static_cast< float * >( buf_tangent.ptr ),
            buf_tangent.shape [ 0 ]
            );

        // Now you can use them
        std :: cout << "Force (Eigen): " << force_eigen.transpose() << std :: endl;
        std :: cout << "Tangent (Eigen): " << tangent_eigen.transpose() << std :: endl;

        // std::cout << "Force: " << py::str(force).cast<std::string>() << std::endl;
    } catch ( const std :: exception &e ) {
        std :: cerr << "Python error: " << e.what() << std :: endl;
    }

    return 0;
}
