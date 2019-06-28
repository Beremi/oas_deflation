from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(ext_modules=cythonize(
           "point_generators_cython.pyx",                 # our Cython source
           # sources=["Rectangle.cpp"],  # additional source file(s)
           # language="c++",             # generate C++ code - old Cython
      ),
      include_dirs = [np.get_include()],
)
