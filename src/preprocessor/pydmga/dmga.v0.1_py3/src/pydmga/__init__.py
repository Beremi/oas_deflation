"""
This is a Python interface to DMGAlpha library. It has more functions and it
is much simpler than the bare C++ counterpart. It uses dmga2py module to communicate
with C++ library dmga. 
"""

from . import container
from . import diagram
from . import geometry
from . import model
from . import shape

__all__ = ["container", "diagram", "geometry", "model", "shape"]
