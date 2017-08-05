
''' Curver is a program for computations in the curve complex. '''

import pkg_resources
import curver.kernel
from curver.load import load
from numbers import Integral as IntegerType

__version__ = pkg_resources.require('curver')[0].version

# Set up really short names for the most commonly used classes and functions by users.
create_triangulation = curver.kernel.create_triangulation
norm = curver.kernel.norm

AbortError = curver.kernel.AbortError
AssumptionError = curver.kernel.AssumptionError

