
''' Curver is a program for computations in the curve complex. '''

from curver.version import __version__

import curver.kernel

from curver.load import load

from numbers import Integral as IntegerType

# Set up really short names for the most commonly used classes and functions by users.
create_triangulation = curver.kernel.create_triangulation
norm = curver.kernel.norm

AbortError = curver.kernel.AbortError
AssumptionError = curver.kernel.AssumptionError

