
''' Curver is a program for computations in the curve complex. '''

from numbers import Integral as IntegerType  # noqa: F401
import warnings

import curver.kernel
from curver.load import load  # noqa: F401

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import pkg_resources  # Suppress 'UserWarning: Module curver was already imported from ...'
    __version__ = pkg_resources.require('curver')[0].version

# Set up really short names for the most commonly used classes and functions by users.
create_triangulation = curver.kernel.create_triangulation
triangulation_from_sig = curver.kernel.triangulation_from_sig

AbortError = curver.kernel.AbortError
AssumptionError = curver.kernel.AssumptionError

def show(*items):
    import curver.application
    curver.application.start(*items)

