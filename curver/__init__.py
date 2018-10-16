
''' Curver is a program for computations in the curve complex. '''

from numbers import Integral as IntegerType  # noqa: F401

import curver.kernel
from curver.load import load  # noqa: F401

import pkg_resources
__version__ = pkg_resources.get_distribution('curver').version

# Set up really short names for the most commonly used classes and functions by users.
create_triangulation = curver.kernel.create_triangulation
triangulation_from_sig = curver.kernel.triangulation_from_sig

def show(*items):
    import curver.application
    curver.application.start(*items)

