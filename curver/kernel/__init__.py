
''' The curver kernel.

Some of the functions and methods have assumptions on them. These are denoted in the docstrings
by "Assumes that ..." meaning that:

    - If the assumptions are met then this function is guaranteed to terminate correctly.
    - If not then a curver.AssumptionError will be raised. '''

# flake8: noqa  # Can't check this file since the imports are never used.
from .arc import Arc, MultiArc
from .crush import Crush, Lift
from .curve import Curve, MultiCurve
from .curvegraph import CurveGraph
from .encoding import Encoding, Mapping, MappingClass
from .error import AssumptionError, AbortError
from .homology import HomologyClass
from .lamination import Lamination
from .mappingclassgroup import MappingClassGroup
from .moves import Move, FlipGraphMove, EdgeFlip, Isometry
from .partition import CurvePartitionGraph
from .permutation import Permutation
from .structures import UnionFind, StraightLineProgram
from .triangulation import Edge, Triangle, Triangulation, norm
from .twist import Twist, HalfTwist
from .utilities import matrix_vector

# Set up shorter names for all of the different classes.
MCG = MappingClassGroup  # Alias.

# Functions that help with construction.
create_triangulation = Triangulation.from_tuple
triangulation_from_sig = Triangulation.from_sig

