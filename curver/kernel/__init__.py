
''' The curver kernel.

Some of the functions and methods have assumptions on them. These are denoted in the docstrings
by "Assumes that ..." meaning that:

    - If the assumptions are met then this function is guaranteed to terminate correctly.
    - If not then a curver.AssumptionError will be raised. '''

from . import arc
from . import crush
from . import curve
from . import curvegraph
from . import encoding
from . import error
from . import homology
from . import lamination
from . import mappingclassgroup
from . import moves
from . import partition
from . import permutation
from . import structures
from . import triangulation
from . import twist
from . import utilities

# Set up shorter names for all of the different classes.
Edge = triangulation.Edge
Triangle = triangulation.Triangle
Triangulation = triangulation.Triangulation
Encoding = encoding.Encoding
Mapping = encoding.Mapping
MappingClass = encoding.MappingClass
AbortError = error.AbortError
AssumptionError = error.AssumptionError
MappingClassGroup = mappingclassgroup.MappingClassGroup
MCG = MappingClassGroup  # Alias.
HomologyClass = homology.HomologyClass
Lamination = lamination.Lamination
MultiCurve = curve.MultiCurve
Curve = curve.Curve
MultiArc = arc.MultiArc
Arc = arc.Arc
Move = moves.Move
FlipGraphMove = moves.FlipGraphMove
Isometry = moves.Isometry
EdgeFlip = moves.EdgeFlip
Crush = crush.Crush
Lift = crush.Lift
CurveGraph = curvegraph.CurveGraph
Twist = twist.Twist
HalfTwist = twist.HalfTwist
CurvePartitionGraph = partition.CurvePartitionGraph
Permutation = permutation.Permutation

norm = triangulation.norm
UnionFind = structures.UnionFind
StraightLineProgram = structures.StraightLineProgram
matrix_vector = utilities.matrix_vector

# Functions that help with construction.
create_triangulation = Triangulation.from_tuple
triangulation_from_sig = Triangulation.from_sig

