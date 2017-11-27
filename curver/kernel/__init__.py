
''' The curver kernel.

Some of the functions and methods have assumptions on them. These are denoted in the docstrings
by "Assumes that ..." meaning that:

    - If the assumptions are met then this function is guaranteed to terminate correctly.
    - If not then a curver.AssumptionError will be raised. '''

from . import arc
from . import crush
from . import curve
from . import curvecomplex
from . import encoding
from . import error
from . import homology
from . import lamination
from . import mappingclassgroup
from . import moves
from . import partition
from . import permutation
from . import traintrack
from . import triangulation
from . import twist
from . import utilities

# Set up shorter names for all of the different classes.
Edge = triangulation.Edge
Triangle = triangulation.Triangle
Triangulation = triangulation.Triangulation
Encoding = encoding.Encoding
MappingClass = encoding.MappingClass
AbortError = error.AbortError
AssumptionError = error.AssumptionError
MappingClassGroup = mappingclassgroup.MappingClassGroup
MCG = MappingClassGroup  # Alias.
HomologyClass = homology.HomologyClass
Lamination = lamination.Lamination
Shortenable = lamination.Shortenable
MultiCurve = curve.MultiCurve
Curve = curve.Curve
MultiArc = arc.MultiArc
Arc = arc.Arc
Move = moves.Move
Isometry = moves.Isometry
EdgeFlip = moves.EdgeFlip
TrainTrack = traintrack.TrainTrack
Crush = crush.Crush
Lift = crush.Lift
CurveComplex = curvecomplex.CurveComplex
Twist = twist.Twist
HalfTwist = twist.HalfTwist
CurvePartitionGraph = partition.CurvePartitionGraph
Permutation = permutation.Permutation

norm = triangulation.norm
UnionFind = utilities.UnionFind
memoize = utilities.memoize  # Unused?

# Functions that help with construction.
create_triangulation = Triangulation.from_tuple
triangulation_from_sig = Triangulation.from_sig

