
''' The curver kernel.

Some of the functions and methods have assumptions on them. We denote
in a functions docstring:
	Assumes that ...
		If the assumptions are met then this function is guaranteed to terminate correctly.
		If not then this function will either:
			terminate correctly, OR
			a curver.AssumptionError will be raised.
	
	Assumes (and checks) that ...
		If the assumptions are met then this function is guaranteed to terminate correctly.
		If not then this a curver.AssumptionError will be raised. '''

from . import arc
from . import curve
from . import curvecomplex
from . import crush
from . import encoding
from . import error
from . import equippedtriangulation
from . import homology
from . import lamination
from . import moves
from . import triangulation
from . import traintrack
from . import utilities

# Set up shorter names for all of the different classes.
Edge = triangulation.Edge
Triangle = triangulation.Triangle
Triangulation = triangulation.Triangulation
Encoding = encoding.Encoding
AbortError = error.AbortError
AssumptionError = error.AssumptionError
EquippedTriangulation = equippedtriangulation.EquippedTriangulation
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
Spiral = moves.Spiral
TrainTrack = traintrack.TrainTrack
Crush = crush.Crush
Lift = crush.Lift
CurveComple = curvecomplex.CurveComplex

norm = triangulation.norm
UnionFind = utilities.UnionFind
memoize = utilities.memoize  # Unused?

# Functions that help with construction.
create_triangulation = Triangulation.from_tuple

