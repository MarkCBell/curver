
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

from . import encoding
from . import error
from . import equippedtriangulation
from . import lamination
from . import moves
from . import triangulation

# Set up shorter names for all of the different classes.
Edge = triangulation.Edge
Triangle = triangulation.Triangle
Triangulation = triangulation.Triangulation
Encoding = encoding.Encoding
AbortError = error.AbortError
AssumptionError = error.AssumptionError
EquippedTriangulation = equippedtriangulation.EquippedTriangulation
Lamination = lamination.Lamination
MultiCurve = lamination.MultiCurve
Curve = lamination.Curve
MultiArc = lamination.MultiArc
Arc = lamination.Arc
Move = moves.Move
Isometry = moves.Isometry
EdgeFlip = moves.EdgeFlip
Spiral = moves.Spiral

norm = triangulation.norm

# Functions that help with construction.
create_triangulation = Triangulation.from_tuple

