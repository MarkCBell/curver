
''' The curver GUI application. '''

from . import main
from . import pieces
from . import options

# Set up shorter names for all of the different classes and some common constructors.
start = main.start
Options = options.Options
CurverApplication = main.CurverApplication
CanvasVertex = pieces.CanvasVertex
CanvasEdge = pieces.CanvasEdge
CanvasTriangle = pieces.CanvasTriangle
CurveComponent = pieces.CurveComponent

lines_intersect = pieces.lines_intersect
interpolate = pieces.interpolate

