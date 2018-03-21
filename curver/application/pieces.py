
import curver

DEFAULT_OBJECT_COLOUR = 'black'
DEFAULT_VERTEX_COLOUR = 'black'
DEFAULT_EDGE_COLOUR = 'black'
DEFAULT_TRIANGLE_COLOUR = 'gray80'
DEFAULT_CURVE_COLOUR = 'grey40'

ARROW_FRAC = 0.55

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]

def intersection(A, d, B, d2):
    # Find the intersection parameters of A + t d and B + t2 d2
    # Want t & t2 such that:
    #   A.x + t d.x = B.x + t2 d2.x
    #   A.y + t d.y = B.y + t2 d2.y
    # So:
    #   (d.x -d2.x) (t ) = (B.x - A.x)
    #   (d.y -d2.y) (t2)   (B.y - A.y)
    # The inverse of this matrix is:
    #   (-d2.y d2.x)
    #   ( -d.y  d.x) / det
    # where:
    det = d2.cross(d)
    # So:
    t = ((B.x - A.x) * -d2.y + (B.y - A.y) * d2.x) / det
    t2 = ((B.x - A.x) * -d.y + (B.y - A.y) * d.x) / det
    
    return t, t2

def interpolate(A, B, C, r, s):
    # Given points A, B, C and parameters r, s
    # Let X := rB + (1-r)A and
    # Y := sB + (1-s)C
    d = A - B
    d2 = C - B
    
    X = B + r*d
    Y = B + s*d2
    
    centroid = (A + B + C) / 3
    d1a = intersection(X, d.rotate(90), B, centroid - B)[0]
    d1b = intersection(X, d.rotate(90), A, centroid - A)[0]
    t = min([x for x in [d1a, d1b] if x > 0]) / 2
    
    d2a = intersection(Y, d2.rotate(-90), B, centroid - B)[0]
    d2b = intersection(Y, d2.rotate(-90), C, centroid - C)[0]
    t2 = min([x for x in [d2a, d2b] if x > 0]) / 2
    
    P = X + t * d.rotate(90)
    Q = Y + t2 * d2.rotate(-90)
    
    return X, P, Q, Y

class DrawableObject(object):
    def __init__(self, canvas, vertices, options):
        self.options = options
        self.canvas = canvas
        self.vertices = vertices
        self.colour = DEFAULT_OBJECT_COLOUR
        self.drawn = None
    
    def __repr__(self):
        return str(self)
    # Note that this means that CanvasTriangle will NOT have the same convention as AbstractTriangle,
    # there iterating and index accesses return edges.
    def __getitem__(self, index):
        return self.vertices[index % len(self)]
    
    def __iter__(self):
        return iter(self.vertices)
    
    def __len__(self):
        return len(self.vertices)
    
    def set_colour(self, colour=None):
        self.colour = colour
        self.canvas.itemconfig(self.drawn, fill=self.colour)
    
    def centre(self):
        return sum([v.vector for v in self.vertices], curver.application.Vector2(0, 0)) * (1.0 / len(self))
    
    def update(self):
        self.canvas.coords(self.drawn, *[c for v in self for c in v])

class CanvasVertex(DrawableObject):
    def __init__(self, canvas, vector, options):
        super(CanvasVertex, self).__init__(canvas, [self], options)
        self.colour = DEFAULT_VERTEX_COLOUR
        self.vector = vector
        self.drawn = self.canvas.create_oval(
            [p + scale*self.options.dot_size for scale in [-1, 1] for p in self],
            outline=self.colour, fill=self.colour, tag='oval'
            )
    
    def __str__(self):
        return str(self.vector)
    
    def __sub__(self, other):
        return self.vector - other.vector
    
    # We have to redo these manually.
    def __iter__(self):
        return iter(self.vector)
    
    def __contains__(self, point):
        return all(abs(c - v) < self.options.epsilon for c, v in zip(point, self))
    
    def update(self):
        self.canvas.coords(self.drawn, *[p + scale*self.options.dot_size for scale in [-1, 1] for p in self])

class CanvasEdge(DrawableObject):
    def __init__(self, canvas, vertices, label, colour, options):
        super(CanvasEdge, self).__init__(canvas, vertices, options)
        self.label = label
        self.colour = DEFAULT_EDGE_COLOUR if colour is None else colour
        m = (1-ARROW_FRAC)*self.vertices[0].vector + ARROW_FRAC*self.vertices[1].vector
        self.drawn = [  # Need two lines really so arrows work correctly.
            self.canvas.create_line(
                [c for v in [self.vertices[0], m] for c in v],
                width=self.options.line_size,
                fill=self.colour,
                tags=['line', 'line_start'],
                arrowshape=self.options.arrow_shape
            ),
            self.canvas.create_line(
                [c for v in self.vertices for c in v],
                width=self.options.line_size,
                fill=self.colour,
                tags=['line', 'line_end'],
                arrowshape=self.options.arrow_shape
            )
            ]
        
        self.in_triangles = []
    
    def __str__(self):
        return str(self.vertices)
    
    def hide(self, hide=False):
        for drawn in self.drawn:
            self.canvas.itemconfig(drawn, state='hidden' if hide else 'normal')
    
    def update(self):
        m = (1-ARROW_FRAC)*self.vertices[0].vector + ARROW_FRAC*self.vertices[1].vector
        self.canvas.coords(self.drawn[0], *[c for v in [self.vertices[0], m] for c in v])
        self.canvas.coords(self.drawn[1], *[c for v in self.vertices for c in v])

class CanvasTriangle(DrawableObject):
    def __init__(self, canvas, edges, options):
        super(CanvasTriangle, self).__init__(canvas, list(set(v for e in edges for v in e)), options)
        self.colour = DEFAULT_TRIANGLE_COLOUR
        self.edges = edges
        
        # We reorder the vertices to guarantee that the vertices are cyclically ordered anticlockwise in the plane.
        if (self[1] - self[0]).cross(self[2] - self[0]) > 0: self.vertices = [self[0], self[2], self[1]]
        # Now we reorder the edges such that edges[i] does not meet vertices[i].
        self.edges = [edge for vertex in self for edge in self.edges if vertex not in edge.vertices]
        
        # And check to make sure everyone made it through alive.
        assert len(self.edges) == 3
        assert self[0] != self[1] and self[1] != self[2] and self[2] != self[0]
        assert self.edges[0] != self.edges[1] and self.edges[1] != self.edges[2] and self.edges[2] != self.edges[0]
        
        self.drawn = self.canvas.create_polygon([c for v in self for c in v], fill=self.colour, tag='polygon')
        # Add this triangle to each edge involved.
        for edge in self.edges:
            edge.in_triangles.append(self)
    
    def __str__(self):
        return str(self.vertices)

class CurveComponent(DrawableObject):
    def __init__(self, canvas, vertices, options, thin=True, smooth=False):
        super(CurveComponent, self).__init__(canvas, vertices, options)
        self.colour = DEFAULT_CURVE_COLOUR
        if thin:
            self.drawn = self.canvas.create_line(
                [c for v in self.vertices for c in v],
                width=self.options.line_size,
                fill=self.colour,
                tag='curve',
                smooth=smooth,
                splinesteps=50
                )
        else:
            self.drawn = self.canvas.create_polygon(
                [c for v in self.vertices for c in v],
                fill=self.colour,
                tag='curve',
                outline=self.colour,
                smooth=smooth,
                splinesteps=50
                )

