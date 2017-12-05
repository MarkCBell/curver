
''' A module for representing (multi)arcs on triangulations. '''

import curver
from curver.kernel.lamination import Shortenable  # Special import needed for subclassing.

class MultiArc(Shortenable):
    ''' A Lamination in which every component is an Arc. '''
    def is_multicurve(self):
        return False
    def is_multiarc(self):
        return True
    
    def shorten_strategy(self, edge):
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        if not self.triangulation.is_flippable(edge):
            return 0
        if self.dual_weight(edge) < 0:
            return 1
        
        return 0
    
    def boundary(self):
        ''' Return the multicurve which is the boundary of a regular neighbourhood of this multiarc. '''
        
        short, conjugator = self.shorten()
        # short is a subset of the edges of the triangulation it is defined on.
        # So its geometric vector is non-positive.
        
        used = [weight < 0 for weight in short]
        
        to_check = [triangle for triangle in short.triangulation if sum(1 for index in triangle.indices if used[index]) == 2]
        while to_check:
            triangle = to_check.pop()
            if sum(1 for index in triangle.indices if used[index]) == 2:
                for edge in triangle:
                    if not used[edge.index]:
                        used[edge.index] = True
                        to_check.append(short.triangulation.triangle_lookup[~edge.label])
        
        marked = set()
        components = []
        for edge in short.triangulation.edges:
            corner = short.triangulation.corner_lookup[edge.label]
            if not used[edge.index] and used[corner[2].index]:  # Start on a good edge.  and edge not not in marked:
                start_edge = edge
                geometric = [0] * short.zeta
                while edge not in marked:
                    marked.add(edge)
                    geometric[edge.index] += 1
                    corner = short.triangulation.corner_lookup[edge.label]
                    if not used[corner[2].index]:
                        edge = ~corner[2]
                    else:
                        edge = ~corner[1]
                components.append(curver.kernel.Curve(short.triangulation, geometric))
        
        boundary = short.triangulation.disjoint_sum(components)
        return conjugator.inverse()(boundary)
    
    def is_polygonalisation(self):
        ''' Return if this MultiArc is a polygonalisation, that is, if it cuts the surface into polygons. '''
        short, _ = self.shorten()
        
        avoid = set(index for index in short.triangulation.indices if short(index) < 0)  # All of the edges used.
        dual_tree = short.triangulation.dual_tree(avoid=avoid)
        
        return all(dual_tree[index] or index in avoid for index in short.triangulation.indices)
    
    def explore_ball(self, radius):
        ''' Extend this MultiArc to a triangulation and return all triangulations within the ball of the given radius of that one.
        
        Runs in exp(radius) time.
        Note that this is only well-defined if this multiarc is filling. '''
        
        assert(self.is_filling())
        
        short, conjugator = self.shorten()
        
        triangulations = set()
        for encoding in short.triangulation.all_encodings(radius):
            T = encoding.target_triangulation.as_lamination()
            triangulations.add(conjugator.inverse()(encoding.inverse()(T)))
        
        return triangulations
    
    def topological_type(self):
        ''' Return the topological type of this multiarc.
        
        Two multiarcs are in the same mapping class group orbit if and only their topological types are equal.
        These are labelled graphs and so equal means 'label isomorphic', so we return a custom class that uses networkx.is_isomorphic to determine equality. '''
        
        return NotImplemented  # TODO: 3) Implement.

class Arc(MultiArc):
    ''' A MultiArc with a single component. '''
    def components(self):
        return {self: 1}
    
    def parallel(self):
        ''' Return an edge that this arc is parallel to. '''
        assert(self.is_short())
        
        return min([edge for edge in self.triangulation.edges if self(edge) < 0], key=lambda e: e.label)
    
    def vertices(self):
        ''' Return the pair of vertices that this arc connects from / to. '''
        
        vertices = []
        for vertex in self.triangulation.vertices:
            for edge in vertex:
                if self(edge) < 0: vertices.append(vertex)
                if self.side_weight(edge) == -1: vertices.append(vertex)
                if self.side_weight(edge) == -2: vertices += [vertex, vertex]
        assert(len(vertices) == 2)
        return vertices
    
    def connects_distinct_vertices(self):
        ''' Return whether this arc connects between distict vertices of its underlying triangulation. '''
        
        return len(set(self.vertices())) == 2
    
    def encode_halftwist(self, power=1):
        ''' Return an Encoding of a right half twist about a regular neighbourhood of this arc, raised to the given power.
        
        Assumes that this arc connects between distinct vertices. '''
        
        if not self.connects_distinct_vertices():  # Check where it connects.
            raise curver.AssumptionError('Arc connects a vertex to itself.')
        
        short, conjugator = self.shorten()
        
        edge = short.parallel()
        
        return conjugator.inverse() * curver.kernel.HalfTwist(short, power).encode() * conjugator
    
    def intersection(self, lamination):
        ''' Return the geometric intersection between self and the given lamination. '''
        
        assert(isinstance(lamination, curver.kernel.Lamination))
        assert(lamination.triangulation == self.triangulation)
        
        short, conjugator = self.shorten()
        short_lamination = conjugator(lamination)
        
        arc = short.parallel()
        
        return max(short_lamination(arc), 0)

