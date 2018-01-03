
''' A module for representing (multi)arcs on triangulations. '''

import curver
from curver.kernel.lamination import Shortenable  # Special import needed for subclassing.
from curver.kernel.utilities import memoize  # Special import needed for decorating.

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
        
        used = set(edge for edge in short.triangulation.edges if short(edge) < 0)
        
        # Also use any edge that is in a triangle where two sides are used.
        # Note that this does not change the boundary.
        to_check = [triangle for triangle in short.triangulation if sum(1 for edge in triangle if edge in used) == 2]
        while to_check:
            triangle = to_check.pop()
            for edge in triangle:
                if edge not in used:
                    used.add(edge)
                    used.add(~edge)
                    neighbour = short.triangulation.triangle_lookup[~edge.label]
                    if sum(1 for edge in neighbour if edge in used) == 2:
                        to_check.append(neighbour)
                    break
        
        # Now build each component by walking around the outside of the used edges.
        passed_through = set()
        geometric = [0] * short.zeta
        for edge in short.triangulation.edges:
            corner = short.triangulation.corner_lookup[edge.label]
            if edge not in passed_through and edge not in used and corner[2] in used:  # Start on a good edge.
                while edge not in passed_through:
                    geometric[edge.index] += 1
                    passed_through.add(edge)
                    corner = short.triangulation.corner_lookup[edge.label]
                    if ~corner[2] not in used:
                        edge = ~corner[2]
                    else:
                        edge = ~corner[1]
        
        boundary = short.triangulation(geometric)
        return conjugator.inverse()(boundary)
    
    def is_polygonalisation(self):
        ''' Return if this MultiArc is a polygonalisation, that is, if it cuts the surface into polygons. '''
        short, _ = self.shorten()
        
        avoid = set(index for index in short.triangulation.indices if short(index) < 0)  # All of the edges used.
        dual_tree = short.triangulation.dual_tree(avoid=avoid)
        
        return all(dual_tree[index] or index in avoid for index in short.triangulation.indices)
    
    def is_triangulation(self):
        ''' Return if this MultiArc is a triangulation. '''
        short, _ = self.shorten()
        
        return all(weight == -1 for weight in short)
    
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
        
        return NotImplemented  # TODO: 2) Implement.

class Arc(MultiArc):
    ''' A MultiArc with a single component. '''
    @memoize(fast=True)
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
        
        return conjugator.inverse() * curver.kernel.HalfTwist(short, power).encode() * conjugator
    
    def intersection(self, lamination):
        ''' Return the geometric intersection between self and the given lamination. '''
        
        assert(isinstance(lamination, curver.kernel.Lamination))
        assert(lamination.triangulation == self.triangulation)
        
        short, conjugator = self.shorten()
        short_lamination = conjugator(lamination)
        
        arc = short.parallel()
        
        return max(short_lamination(arc), 0)

