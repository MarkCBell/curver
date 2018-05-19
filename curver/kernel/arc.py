
''' A module for representing (multi)arcs on triangulations. '''

import curver
from curver.kernel.lamination import Lamination  # Special import needed for subclassing.
from curver.kernel.decorators import memoize, topological_invariant, ensure  # Special import needed for decorating.

class MultiArc(Lamination):
    ''' A Lamination in which every component is an Arc. '''
    def is_multicurve(self):
        return False
    def is_multiarc(self):
        return True
    def is_short(self):
        return all(weight <= 0 for weight in self)
    
    def vertices(self):
        ''' Return set of vertices that the components of this MultiArc connects to. '''
        
        return set(vertex for vertex in self.triangulation.vertices if any(self(edge) < 0 or self.side_weight(edge) < 0 for edge in vertex))
    
    def boundary(self):
        ''' Return the multicurve which is the boundary of a regular neighbourhood of this multiarc. '''
        
        conjugator = self.shorten()
        short = conjugator(self)
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
                    neighbour = short.triangulation.triangle_lookup[~edge]
                    if sum(1 for edge in neighbour if edge in used) == 2:
                        to_check.append(neighbour)
                    break
        
        # Now build each component by walking around the outside of the used edges.
        passed_through = set()
        geometric = [0] * short.zeta
        for edge in short.triangulation.edges:
            corner = short.triangulation.corner_lookup[edge]
            if edge not in passed_through and edge not in used and corner[2] in used:  # Start on a good edge.
                while edge not in passed_through:
                    geometric[edge.index] += 1
                    passed_through.add(edge)
                    corner = short.triangulation.corner_lookup[edge]
                    if ~corner[2] not in used:
                        edge = ~corner[2]
                    else:
                        edge = ~corner[1]
        
        boundary = short.triangulation(geometric)
        return conjugator.inverse()(boundary)
    
    @topological_invariant
    def is_polygonalisation(self):
        ''' Return if this MultiArc is a polygonalisation, that is, if it cuts the surface into polygons. '''
        short = self.shorten()(self)
       
        avoid = set(index for index in short.triangulation.indices if short(index) < 0)  # All of the edges used.
        dual_tree = short.triangulation.dual_tree(avoid=avoid)
        
        return all(dual_tree[index] or index in avoid for index in short.triangulation.indices)
    
    @topological_invariant
    def is_triangulation(self):
        ''' Return if this MultiArc is a triangulation. '''
        short = self.shorten()(self)
        
        return all(weight == -1 for weight in short)
    
    def is_minimal(self):
        ''' Return whether this multiarc is minimal.
        
        A multiarc is minimal if its weight is as small as possible.
        
        Note that minimal ==> short. '''
        
        return all(weight <= 0 for weight in self)
    
    @ensure(lambda data: data.result(data.self).is_minimal())
    def minimise(self):
        ''' Return an encoding which maps this multiarc to a minimal one. '''
        
        return self.shorten()
    
    def explore_ball(self, radius):
        ''' Extend this MultiArc to a triangulation and return all triangulations within the ball of the given radius of that one.
        
        Runs in exp(radius) time.
        Note that this is only well-defined if this multiarc is filling. '''
        
        assert self.is_filling()
        
        conjugator = self.shorten()
        short = conjugator(self)
        
        triangulations = set()
        for encoding in short.triangulation.all_encodings(radius):
            T = encoding.target_triangulation.as_lamination()
            triangulations.add(conjugator.inverse()(encoding.inverse()(T)))
        
        return triangulations
    
    # @topological_invariant
    def topological_type(self):
        ''' Return the topological type of this multiarc.
        
        Two multiarcs are in the same mapping class group orbit if and only their topological types are equal.
        These are labelled graphs and so equal means 'label isomorphic', so we return a custom class that uses networkx.is_isomorphic to determine equality. '''
        
        return NotImplemented  # TODO: 2) Implement.

class Arc(MultiArc):
    ''' A MultiArc with a single component. '''
    @memoize
    def components(self):
        return {self: 1}
    
    def parallel(self):
        ''' Return an edge that this arc is parallel to.
        
        Note that this is only defined for short arcs. '''
        
        assert self.is_short()
        
        [(component, (multiplicity, edge))] = self.parallel_components().items()  # pylint: disable=unbalanced-tuple-unpacking
        assert component == self  # Sanity.
        assert multiplicity == 1  # Sanity.
        
        return edge
    
    @topological_invariant
    def connects_distinct_vertices(self):
        ''' Return whether this arc connects between distict vertices of its underlying triangulation. '''
        
        return len(self.vertices()) == 2
    
    def encode_halftwist(self, power=1):
        ''' Return an Encoding of a right half twist about a regular neighbourhood of this arc, raised to the given power.
        
        Assumes that this arc connects between distinct vertices. '''
        
        if not self.connects_distinct_vertices():  # Check where it connects.
            raise curver.AssumptionError('Arc connects a vertex to itself.')
        
        conjugator = self.shorten()
        short = conjugator(self)
        
        return conjugator.inverse() * curver.kernel.HalfTwist(short, power).encode() * conjugator

