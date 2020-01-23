
''' A module for representing (multi)arcs on triangulations. '''

from itertools import combinations
import networkx

import curver
from curver.kernel.lamination import IntegralLamination  # Special import needed for subclassing.
from curver.kernel.decorators import memoize, topological_invariant  # Special import needed for decorating.

class MultiArc(IntegralLamination):
    ''' An IntegralLamination in which every component is an Arc. '''
    def is_short(self):
        return all(weight <= 0 for weight in self)
    
    def vertices(self):
        ''' Return set of vertices that the components of this MultiArc connects to. '''
        
        return set(vertex for vertex in self.triangulation.vertices if any(self(edge) < 0 or self.left_weight(edge) < 0 for edge in vertex))
    
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
    def is_triangulation(self):
        ''' Return whether this MultiArc is a triangulation. '''
        short, _ = self.shorten()
        
        return all(weight == -1 for weight in short)
    
    def explore_ball(self, radius):
        ''' Extend this MultiArc to a triangulation and return all triangulations within the ball of the given radius of that one.
        
        Runs in exp(radius) time.
        Note that this is only well-defined if this multiarc is filling. '''
        
        assert self.is_filling()
        
        short, conjugator = self.shorten()
        
        triangulations = set()
        for encoding in short.triangulation.all_encodings(radius):
            T = encoding.target_triangulation.as_lamination()
            triangulations.add(conjugator.inverse()(encoding.inverse()(T)))
        
        return triangulations
    
    def all_disjoint_arcs(self):
        ''' Yield all arcs that are disjoint from this multiarc.
        
        Note that is only well defined if self is filling. '''
        
        assert self.is_filling()
        
        short, conjugator = self.shorten()
        conjugator_inv = conjugator.inverse()
        
        # All arcs that meet 0 edges.
        for arc in short.triangulation.edge_arcs():
            yield conjugator_inv(arc)
        
        # All arcs that meet 1 edge.
        for index in short.triangulation.indices:
            if short(index) == 0:
                arc = short.triangulation([1 if i == index else 0 for i in range(short.triangulation.zeta)])
                yield conjugator_inv(arc)
        
        # All arcs that meet at least two edges.
        edges = [(short.triangulation.triangle_lookup[edge], short.triangulation.triangle_lookup[~edge]) for edge in short.triangulation.positive_edges if short(edge) == 0]
        G = networkx.Graph(edges)
        for t1, t2 in combinations(G.nodes(), r=2):
            paths = networkx.all_simple_paths(G, t1, t2)
            for path in paths:
                if len(path) > 2:  # Not covered by above case.
                    # Convert back to a sequence of edges.
                    cut_sequence = [edge.index for p1, p2 in zip(path, path[1:]) for edge in p1 if ~edge in p2]
                    arc = short.triangulation.lamination_from_cut_sequence(cut_sequence)
                    yield conjugator_inv(arc)
    
    def all_disjoint_multiarcs(self):
        ''' Yield all multiarcs that are disjoint from this one.
        
        Assumes that this multiarc is filling. '''
        
        arcs = list(self.all_disjoint_arcs())  # Checks is filling.
        
        G = networkx.Graph()
        G.add_nodes_from(arcs)
        G.add_edges_from([(a_1, a_2) for a_1, a_2 in combinations(arcs, r=2) if a_1.intersection(a_2) == 0])
        
        for clique in networkx.enumerate_all_cliques(G):
            yield self.triangulation.disjoint_sum(clique)
    
    def all_disjoint_triangulations(self):
        ''' Yield all multiarcs that are triangulations that are disjoint from self.
        
        Note that these must all therefore contain this multiarc.
        
        Assumes that this multiarc is filling. '''
        
        arcs = list(self.all_disjoint_arcs())  # Checks is filling.
        
        G = networkx.Graph()
        G.add_nodes_from(arcs)
        G.add_edges_from([(a_1, a_2) for a_1, a_2 in combinations(arcs, r=2) if a_1.intersection(a_2) == 0])
        
        for clique in networkx.find_cliques(G):
            yield self.triangulation.disjoint_sum(clique)

class Arc(MultiArc):
    ''' A MultiArc with a single component. '''
    @memoize
    def components(self):
        return {self: 1}
    
    def parallel(self):
        ''' Return an edge that this arc is parallel to.
        
        Note that this is only defined for short arcs. '''
        
        assert self.is_short()
        
        [(component, (multiplicity, edge))] = self.parallel_components().items()
        assert component == self  # Sanity.
        assert multiplicity == 1  # Sanity.
        
        return edge
    
    @topological_invariant
    def connects_distinct_vertices(self):
        ''' Return whether this arc connects between distict vertices of its underlying triangulation. '''
        
        return len(self.vertices()) == 2
    
    def encode_halftwist(self, power=1):
        ''' Return an Encoding of a right half twist about a regular neighbourhood of this arc, raised to the given power.
        
        This arc must connects between distinct vertices. '''
        
        if not self.connects_distinct_vertices():  # Check where it connects.
            raise ValueError('Arc connects a vertex to itself')
        
        if power == 0:  # Boring case.
            return self.triangulation.id_encoding()
        
        short, conjugator = self.shorten()
        
        return conjugator.inverse() * curver.kernel.create.halftwist(short, power).encode() * conjugator

