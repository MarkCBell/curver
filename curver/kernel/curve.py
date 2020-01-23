
''' A module for representing (multi)curves on triangulations. '''

from fractions import Fraction
from collections import Counter
import networkx
import numpy as np

import curver
from curver.kernel.lamination import IntegralLamination  # Special import needed for subclassing.
from curver.kernel.decorators import memoize, topological_invariant  # Special import needed for decorating.

class MultiCurve(IntegralLamination):
    ''' An IntegralLamination in which every component is a Curve. '''
    def boundary(self):
        return 2*self
    def is_filling(self):
        return False
    
    def encode_twist(self, power=1):
        ''' Return an Encoding of a right Dehn (multi)twist about the components of this multicurve, raised to the given power. '''
        
        if self.is_peripheral() or power == 0:  # Boring case.
            return self.triangulation.id_encoding()
        
        short, conjugator = self.shorten()
        
        h = conjugator
        for curve, multiplicity in short.components().items():
            if not curve.is_peripheral():
                h = curver.kernel.create.twist(curve, power * multiplicity).encode() * h
        h = conjugator.inverse() * h
        
        return h
    
    def vertex_cycles(self):
        ''' Yield the vertex cycles of this multicurve.
        
        These are the curves that use the same normal arcs as this multicurve but only pass through each edge at most twice.
        Be careful as there are often a *lot* of them. '''
        
        def connected_to(edge):
            ''' Yield the edges you can reach by travelling out of the given edge. '''
            corner = self.triangulation.corner_lookup[edge]
            if self.left_weight(edge): yield ~corner[2]  # Can turn left.
            if self.right_weight(edge): yield ~corner[1]  # Can turn right.
        
        # Build graph.
        edges = [(edge, edgy) for edge in self.triangulation.edges for edgy in connected_to(edge)]
        G = networkx.DiGraph(edges)
        
        for cycle in networkx.simple_cycles(G):
            curve = self.triangulation.lamination_from_cut_sequence(cycle)
            if isinstance(curve, curver.kernel.Curve):
                yield curve
    
    def vertex_cycle(self):
        ''' Return a vertex cycle of this multicurve. '''
        
        return next(iter(self.vertex_cycles()))
    
    def crush(self):
        ''' Return the crush map associated to this MultiCurve. '''
        
        short, conjugator = self.shorten()
        
        crush = short.triangulation.id_encoding()
        for curve in short.components():
            if not curve.is_peripheral():
                next_crush = crush(curve).crush()  # Map forward under crushes first.
                crush = next_crush * crush
        
        _, post_conjugator = crush(conjugator(self.triangulation.as_lamination())).shorten()
        
        return post_conjugator * crush * conjugator
    
    def boundary_union(self, other):
        ''' Return \\partial N(self \\cup other). '''
        assert isinstance(other, IntegralLamination)
        
        crush = self.crush()
        lift = crush.inverse()
        other_prime = crush(other)
        m_prime = other_prime.boundary()
        return lift(m_prime)  # = m.
    
    def fills_with(self, other):
        ''' Return whether self \\cup other fills. '''
        assert isinstance(other, IntegralLamination)
        
        if any(component.intersection(other) == 0 for component in self.components()):
            return False
        
        crush = self.crush()
        other_prime = crush(other)
        return other_prime.is_filling()
    
    @topological_invariant
    def is_separating(self):
        ''' Return whether this multicurve separates S.
        
        That is, whether S - self has more components than S. '''
        
        crush = self.crush()
        return len(crush.target_triangulation.components()) > len(self.triangulation.components())

class Curve(MultiCurve):
    ''' A MultiCurve with a single component. '''
    @memoize
    def components(self):
        return {self: 1}
    
    def parallel(self):
        ''' Return an edge that this curve is parallel to.
        
        Note that this is only defined for short, non-peripheral curves. '''
        
        assert not self.is_peripheral()
        assert self.is_short()
        
        [(component, (multiplicity, edge))] = self.parallel_components().items()
        assert component == self  # Sanity.
        assert multiplicity == 1  # Sanity.
        
        return edge
    
    def is_short(self):
        return self.is_peripheral() or len(self.parallel_components()) == 1
    
    def encode_twist(self, power=1):
        ''' Return an Encoding of a right Dehn twist about this curve, raised to the given power. '''
        
        if self.is_peripheral() or power == 0:  # Boring case.
            return self.triangulation.id_encoding()
        
        short, conjugator = self.shorten()
        
        return conjugator.inverse() * curver.kernel.create.twist(short, power).encode() * conjugator
    
    def slope(self, lamination):
        ''' Return the slope of the given lamination about this curve.
        
        This is a Fraction that increases by one each time a right Dehn twist about
        this curve is performed unless -1 <= slope <= 1.
        
        This curve must be non-peripheral and intersect the given lamination. '''
        
        assert isinstance(lamination, curver.kernel.Lamination)
        
        if self.is_peripheral():
            raise ValueError('Curve is peripheral')
        
        short, conjugator = self.shorten()
        short_lamination = conjugator(lamination)
        
        # Get some edges.
        a = short.parallel()
        v = short.triangulation.vertex_lookup[a]  # = short.triangulation.vertex_lookup[~a].
        
        v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
        around_v = curver.kernel.utilities.minimal((short_lamination.left_weight(edgy) for edgy in v_edges), lower_bound=0)
        out_v = sum(max(-short_lamination.left_weight(edge), 0) for edge in v_edges) + sum(max(-short_lamination(edge), 0) for edge in v_edges[1:])
        
        denominator = max(short_lamination(a), 0) - 2 * around_v + out_v  # = short.intersection(short_lamination)
        if denominator == 0:
            raise ValueError('Slope is undefined when self is disjoint from lamination')
        
        twisting = curver.kernel.utilities.minimal((short_lamination.left_weight(edgy) - around_v for edgy in v_edges[1:-1]), lower_bound=0)
        
        numerator = twisting
        
        sign = -1 if short_lamination.left_weight(a) - around_v > 0 or short_lamination.right_weight(a) < 0 else +1
        
        return Fraction(sign * numerator, denominator)  # + (1 if sign < 0 and not short.is_isolating() else 0)  # Curver is right biased on non-isolating curves.
    
    def relative_twisting(self, b, c):
        ''' Return the relative twisting number of b about self relative to c.
        
        This is the number of (right) Dehn twists about self that must be applied to b in order to minimise its intersection with c.
        
        This curve must be non-peripheral and intersect the given laminations. '''
        
        assert isinstance(b, curver.kernel.IntegralLamination)
        assert isinstance(c, curver.kernel.IntegralLamination)
        if self.is_peripheral():
            raise ValueError('Curve is peripheral')
        
        ab = self.intersection(b)
        if ab == 0: raise ValueError('Relative slope is undefined when self and b are disjoint')
        ac = self.intersection(c)  # Faster than c.intersection(a) since we know a is a curve.
        if ac == 0: raise ValueError('Relative slope is undefined when self and c are disjoint')
        bc = b.intersection(c)
        
        f_lower = self.encode_twist(power=-2*bc)(b).intersection(c)  # f(-2*bc).
        f_upper = self.encode_twist(power=2*bc)(b).intersection(c)  # f(2*bc).
        
        return Fraction(f_lower - f_upper, 2*ab*ac)  # No division by zero thanks to our previous checks.
    
    def crush(self):
        ''' Return the crush map associated to this Curve. '''
        
        if self.is_peripheral():  # Boring case.
            return self.triangulation.id_encoding()
        
        short, conjugator = self.shorten()
        
        # Use the following for reference:
        #             #<----------#                #  #-----------#  #
        #            /|     a    ^|               /|  |     a    /  /|
        #           / |         / |              / |  |         /  / |
        #          /  |        /  |             /  |  |        /  /  |
        #         /   |       /   |            /   |  |       /  /   |
        #        /    |b    e/    |   ===>>   /    |  |b   ~b/  /    |
        #       /   ~b|     /~e   |          /    e|  |     /  /~e   |
        #      /      |    /      |         /      |  |    /  /      |
        #     /=======|===/=======|        /       |  |   /  /       |
        #    /        |  /        |       /        |  |  /  /        |
        #   /         | /         |      /         |  | /  /         |
        #  /          V/          |     /          |  |/  /          |
        # #-----------#-----------#    #-----------#  #  #-----------#
        # Where === is the short curve and a is its parallel edge.
        
        a = short.parallel()
        a, b, e = short.triangulation.corner_lookup[a]
        
        # Build the new triangulation.
        edge_map = dict((edge, curver.kernel.Edge(edge.label)) for edge in short.triangulation.edges)
        # Remap some edges.
        edge_map[e] = curver.kernel.Edge(~b.label)
        edge_map[~b] = curver.kernel.Edge(e.label)
        
        new_triangulation = curver.kernel.Triangulation([curver.kernel.Triangle([edge_map[edgy] for edgy in triangle]) for triangle in short.triangulation])
        
        # Build the lifting matrix back.
        v = short.triangulation.vertex_lookup[a]  # = short.triangulation.vertex_lookup[~a].
        indices = Counter([edge.index for edge in curver.kernel.utilities.cyclic_slice(v, a, ~a)[1:]])  # The indices that appear walking around v from a to ~a. Note need to exclude the initial a.
        matrix = np.array([[indices[j] if i == b.index else 1 if (i == e.index and j == b.index) else 1 if i == j else 0 for i in range(self.zeta)] for j in range(self.zeta)], dtype=object)
        
        crush = curver.kernel.create.crush(short.triangulation, new_triangulation, short, matrix).encode()
        
        _, post_conjugator = crush(conjugator(self.triangulation.as_lamination())).shorten()
        
        return post_conjugator * crush * conjugator

