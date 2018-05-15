
''' A module for representing (multi)curves on triangulations. '''

from fractions import Fraction
from collections import Counter, defaultdict
import networkx
import numpy as np

import curver
from curver.kernel.lamination import Lamination  # Special import needed for subclassing.
from curver.kernel.decorators import memoize, topological_invariant, ensure  # Special import needed for decorating.

class MultiCurve(Lamination):
    ''' A Lamination in which every component is a Curve. '''
    def is_multicurve(self):
        return True
    def is_multiarc(self):
        return False
    def boundary(self):
        return 2*self
    def is_filling(self):
        return False
    
    def encode_twist(self, power=1):
        ''' Return an Encoding of a right Dehn (multi)twist about the components of this multicurve, raised to the given power. '''
        
        h = self.triangulation.id_encoding()
        for curve, multiplicity in self.components().items():
            h = curve.encode_twist(power * multiplicity) * h
        
        return h
    
    def vertex_cycles(self):
        ''' Yield the vertex cycles of this multicurve.
        
        These are the curves that use the same normal arcs as this multicurve but only pass through each edge at most twice.
        Be careful as there are often a *lot* of them. '''
        
        def connected_to(edge):
            ''' Yield the edges you can reach by travelling out of the given edge. '''
            corner = self.triangulation.corner_lookup[edge]
            if self.dual_weight(corner[1]): yield ~corner[2]
            if self.dual_weight(corner[2]): yield ~corner[1]
        
        # Build graph.
        edges = [(edge, edgy) for edge in self.triangulation.edges for edgy in connected_to(edge)]
        G = networkx.DiGraph(edges)
        
        for cycle in networkx.simple_cycles(G):
            curve = self.triangulation.lamination_from_cut_sequence(cycle)
            if isinstance(curve, curver.kernel.Curve):
                yield curve
        
        return
    
    def vertex_cycle(self):
        ''' Return a vertex cycle of this multicurve. '''
        
        return next(iter(self.vertex_cycles()))
    
    def crush(self):
        ''' Return the crush map associated to this MultiCurve. '''
        
        g = self.triangulation.id_encoding()
        for curve in self.components():
            h = g(curve).crush()  # Map forward under crushes first.
            g = h * g
        
        return g
    
    def boundary_union(self, other):
        ''' Return \\partial N(self \\cup other). '''
        assert isinstance(other, Lamination)
        
        crush = self.crush()
        lift = crush.inverse()
        other_prime = crush(other)
        m_prime = other_prime.boundary()
        return lift(m_prime)  # = m.
    
    def fills_with(self, other):
        ''' Return whether self \\cup other fills. '''
        assert isinstance(other, Lamination)
        
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
    
    @topological_invariant
    def topological_type(self):
        ''' Return the topological type of this multicurve.
        
        Two multicurves are in the same mapping class group orbit if and only their topological types are equal.
        These are labelled graphs and so equal means 'label isomorphic', so we return a CurvePartitionGraph class that uses networkx.is_isomorphic to determine equality. '''
        
        # We build the CurvePartitionGraph as follows:
        #  1) Crush along all curve components to obtain S'.
        #  2) Create a graph with a vertex for each component of S'.
        #  3) Label each vertex with the topological type of its component.
        #  4) Connect two vertices with an edge for each curve that you crushed along.
        #  5) Label each edge with the multiplicity of the corresponding curve.
        # Then two multicurves are in the same mapping class group orbit iff there is a label-preserving isomorphism between their graphs.
        
        components = self.components()  # The components of this multicurves.
        crush = self.crush()
        lift = crush.inverse()
        triangulation = crush.target_triangulation
        
        graph = networkx.MultiGraph()
        half_edges = defaultdict(list)
        for index, (component, (g, v)) in enumerate(triangulation.surface().items()):
            graph.add_node(index, genus=g, vertices=v)
            
            for vertex in triangulation.vertices:
                if vertex[0] in component:
                    curve = triangulation.curve_from_cut_sequence(vertex)
                    lifted_curve = lift(curve)
                    if lifted_curve in components:
                        half_edges[lifted_curve].append(index)
        
        dummy_index = len(graph)
        graph.add_node(dummy_index, genus=0, vertices=0)  # Dummy node for peripheral components.
        
        for curve, nodes in half_edges.items():
            if len(nodes) == 2:
                graph.add_edge(nodes[0], nodes[1], weight=components[curve])
            else:  # len(nodes) == 1:
                graph.add_edge(nodes[0], dummy_index, weight=components[curve])
        
        return curver.kernel.CurvePartitionGraph(self, graph)

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
        
        [(component, (multiplicity, edge))] = self.parallel_components().items()  # pylint: disable=unbalanced-tuple-unpacking
        assert component == self  # Sanity.
        assert multiplicity == 1  # Sanity.
        
        return edge
    
    def is_short(self):
        return self.is_peripheral() or len(self.parallel_components()) == 1
    
    def is_minimal(self):
        ''' Return whether this curve is minimal.
        
        A curve is minimal if its weight is as small as possible.
        
        Note that minimal ==> short. '''
        
        if self.is_peripheral():
            return self.weight() == 1
        elif self.is_isolating():
            shapes = [len([edge for edge in triangle if self.dual_weight(edge) > 0]) for triangle in self.triangulation]
            return all(weight == 0 or weight == 2 for weight in self) and shapes.count(1) == 1 and shapes.count(2) == 0  # Exactly one corridor and no bipods.
        else:
            return self.weight() == 2
    
    @ensure(lambda data: data.result(data.self).is_minimal())
    def minimise(self):
        ''' Return an encoding which maps this curve to a minimal one. '''
        
        def minimise_strategy(self, edge):
            ''' Return a float in [0, 1] describing how good flipping this edge is for making this curve minimal. '''
            
            if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
            
            if not self.triangulation.is_flippable(edge): return 0
            
            ai, bi, ci, di, ei = [self(edgy) for edgy in self.triangulation.square(edge)]
            
            if max(ai + ci, bi + di) - ei < ei:  # This flip drops the total weight.
                return 1
            
            return 0
        
        conjugator = self.shorten()
        short = conjugator(self)
        
        # Lemma: There path in the flip graph from a short representative to a minimal one in which the weight strictly decreases at each step.
        
        minimal = short
        while True:
            edge = curver.kernel.utilities.maximum(minimal.triangulation.edges, key=lambda edge: minimise_strategy(minimal, edge), upper_bound=1)
            if minimise_strategy(minimal, edge) == 0: break
            
            move = minimal.triangulation.encode_flip(edge)  # edge is always flippable.
            conjugator = move * conjugator
            minimal = move(minimal)
        
        return conjugator
    
    @topological_invariant
    def is_isolating(self):
        ''' Return if this curve is isolating, that is, if it is non-peripheral and a component of S - self does not contain a puncture. '''
        
        if self.is_peripheral():
            return False
        
        short = self.shorten()(self)
        return all(weight % 2 == 0 for weight in short)
    
    def encode_twist(self, power=1):
        ''' Return an Encoding of a right Dehn twist about this curve, raised to the given power. '''
        
        if self.is_peripheral():  # Boring case.
            return self.triangulation.id_encoding()
        
        conjugator = self.shorten()
        short = conjugator(self)
        
        return conjugator.inverse() * curver.kernel.Twist(short, power).encode() * conjugator
    
    def slope(self, lamination):
        ''' Return the slope of the given lamination about this curve.
        
        This is a Fraction that increases by one each time a right Dehn twist about
        this curve is performed unless -1 <= slope <= 1.
        
        Assumes that this curve and the given lamination intersect.
        This curve must be non-peripheral. '''
        
        assert not self.is_peripheral()
        
        conjugator = self.shorten()
        short = conjugator(self)
        short_lamination = conjugator(lamination)
        
        denominator = short.intersection(short_lamination)
        if denominator == 0:
            raise curver.AssumptionError('Slope is undefined when curves are disjoint.')
        
        # Get some edges.
        a = short.parallel()
        v = short.triangulation.vertex_lookup[a]  # = short.triangulation.vertex_lookup[~a].
        _, b, e = short.triangulation.corner_lookup[a]
        
        v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
        around_v = min(max(short_lamination.side_weight(edge), 0) for edge in v_edges)
        twisting = min(max(short_lamination.side_weight(edge) - around_v, 0) for edge in v_edges[1:-1])
        
        numerator = twisting
        
        sign = -1 if short_lamination.side_weight(a) > around_v or short_lamination.dual_weight(e) < 0 else +1
        
        return Fraction(sign * numerator, denominator) + (1 if sign < 0 and not short.is_isolating() else 0)  # Curver is right biased on non-isolating curves.
    
    def relative_twisting(self, b, c):
        ''' Return the relative twisting number of b about self relative to c.
        
        This is the number of (right) Dehn twists about self that must be applied to b in order to minimise its intersection with c.
        
        Assumes that this curve and the given laminations intersect.
        This curve must be non-peripheral. '''
        
        assert isinstance(b, curver.kernel.Lamination)
        assert isinstance(c, curver.kernel.Lamination)
        assert not self.is_peripheral()
        
        ab = self.intersection(b)
        if ab == 0: raise curver.AssumptionError('Relative slope is undefined when self and b are disjoint.')
        ac = self.intersection(c)  # Faster than c.intersection(a) since we know a is a curve.
        if ac == 0: raise curver.AssumptionError('Relative slope is undefined when self and c are disjoint.')
        bc = b.intersection(c)
        
        f_lower = self.encode_twist(power=-2*bc)(b).intersection(c)  # f(-2*bc).
        f_upper = self.encode_twist(power=2*bc)(b).intersection(c)  # f(2*bc).
        
        return Fraction(f_lower - f_upper, 2*ab*ac)  # No division by zero thanks to our previous checks.
    
    def crush(self):
        ''' Return the crush map associated to this Curve. '''
        
        if self.is_peripheral():  # Boring case.
            return self.triangulation.id_encoding()
        
        conjugator = self.shorten()
        short = conjugator(self)
        
        # Use the following for reference:
        #             #<----------#                #  #-----------#  #
        #            /|     a    ^|               /|  |     a    /  /|
        #           / |         / |              / |  |         /  / |
        #          /  |        /  |             /  |  |        /  /  |
        #         /   |       /   |            /   |  |       /  /   |
        #        /    |b    e/    |   ===>>   /    |  |b   ~b/  /    |
        #       /   ~b|     /~e   |          /    e|  |     /  /~e   |
        #      /      |    /      |         /      |  |    /  /      |
        #     /       |   /       |        /       |  |   /  /       |
        #    /        |  /        |       /        |  |  /  /        |
        #   /         | /         |      /         |  | /  /         |
        #  /          V/          |     /          |  |/  /          |
        # #-----------#-----------#    #-----------#  #  #-----------#
        # Where a is parallel to short.
        
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
        matrix = np.matrix([[indices[j] if i == b.index else 1 if (i == e.index and j == b.index) else 1 if i == j else 0 for i in range(self.zeta)] for j in range(self.zeta)], dtype=object)
        
        crush = curver.kernel.Crush(short.triangulation, new_triangulation, short, matrix).encode()
        return crush * conjugator

