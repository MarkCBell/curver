
''' A module for representing laminations on Triangulations. '''

from itertools import permutations, groupby, product, chain
from six.moves.queue import Queue

import curver
from curver.kernel.decorators import memoize, topological_invariant, ensure  # Special import needed for decorating.

class Lamination(object):
    ''' This represents a lamination on a triangulation.
    
    Users should create these via Triangulation(...) or Triangulation.lamination(...). '''
    def __init__(self, triangulation, geometric):
        assert isinstance(triangulation, curver.kernel.Triangulation)
        
        self.triangulation = triangulation
        self.zeta = self.triangulation.zeta
        self.geometric = geometric
        
        # Store some additional weights that are often used.
        self._dual = dict()
        self._left = dict()
        self._right = dict()
        for triangle in self.triangulation:
            i, j, k = triangle  # Edges.
            a, b, c = self.geometric[i.index], self.geometric[j.index], self.geometric[k.index]
            af, bf, cf = max(a, 0), max(b, 0), max(c, 0)  # Correct for negatives.
            correction = min(af + bf - cf, bf + cf - af, cf + af - bf, 0)
            self._dual[i] = self._right[j] = self._left[k] = curver.kernel.utilities.half(bf + cf - af + correction)
            self._dual[j] = self._right[k] = self._left[i] = curver.kernel.utilities.half(cf + af - bf + correction)
            self._dual[k] = self._right[i] = self._left[j] = curver.kernel.utilities.half(af + bf - cf + correction)
    
    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.triangulation, self.geometric)
    def __str__(self):
        return '%s %s on %s' % (self.__class__.__name__, '[' + ', '.join(str(weight) for weight in self.geometric) + ']', self.triangulation)
    def __iter__(self):
        return iter(self.geometric)
    def __call__(self, edge):
        ''' Return the geometric measure assigned to item. '''
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        return self.geometric[edge.index]
    def __bool__(self):
        return not self.is_empty()
    def __nonzero__(self):  # For Python2.
        return self.__bool__()
    def __eq__(self, other):
        if not isinstance(other, Lamination): return False
        return self.triangulation == other.triangulation and self.geometric == other.geometric
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(tuple(self.geometric))
    def __add__(self, other):
        # Haken sum.
        if isinstance(other, Lamination):
            if other.triangulation != self.triangulation:
                raise ValueError('Laminations must be on the same triangulation to add them')
            
            geometric = [x + y for x, y in zip(self.geometric, other.geometric)]
            return self.triangulation(geometric)  # Have to promote.
        else:
            return NotImplemented
    def __radd__(self, other):
        return self + other  # Commutative.
    def __mul__(self, other):  # FIXME: Make work for non-integrals.
        assert other >= 0
        
        # TODO: 3) Could save components if they have already been computed.
        geometric = [other * x for x in self]
        
        # In some easy cases we use short-cuts to avoid promote.
        if other == 0:
            return IntegralLamination(self.triangulation, geometric)
        elif other == 1:
            return self.__class__(self.triangulation, geometric)
        elif isinstance(other, curver.IntegerType) and isinstance(self, curver.kernel.MultiArc):  # or Arc.
            return curver.kernel.MultiArc(self.triangulation, geometric)
        elif isinstance(other, curver.IntegerType) and isinstance(self, curver.kernel.MultiCurve):  # or Curve.
            return curver.kernel.MultiCurve(self.triangulation, geometric)
        elif isinstance(other, curver.IntegerType) and isinstance(self, curver.kernel.IntegralLamination):
            return curver.kernel.IntegralLamination(self.triangulation, geometric)
        else:
            return self.triangulation(geometric)  # Have to promote.
    def __rmul__(self, other):
        return self * other  # Commutative.
    
    def __reduce__(self):
        return (self.__class__, (self.triangulation, self.geometric))
    
    def weight(self):
        ''' Return the geometric intersection of this lamination with its underlying triangulation. '''
        
        return sum(max(weight, 0) for weight in self)
    
    def dual_weight(self, edge):
        ''' Return the number of component of this lamination dual to the given edge.
        
        Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
        
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        return self._dual[edge]
    
    def left_weight(self, edge):
        ''' Return the number of component of this lamination dual to the left of the given edge.
        
        Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
        
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        return self._left[edge]
    
    def right_weight(self, edge):
        ''' Return the number of component of this lamination dual to the right the given edge.
        
        Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
        
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        return self._right[edge]
    
    def is_integral(self):
        ''' Return whether this lamination is integral. '''
        
        return all(weight == int(weight) for weight in self) and all(self.dual_weight(edge) == int(self.dual_weight(edge)) for edge in self.triangulation.edges)
    
    def promote(self):
        ''' Return this lamination in its finest form. '''
        
        if not self.is_integral():
            return self
        
        temp = IntegralLamination(self.triangulation, self.geometric)
        short, _ = temp.shorten()  # Shorten returns a short lamination of the correct class.
        promoted = short.__class__(self.triangulation, self.geometric)
        
        # Move cache across.
        try:
            promoted._cache = temp._cache  # pylint: disable=attribute-defined-outside-init
        except AttributeError:
            pass  # No cache.
        
        return promoted
    
    @topological_invariant
    def is_empty(self):
        ''' Return whether this lamination has no components. '''
        
        return not any(self)  # self.num_components() == 0
    
    @topological_invariant
    def is_peripheral(self):
        ''' Return whether this lamination consists entirely of peripheral components. '''
        
        return self.peripheral(promote=False) == self
    
    def peripheral(self, promote=True):  # pylint: disable=unused-argument
        ''' Return the lamination consisting of the peripheral components of this Lamination. '''
        
        return self.triangulation.disjoint_sum(dict((component, multiplicity) for component, (multiplicity, _) in self.peripheral_components().items()))  # Promotes are free!
    
    @topological_invariant
    def is_non_peripheral(self):
        ''' Return whether this lamination does not have any peripheral components. '''
        
        return not self.peripheral(promote=False)
    
    def non_peripheral(self, promote=True):
        ''' Return the lamination consisting of the non-peripheral components of this Lamination. '''
        
        geometric = [x - y for x, y in zip(self, self.peripheral(promote=False))]
        
        return self.triangulation(geometric, promote)  # Have to promote.
    
    def trace_curve(self, edge, intersection, max_length):
        ''' Return the curve obtained by following along this lamination and closing up when you get back to this edge.
        
        We start at the given edge and intersection point and only go for at most max_length.
        A ValueError is raised if:
        
         - we do not get back to the starting edge within this number of steps,
         - the lamination terminates into a vertex, or
         - upon returning to start_edge we cannot close up without creating intersections. '''
        
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        tilde_upper = self(edge) + 1
        tilde_intersection = self(edge) - intersection
        tilde_lower = -1
        
        #   ~upper ~intersection ~lower
        #       V       V          V
        # *<----|||||---|----||-||-|----*
        #               ^
        #         intersection
        
        start_edge = edge
        
        assert 0 <= intersection < self(edge)  # Sanity.
        dual_weights = dict((edge, self.dual_weight(edge)) for edge in self.triangulation.edges)
        trace = [edge]
        for _ in range(max_length):
            x, y, z = self.triangulation.corner_lookup[~edge]
            # Move onto next edge.
            if intersection < dual_weights[z]:  # Turn right.
                edge, intersection = y, intersection  # pylint: disable=self-assigning-variable
            elif dual_weights[x] < 0 and dual_weights[z] <= intersection < dual_weights[z] - dual_weights[x]:  # Terminate.
                raise ValueError('Lamination does not trace to a curve')
            else:  # Turn left.
                edge, intersection = z, self(z) - self(x) + intersection
            
            if edge == start_edge:
                tilde_return = self(edge) - intersection
                if tilde_lower < tilde_return < tilde_upper:
                    return self.triangulation.curve_from_cut_sequence(trace)
                else:
                    raise ValueError('Curve does not close up without intersection')
            if edge == ~start_edge:  # Move the bound in.
                if intersection < tilde_intersection:
                    tilde_lower = max(tilde_lower, intersection)
                elif intersection > tilde_intersection:
                    tilde_upper = min(tilde_upper, intersection)
            
            trace.append(edge)
            assert 0 <= intersection < self(edge)  # Sanity.
        
        raise ValueError('Curve does not close up in {} steps'.format(max_length))
    
    @memoize
    def peripheral_components(self):
        ''' Return a dictionary mapping component to (multiplicity, vertex) for each component of self that is peripheral around a vertex. '''
        
        components = dict()
        for vertex in self.triangulation.vertices:
            multiplicity = curver.kernel.utilities.minimal((self.left_weight(edge) for edge in vertex), lower_bound=0)
            if multiplicity > 0:
                component = self.triangulation.curve_from_cut_sequence(vertex)
                components[component] = (multiplicity, vertex)
        
        return components
    
    @memoize
    def parallel_components(self):
        ''' Return a dictionary mapping component to (multiplicity, edge) for each component of self that is parallel to an edge. '''
        
        components = dict()
        for edge in self.triangulation.edges:
            if edge.sign() == +1:  # Don't double count.
                multiplicity = -self(edge)
                if multiplicity > 0:
                    component = self.triangulation.edge_arc(edge)
                    components[component] = (multiplicity, edge)
            
            if self.triangulation.vertex_lookup[edge] == self.triangulation.vertex_lookup[~edge]:
                v = self.triangulation.vertex_lookup[edge]  # = self.triangulation.vertex_lookup[~edge].
                v_edges = curver.kernel.utilities.cyclic_slice(v, edge, ~edge)  # The set of edges that come out of v from edge round to ~edge.
                if len(v_edges) > 2:
                    around_v = curver.kernel.utilities.minimal((self.left_weight(edgy) for edgy in v_edges), lower_bound=0)
                    twisting = curver.kernel.utilities.minimal((self.left_weight(edgy) - around_v for edgy in v_edges[1:-1]), lower_bound=0)
                    
                    if self.left_weight(v_edges[0]) == around_v and self.left_weight(v_edges[-1]) == around_v:
                        multiplicity = twisting
                        
                        if multiplicity > 0:
                            component = self.triangulation.curve_from_cut_sequence(v_edges[1:])
                            components[component] = (multiplicity, edge)
        
        return components

class IntegralLamination(Lamination):
    ''' This represents a lamination in which all weights are integral. '''
    
    def skeleton(self):
        ''' Return the lamination obtained by collapsing parallel components. '''
        
        return self.triangulation.disjoint_sum(list(self.components()))
    
    def peek_component(self):
        ''' Return one component of this Lamination. '''
        
        return next(iter(self.components()))
    
    def intersection(self, *laminations):
        ''' Return the geometric intersection number between this lamination and the given one(s).
        
        If multiple laminations are given then ``sum(i(self, lamination) for laminations)`` is returned. '''
        
        assert all(isinstance(lamination, Lamination) for lamination in laminations)
        assert all(lamination.triangulation == self.triangulation for lamination in laminations)
        
        short, conjugator = self.shorten()
        short_laminations = [conjugator(lamination) for lamination in laminations]
        
        intersection = 0
        
        # Peripheral components.
        for _, (multiplicity, vertex) in short.peripheral_components().items():
            for lamination in laminations:
                intersection += multiplicity * sum(max(-lamination(edge), 0) + max(-lamination.left_weight(edge), 0) for edge in vertex)
        
        # Parallel components.
        for component, (multiplicity, p) in short.parallel_components().items():
            if isinstance(component, curver.kernel.Arc):
                for short_lamination in short_laminations:
                    intersection += multiplicity * max(short_lamination(p), 0)
            else:  # isinstance(component, curver.kernel.Curve):
                v = short.triangulation.vertex_lookup[p]  # = self.triangulation.vertex_lookup[~p].
                v_edges = curver.kernel.utilities.cyclic_slice(v, p, ~p)  # The set of edges that come out of v from p round to ~p.
                
                for short_lamination in short_laminations:
                    around_v = curver.kernel.utilities.minimal((short_lamination.left_weight(edgy) for edgy in v_edges), lower_bound=0)
                    out_v = sum(max(-short_lamination.left_weight(edge), 0) for edge in v_edges) + sum(max(-short_lamination(edge), 0) for edge in v_edges[1:])
                    # around_v > 0 ==> out_v == 0; out_v > 0 ==> around_v == 0.
                    intersection += multiplicity * (max(short_lamination(p), 0) - 2 * around_v + out_v)
        
        return intersection
    
    def no_common_component(self, lamination):
        ''' Return that self does not share any components with the given IntegralLamination. '''
        
        assert isinstance(lamination, IntegralLamination)
        
        self_components = self.components()
        return not any(component in self_components for component in lamination.components())
    
    @topological_invariant
    def num_components(self):
        ''' Return the total number of components. '''
        
        return sum(self.components().values())
    
    def sublaminations(self):
        ''' Return all sublaminations that appear within self. '''
        
        components = self.components()
        return [self.triangulation.disjoint_sum(sub) for i in range(len(components)) for sub in permutations(components, i+1)]  # Powerset.
    
    def multiarc(self):
        ''' Return the maximal MultiArc contained within this lamination. '''
        
        return self.triangulation.disjoint_sum(dict((component, multiplicity) for component, multiplicity in self.components().items() if isinstance(component, curver.kernel.Arc)))
    
    def multicurve(self):
        ''' Return the maximal MultiCurve contained within this lamination. '''
        
        return self.triangulation.disjoint_sum(dict((component, multiplicity) for component, multiplicity in self.components().items() if isinstance(component, curver.kernel.Curve)))
    
    def boundary(self):
        ''' Return the boundary of a regular neighbourhood of this lamination. '''
        
        if self.is_empty():
            return self
        
        return self.triangulation.disjoint_sum([self.multiarc().boundary(), self.multicurve().boundary()])
    
    @topological_invariant
    def is_filling(self):
        ''' Return whether this IntegralLamination fills the surface, that is, if it intersects all curves on the surface.
        
        Note that this is equivalent to:
            - it meets every non S_{0,3} component of the surface, and
            - its boundary is peripheral.
        
        Furthermore, if any component of this lamination is a non-peripheral curve then it cannot fill. '''
        
        if any(isinstance(component, curver.kernel.Curve) and not component.is_peripheral() for component in self.components()):
            return False
        
        for component in self.triangulation.components():
            V, E = len([vertex for vertex in self.triangulation.vertices if vertex[0] in component]), len(component) // 2
            if (V, E) != (3, 3):  # component != S_{0, 3}:
                if all(self(edge) == 0 for edge in component):
                    return False
        
        return self.boundary().is_peripheral()
    
    @topological_invariant
    def is_polygonalisation(self):
        ''' Return whether this IntegralLamination is a polygonalisation, that is, if it cuts the surface into polygons. '''
        
        if any(isinstance(component, curver.kernel.Curve) for component in self.components()):
            return False
        
        short, _ = self.shorten()
       
        avoid = set(index for index in short.triangulation.indices if short(index) < 0)  # All of the edges used.
        dual_tree = short.triangulation.dual_tree(avoid=avoid)
        
        return dual_tree.union(avoid) == set(short.triangulation.indices)  # Every edge is in avoid or dual_tree.
    
    def fills_with(self, other):  # pylint: disable=no-self-use
        ''' Return whether self \\cup other fills. '''
        assert isinstance(other, IntegralLamination)
        
        return NotImplemented  # TODO: 2) Implement! (And remove the PyLint disable when done.)
    
    # @topological_invariant
    def topological_type(self, closed=False):  # pylint: disable=no-self-use, unused-argument
        ''' Return the topological type of this lamination..
        If the closed flag is set the the object returned records the topological type of the lamination after applying the forgetful map.
        
        Two laminations are in the same mapping class group orbit if and only their topological types are equal. '''
        
        # We build the topological type as follows:
        #  1) Add boundary(multiarcs) to self.
        #  2) Crush along all curve components to obtain S'.
        #  3) Create a graph with a vertex for each component of S'.
        #  4) Label each vertex with the genus of its component.
        #  5) Connect two vertices with an edge for each curve that you crushed along.
        #  6) Label each edge with the multiplicity of the corresponding curve.
        #  7) Let a := crush(self.multiarc()).
        #  8) Label each vertex with the topological type of a restricted to the corresponding component.
        
        # Then two laminations are in the same mapping class group orbit iff there is a label-preserving isomorphism between their graphs.
        
        if closed and self.multiarc():
            raise ValueError('Cannot apply the forgetful map to a lamination with an Arc component')
        
        short, _ = self.shorten()
        
        boundary = short.multiarc().boundary().non_peripheral()
        lamination = short.triangulation.disjoint_sum([short, boundary])
        crush = lamination.multicurve().crush() if lamination.multicurve() else lamination.triangulation.id_encoding()  # Crush along the curve components.
        lift = crush.inverse()
        triangulation = crush.target_triangulation
        image = crush(short.multiarc())
        
        # Build nodes.
        nodes = list(triangulation.components())
        node_labels = dict((component, S.g) for component, S in triangulation.surface().items())
        
        if not closed:  # Add dummy node, represented by an empty tuple whose node_labels is -1.
            dummy = tuple()
            nodes.append(dummy)
            node_labels[dummy] = -1
        
        duplicated_boundary = set(component for component, multiplicity in boundary.components().items() if multiplicity == 2)
        real_nodes = [node for node in nodes if node not in duplicated_boundary]
        
        for component in duplicated_boundary:
            nodes.append(component)
            node_labels[component] = -2
        
        # Sort nodes by node_labels.
        nodes.sort(key=node_labels.get)
        # Determine their labels.
        best_node_labels = [node_labels[node] for node in nodes]  # We know the best node labels right away.
        
        # Useful lookup maps.
        node_lookup = dict((node, index) for index, node in enumerate(nodes))
        edge_node_map = dict((edge, node) for node in real_nodes for edge in node)
        vertex_node_map = dict((vertex, edge_node_map[vertex[0]]) for vertex in triangulation.vertices)
        
        # Write down all the links.
        links = []  # List of (vertex1, vertex2, label).
        vertex_paired_node_map = dict()  # We will build this in a minute.
        short_components = short.components()  # The components of this lamination.
        half_links = dict()
        for vertex in triangulation.vertices:
            curve = lift(triangulation.curve_from_cut_sequence(vertex))
            if curve in duplicated_boundary:
                links.append((vertex_node_map[vertex], curve, short_components.get(curve, 0)))
                vertex_paired_node_map[vertex] = curve
            elif curve not in half_links:
                half_links[curve] = vertex
            else:  # Found the other half.
                vertex2 = half_links[curve]
                links.append((vertex_node_map[vertex], vertex_node_map[vertex2], short_components.get(curve, 0)))
                del half_links[curve]
                vertex_paired_node_map[vertex] = vertex_node_map[vertex2]
                vertex_paired_node_map[vertex2] = vertex_node_map[vertex]
        if not closed:
            for curve, vertex in half_links.items():  # Add hanging edges to dummy.
                links.append((vertex_node_map[vertex], dummy, short_components.get(curve, 0)))
                vertex_paired_node_map[vertex] = dummy
        
        # Build link label matrix.
        link_labels = [[list() for _ in range(len(nodes))] for _ in range(len(nodes))]  # The empty matrix of lists.
        for node1, node2, label in links:
            link_labels[node_lookup[node1]][node_lookup[node2]].append(label)
            if node1 != node2:  # Don't add self loops twice.
                link_labels[node_lookup[node2]][node_lookup[node1]].append(label)
        # Sort all entries of the matrix.
        for row in link_labels:
            for entry in row:
                entry.sort()
        
        # Build the edge -> edge ordering map.
        ordering = dict()
        for vertex in triangulation.vertices:
            edges = [edge for edge in vertex if image(edge)]
            for e1, e2 in zip(edges, edges[1:] + edges[:1]):
                ordering[e1] = ~e2
        # Build the image -> vertex map.
        classes = curver.kernel.UnionFind(triangulation.edges)
        for edge in triangulation.edges:
            if not image(edge):
                classes.union(edge, ~edge)
        for triangle in triangulation:
            classes.union(*triangle)
        classes_lookup = dict((edge, cls) for cls in classes for edge in cls)
        disjoint_vertices = [vertex for vertex in triangulation.vertices if all(not image(e) for e in vertex)]
        image_vertex_map = dict((edge, vertex) for vertex in disjoint_vertices for edge in classes_lookup[vertex[0]])
        
        best_link_labels = None
        best_node_markings = None
        for X in product(*(permutations(g) for k, g in groupby(range(len(nodes)), key=lambda i: node_labels[nodes[i]]))):  # pylint: disable=too-many-nested-blocks
            perm = list(chain(*X))
            
            inverse_perm = [None] * len(perm)
            for index, i in enumerate(perm):
                inverse_perm[i] = index
            permuted_link_labels = [link_labels[i][j] for index, i in enumerate(inverse_perm) for j in inverse_perm[index:]]
            
            if best_link_labels is not None and permuted_link_labels >= best_link_labels:
                continue
            best_link_labels = list(permuted_link_labels)
            
            permuted_nodes = [nodes[i] for i in perm]
            permuted_node_lookup = dict((node, index) for index, node in enumerate(permuted_nodes))
            
            # Build marking.
            node_markings = []
            for node in permuted_nodes:
                starting_edges = [edge for edge in node if image(edge)]
                if starting_edges:  # Skip nodes that do not have any markings.
                    best_node_marking = None
                    for starting_edge in starting_edges:
                        node_marking = []
                        edge_marking = {starting_edge: 0, ~starting_edge: 0}
                        to_do = Queue()
                        to_do.put(starting_edge)
                        to_do.put(~starting_edge)
                        done = set()
                        while not to_do.empty():
                            current = to_do.get()
                            while current not in done:  # Wall around this polygon.
                                if current not in edge_marking:
                                    edge_marking[current] = edge_marking[~current] = len(edge_marking) // 2
                                    to_do.put(~current)
                                node_marking.append((edge_marking[current], image(current)))
                                done.add(current)
                                current = ordering[current]
                            
                            # Mark the break between one cycle and the next.
                            node_marking.append((-1, permuted_node_lookup[vertex_paired_node_map[image_vertex_map[current]]] if current in image_vertex_map else -1))
                        
                        if best_node_marking is not None and node_marking >= best_node_marking:
                            continue
                        best_node_marking = node_marking
                    
                    node_markings.append(best_node_marking)
                else:
                    node_markings.append([])
            
            if best_node_markings is not None and node_markings >= best_node_markings:
                continue
            best_node_markings = node_markings
        
        return best_node_labels, best_link_labels, best_node_markings
    
    @memoize
    def components(self):
        ''' Return a dictionary mapping components to their multiplicities. '''
        
        components = dict()
        
        short, conjugator = self.shorten()
        conjugator_inv = conjugator.inverse()
        
        for component, (multiplicity, _) in short.peripheral_components().items():
            components[conjugator_inv(component)] = multiplicity
        
        for component, (multiplicity, _) in short.parallel_components().items():
            components[conjugator_inv(component)] = multiplicity
        
        return components
    
    def is_short(self):
        ''' Return whether this lamination is short.
        
        A lamination is short if all of its non-peripheral components are parallel to edges of
        the triangulation that it is defined on. This makes computing its components very easy. '''
        
        lamination = self.non_peripheral(promote=False)
        
        geometric = list(lamination)
        for component, (multiplicity, _) in lamination.parallel_components().items():
            geometric = [x - y * multiplicity for x, y in zip(geometric, component)]
        lamination = Lamination(lamination.triangulation, geometric)
        
        return lamination.is_empty()
    
    @memoize
    @ensure(lambda data: data.result[0].is_short())
    def shorten(self, drop=0.1):
        ''' Return a pair (s, h) where:
         * s is a short lamination of the correct class, and
         * s = h(self)
        
        In each round, we do not look for an accelerating Dehn twist if a flip can drop the weight by at least `drop`%.
        So if `drop` == 0.0 then acceleration is never done and this returns the Mosher flip sequence.
        As long as `drop` > 0 this method runs in poly(||self||) time.
        
        The original version of this method was based on [Bell16]_ but now a simpler and more efficient technique is used.
        The argument why this version runs in polynomial time follows that of [EricksonNayyeri13]_. '''
        
        assert 0.0 <= drop <= 1.0
        
        peripheral = self.peripheral()  # This is more efficient than moving every peripheral component individually.
        lamination = self.non_peripheral(promote=False)
        conjugator = self.triangulation.id_encoding()
        
        def shorten_strategy(self, edge):
            ''' Return a float in [0, 1] describing how good flipping this edge is for making this lamination short. '''
            
            if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
            
            if not self.triangulation.is_flippable(edge): return 0
            
            ad, bd, cd, dd, ed = [self.dual_weight(edgy) for edgy in self.triangulation.square(edge)]
            
            if ed < 0:  # Non-parallel arc.
                return 1
            
            if ed == 0 and ad > 0 and bd > 0:  # Bipod.
                return 0.5
            
            return 0
        
        arc_components, curve_components = dict(), dict()
        has_arcs = True
        while True:
            # Subtract.
            geometric = list(lamination)
            for component, (multiplicity, edge) in lamination.parallel_components().items():
                if lamination(edge) <= 0:
                    geometric = [x - y * multiplicity for x, y in zip(geometric, component)]
                    if isinstance(component, curver.kernel.Arc):
                        arc_components[edge] = multiplicity
                    else:  # isinstance(component, curver.kernel.Curve):
                        curve_components[edge] = multiplicity
            
            lamination = IntegralLamination(lamination.triangulation, geometric)
            
            if not lamination: break
            
            # The arcs will be dealt with in the first round and once they are gone, they are gone.
            has_arcs = has_arcs and any(lamination(edge) < 0 or lamination.dual_weight(edge) < 0 for edge in lamination.triangulation.edges)
            extra = []  # High priority edges to check.
            while True:
                # Note that if lamination does not have any arcs then the max value that shorten_strategy can return is 0.5.
                # Also triangulation.edges are listed in increasing order so this process is deterministic.
                edge = curver.kernel.utilities.maximum(
                    extra + lamination.triangulation.edges,
                    key=lambda edge: shorten_strategy(lamination, edge),
                    upper_bound=1 if has_arcs else 0.5)
                if shorten_strategy(lamination, edge) == 0: break  # No non-parallel arcs or bipods.
                
                a, b, c, d, e = lamination.triangulation.square(edge)
                move = lamination.triangulation.encode_flip(edge)  # edge is always flippable.
                # Since looking for and applying a twist is expensive, we will not do it if:
                #  * drop == 0,
                #  * lamination has little weight, or
                #  * flipping drops the weight by at least drop%.
                if drop > 0 and 4 * self.zeta < lamination.weight() and (1 - drop) * lamination.weight() < move(lamination).weight() < lamination.weight():
                    try:
                        curve = lamination.trace_curve(edge, lamination.left_weight(edge), 2*self.zeta)
                        slope = curve.slope(lamination)  # Will raise a ValueError if these are disjoint.
                        if abs(slope) > 2:  # Can accelerate and slope is large enough to be efficient.
                            move = curve.encode_twist(power=-int(slope))  # Round towards zero.
                    except ValueError:
                        extra = [c, d]
                else:
                    extra = [c, d]
                
                conjugator = move * conjugator
                lamination = move(lamination)
                peripheral = move(peripheral)
            
            # Now all arcs should be parallel to edges and there should now be no bipods.
            assert all(lamination.left_weight(edge) >= 0 for edge in lamination.triangulation.edges)
            assert all(sum(1 if lamination.left_weight(edge) > 0 else 0 for edge in triangle) != 2 for triangle in lamination.triangulation)
            
            # This is pretty inefficient.
            sequence = []  # This contains each (oriented) edge at most once and so can never contain more than 2*self.zeta elements.
            used_edges = set()
            for starting_edge in lamination.triangulation.edges:
                # Found a good (unused) starting place.
                if starting_edge in used_edges or lamination.left_weight(starting_edge) <= 0 or lamination.right_weight(starting_edge) > 0:
                    continue
                
                edge = starting_edge
                add_sequence = False
                while True:  # Until we get back to the starting point.
                    used_edges.add(edge)
                    if add_sequence:  # Only record the edge in the sequence once we have made a right turn away from the vertex.
                        sequence.append(edge)
                    
                    # Move around to the next edge following the lamination.
                    edge = lamination.triangulation.corner_lookup[~edge].edges[2 if lamination.left_weight(~edge) > 0 else 1]
                    
                    add_sequence = add_sequence or lamination.right_weight(edge) <= 0
                    if edge == starting_edge:
                        break
            
            if sequence:
                multiarc = curver.kernel.IntegralLamination(lamination.triangulation, lamination.triangulation.cut_sequence_intersections(sequence))
                # Since multiarc only intersects edges that lamination does, its shorten will never affect an edge that has been frozen.
                _, sub_conjugator = multiarc.shorten()
                conjugator = sub_conjugator * conjugator
                lamination = sub_conjugator(lamination)
                peripheral = sub_conjugator(peripheral)
        
        # Rebuild the image of self under conjugator from its components.
        short = lamination.triangulation.disjoint_sum(dict(
            [(peripheral, 1)]
            + [(lamination.triangulation.edge_arc(edge), multiplicity) for edge, multiplicity in arc_components.items()]
            + [(lamination.triangulation.edge_curve(edge), multiplicity) for edge, multiplicity in curve_components.items()]
            ))
        
        return short, conjugator

