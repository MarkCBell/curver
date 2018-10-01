
''' A module for representing and manipulating maps between Triangulations. '''

from collections import defaultdict, namedtuple
from fractions import Fraction
from itertools import groupby
import operator
import numpy as np

import curver
from curver.kernel.decorators import ensure, memoize

NT_TYPE_PERIODIC = 'Periodic'
NT_TYPE_REDUCIBLE = 'Reducible'  # Strictly this  means 'reducible and not periodic'.
NT_TYPE_PSEUDO_ANOSOV = 'Pseudo-Anosov'

class Encoding(object):
    ''' This represents a map between two Triangulations.
    
    The map is given by a sequence of Moves which act from right to left. '''
    def __init__(self, sequence):
        assert isinstance(sequence, (list, tuple))
        assert sequence
        # assert all(isinstance(item, curver.kernel.Move) for item in sequence)  # Quadratic.
        
        self.sequence = sequence
        
        self.source_triangulation = self.sequence[-1].source_triangulation
        self.target_triangulation = self.sequence[0].target_triangulation
        self.zeta = self.source_triangulation.zeta
    
    def __repr__(self):
        return str(self)
    def __str__(self):
        return 'Encoding %s' % self.sequence
    def __iter__(self):
        return iter(self.sequence)
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self, value):
        if isinstance(value, slice):
            # It turns out that handling all slices correctly is really hard.
            # We need to be very careful with "empty" slices. As Encodings require
            # non-empty sequences, we have to return just the id_encoding. This
            # ensures the Encoding that we return satisfies:
            #   self == self[:i] * self[i:j] * self[j:]
            # even when i == j.
            
            start = 0 if value.start is None else value.start if value.start >= 0 else len(self) + value.start
            stop = len(self) if value.stop is None else value.stop if value.stop >= 0 else len(self) + value.stop
            if start == stop:
                if 0 <= start < len(self):
                    return self.sequence[start].target_triangulation.id_encoding()
                elif start == len(self):
                    return self.source_triangulation.id_encoding()
                else:
                    raise IndexError('list index out of range')
            elif stop < start:
                raise IndexError('list index out of range')
            else:  # start < stop.
                triangulation = self.sequence[stop-1].source_triangulation
                return triangulation.encode(self.sequence[value])
        elif isinstance(value, curver.IntegerType):
            return self.sequence[value]
        else:
            return NotImplemented
    def package(self):
        ''' Return a small amount of info that self.source_triangulation can use to reconstruct this triangulation. '''
        return [item.package() for item in self]
    def __reduce__(self):
        return (create_encoding, (self.source_triangulation, self.package()))
    
    def __eq__(self, other):
        if isinstance(other, Encoding):
            if self.source_triangulation != other.source_triangulation or self.target_triangulation != other.target_triangulation:
                return False
            
            return all(self(arc.boundary()) == other(arc.boundary()) for arc in self.source_triangulation.edge_arcs())
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        # In fact this hash is perfect unless the surface is S_{1,1}.
        return hash(tuple(entry for arc in self.source_triangulation.edge_arcs() for entry in self(arc.boundary())))
    
    def __call__(self, other):
        if self.source_triangulation != other.triangulation:
            raise ValueError('Cannot apply an Encoding to something on a triangulation other than source_triangulation.')
        
        for item in reversed(self.sequence):
            other = item(other)
        
        return other
    def __mul__(self, other):
        if isinstance(other, Encoding):
            if self.source_triangulation != other.target_triangulation:
                raise ValueError('Cannot compose Encodings over different triangulations.')
            
            if not (isinstance(self, Mapping) and isinstance(other, Mapping)):
                return Encoding(self.sequence + other.sequence)
            else:  # self and other both at least Mappings:
                if self.target_triangulation != other.source_triangulation:
                    return Mapping(self.sequence + other.sequence)
                else:  # self.target_triangulation == other.source_triangulation:
                    return MappingClass(self.sequence + other.sequence)
        elif other is None:
            return self
        else:
            return NotImplemented
    def inverse(self):
        ''' Return the inverse of this encoding. '''
        
        return self.__class__([item.inverse() for item in reversed(self.sequence)])  # Data structure issue.
    def __invert__(self):
        return self.inverse()

class Mapping(Encoding):
    ''' An Encoding where every move is a FlipGraphMove.
    
    Hence this encoding is a sequence of moves in the same flip graph. '''
    def __str__(self):
        return 'Mapping %s' % self.sequence
    def intersection_matrix(self):
        ''' Return the matrix M = {signed_intersection(self(e_i), e'_j)}_{ij}.
        Here e_i and e'_j are the edges of self.source_triangulation and self.target_triangulation respectively.
        
        Except when on S_{1,1}, this uniquely determines self. '''
        
        return np.array([list(self(arc)) for arc in self.source_triangulation.edge_arcs()], dtype=object)
    
    def homology_matrix(self):
        ''' Return a matrix describing the action of this mapping on first homology (relative to the punctures).
        
        The matrix is given with respect to the homology bases of the source and target triangulations. '''
        
        source_basis = self.source_triangulation.homology_basis()
        target_basis = self.target_triangulation.homology_basis()
        
        source_images = [self(hc).canonical() for hc in source_basis]
        
        return np.array([[sum(x * y for x, y in zip(hc, hc2)) for hc in source_images] for hc2 in target_basis], dtype=object)
    def __eq__(self, other):
        if isinstance(other, Encoding):
            if self.source_triangulation != other.source_triangulation or self.target_triangulation != other.target_triangulation:
                return False
            
            tri_lamination = self.source_triangulation.as_lamination()
            return self(tri_lamination) == other(tri_lamination) and \
                all(self(hc) == other(hc) for hc in self.source_triangulation.edge_homologies())  # We only really need this for S_{1,1}.
        else:
            return NotImplemented
    def __hash__(self):
        return hash(tuple(entry for row in self.intersection_matrix().tolist() for entry in row))
    def vertex_map(self):
        ''' Return the dictionary (vertex, self(vertex)) for each vertex in self.source_triangulation.
        
        When self is a MappingClass this is a permutation of the vertices. '''
        
        source_vertices = dict((vertex, self.source_triangulation.curve_from_cut_sequence(vertex)) for vertex in self.source_triangulation.vertices)
        target_vertices_inverse = dict((self.target_triangulation.curve_from_cut_sequence(vertex), vertex) for vertex in self.target_triangulation.vertices)
        
        return dict((vertex, target_vertices_inverse[self(source_vertices[vertex])]) for vertex in self.source_triangulation.vertices)
    
    def simplify(self):
        ''' Return a new Mapping that is equal to self.
        
        This is obtained by combing the image of the source triangulation under self and so is (hopefully) simpler than self since it depends only on the endpoints. '''
        
        conjugator = self(self.source_triangulation.as_lamination()).shorten()
        # conjugator.inverse() is almost self, however the edge labels might not agree.
        potential_closers = [isom.encode() for isom in conjugator.target_triangulation.isometries_to(self.source_triangulation)]
        identity = self.source_triangulation.id_encoding()
        [closer] = [potential_closer for potential_closer in potential_closers if potential_closer * conjugator * self == identity]  # There should only be one.
        
        return (closer * conjugator).inverse()
    
    def flip_mapping(self):
        ''' Return a Mapping equal to self that only uses EdgeFlips and Isometries. '''
        
        return self.__class__([item for move in self for item in move.flip_mapping()])

class MappingClass(Mapping):
    ''' A Mapping from a Triangulation to itself.
    
    That is, one where self.source_triangulation == self.target_triangulation. '''
    def __str__(self):
        return 'MappingClass %s' % self.sequence
    def __pow__(self, k):
        if k == 0:
            return self.source_triangulation.id_encoding()
        elif k > 0:
            return MappingClass(self.sequence * k)
        else:
            return self.inverse()**abs(k)
    
    def is_in_torelli(self):
        ''' Return whether this mapping class is in the Torelli subgroup. '''
        
        homology_matrix = self.homology_matrix()
        return np.array_equal(homology_matrix, np.identity(homology_matrix.shape[0]))
    
    @memoize
    def order(self):
        ''' Return the order of this mapping class.
        
        If this has infinite order then return 0. '''
        
        # There are several tricks that we could use to rule out powers that we need to test:
        #  - If self**i == identity then homology_matrix**i = identity_matrix.
        #  - If self**i == identity then self**(ij) == identity.
        #  - if self**i == identity and self**j == identity then self**gcd(i, j) == identity.
        #
        # But in terms of raw speed there doesn't appear to be anything faster than:
        
        homology_matrix = self.homology_matrix()
        originals = [np.identity(homology_matrix.shape[0]), self.source_triangulation.as_lamination()]
        images = list(originals)
        cmps = [np.array_equal, operator.eq]
        applies = [homology_matrix.dot, self]
        powers = [0, 0]
        for power in range(1, self.source_triangulation.max_order()+1):
            for i in range(len(originals)):  # pylint: disable=consider-using-enumerate
                while powers[i] < power:
                    images[i] = applies[i](images[i])
                    powers[i] += 1
                if not cmps[i](images[i], originals[i]):
                    break
            else:
                return power
        
        return 0
    
    def vertex_permutation(self):
        ''' Return a permutation describing how the vertices of self.source_triangulation (labelled in sorted order) are permuted. '''
        
        vertex_map = self.vertex_map()
        return curver.kernel.Permutation.from_dict(vertex_map, ordering=sorted(vertex_map))
    
    def is_identity(self):
        ''' Return if this mapping class is the identity. '''
        
        return self.order() == 1
    
    def is_periodic(self):
        ''' Return whether this mapping class has finite order. '''
        
        return self.order() > 0
    
    def is_reducible(self):
        ''' Return whether this mapping class is reducible. '''
        
        if self.is_periodic():
            # A periodic mapping class is reducible iff at least one of the components of its quotient orbifold is not a triangle orbifold.
            # The genus of the surface underlying an orbifold.
            genus = lambda orbifold: (2 - orbifold.euler_characteristic - sum(1 - (0 if cone_point.punctured else Fraction(1, cone_point.rotation.denominator)) for cone_point in orbifold.cone_points)) // 2
            return not all(len(orbifold.cone_points) == 3 and genus(orbifold) == 0 for orbifold in self.quotient_orbifold_signature())
        else:
            return not self.positive_asymptotic_translation_length()
    
    def is_pseudo_anosov(self):
        ''' Return whether this mapping class is pseudo-Anosov. '''
        
        return not self.is_periodic() and not self.is_reducible()
    
    def nielsen_thurston_type(self):
        ''' Return the Nielsen--Thurston type of this mapping class. '''
        
        if self.is_periodic():
            return NT_TYPE_PERIODIC
        elif self.is_reducible():
            return NT_TYPE_REDUCIBLE
        else:  # self.is_pesudo_anosov():
            return NT_TYPE_PSEUDO_ANOSOV
    
    def asymptotic_translation_length(self):
        ''' Return the asymptotic translation length of this mapping class on the curve complex.
        
        From Algorithm 6 of [BellWebb16]_. '''
        
        C = curver.kernel.CurveGraph(self.source_triangulation)
        c = self.source_triangulation.edge_arc(0).boundary()  # A "short" curve.
        geodesic = C.geodesic(c, (self**C.M)(c))
        m = geodesic[len(geodesic)//2]  # midpoint
        
        numerator = C.distance(m, (self**C.M)(m))
        denominator = C.M
        return Fraction(numerator, denominator).limit_denominator(C.D)
    
    def positive_asymptotic_translation_length(self):
        ''' Return whether the asymptotic translation length of this mapping class on the curve complex is positive.
        
        This uses Remark 4.7 of [BellWebb16]_ which is based on [GadreTsai11]_ and so is more efficient than doing::
        
            self.asymptotic_translation_length() > 0 '''
        
        C = curver.kernel.CurveGraph(self.source_triangulation)
        c = self.source_triangulation.edge_arc(0).boundary()  # A "short" curve.
        geodesic = C.geodesic(c, (self**C.M2)(c))
        m = geodesic[len(geodesic)//2]  # midpoint
        
        return C.distance(m, (self**C.M2)(m)) > 4
    
    @memoize
    @ensure(lambda data: data.result.is_polygonalisation(), lambda data: data.self(data.result) == data.result)
    def invariant_polygonalisation(self):
        ''' Return a multiarc that is a polygonalisation and is invariant under self.
        
        self must be a periodic mapping class. '''
        
        assert self.is_periodic()
        
        def orbit(a):
            ''' Yield the orbit of a under h (self conjugated by conjugator). '''
            
            # a is rarely a useful image so we will yield it last.
            image = h(a)
            while image != a:
                yield image
                image = h(image)
            yield a
        
        def unicorns(a, b):
            ''' Yield a collection of arcs that includes all unicorn arcs that can be made with a & b.
            
            Assumes that a is short. The general version of this function would begin by shortening a. '''
            
            assert a.is_short()
            assert b.triangulation == a.triangulation
            
            conjugator = b.shorten(drop=0)
            conjugator_inv = conjugator.inverse()
            for i in range(1, len(conjugator)-1):
                prefix_inv = conjugator_inv[:i]
                if isinstance(prefix_inv[-1], curver.kernel.EdgeFlip):
                    arc = prefix_inv.source_triangulation.edge_arc(prefix_inv[-1].edge)
                    yield prefix_inv(arc)
            for arc in b.triangulation.edge_arcs():
                yield arc
        
        def orbit_unicorns(arc):
            ''' Yield a collection of arcs including all unicorn arcs that can be made from arc and h^i(arc) for each i. '''
            
            for image in orbit(arc):
                for unicorn in unicorns(arc, image):
                    yield unicorn
        
        h = self
        conjugator = self.source_triangulation.id_encoding()
        triangulation = h.source_triangulation
        invariant_multiarc = self.source_triangulation.empty_lamination()
        while not invariant_multiarc.is_polygonalisation():  # Loops at most zeta times.
            # Find a new multiarc that is not a component of invariant_multiarc, is disjoint from invariant_multiarc and is invariant under h.
            # Start by finding an arc that does not cut off a disk in S - invariant_multiarc.
            # Since invariant_multiarc is not a polygonalisation, one exists and in fact one of the edges of triangulation must be one since it is short.
            dual_tree = triangulation.dual_tree(avoid={edge for edge in triangulation.positive_edges if invariant_multiarc(edge) < 0})
            arc = triangulation.edge_arc([edge for edge in triangulation.positive_edges if edge.index not in dual_tree and invariant_multiarc(edge) == 0][0])
            
            for unicorn in orbit_unicorns(arc):  # Loops at most zeta^2 * ||self|| times.
                # Perform tests in order of difficulty.
                if unicorn in invariant_multiarc.components():
                    continue
                if invariant_multiarc.intersection(unicorn) != 0:
                    continue
                unicorn_orbit = list(orbit(unicorn))  # Save result for performance.
                if unicorn.intersection(*unicorn_orbit) != 0:
                    continue
                invariant_multiarc = triangulation.disjoint_sum([invariant_multiarc] + unicorn_orbit)
                break
            
            # Reshorten invariant_multiarc.
            next_conjugator = invariant_multiarc.shorten()
            conjugator = next_conjugator * conjugator
            invariant_multiarc = next_conjugator(invariant_multiarc)
            h = next_conjugator * h * next_conjugator.inverse()
            triangulation = h.source_triangulation
        
        return conjugator.inverse()(invariant_multiarc)
    
    @memoize
    def quotient_orbifold_signature(self):
        ''' Return the signature of self.surface() / self.
        
        Since this also records the covering map via the `preimage` and
        `rotation` fields, this is a total conjugacy invariant for periodic
        mapping classes by Theorem 9 of [Mosher07]_.
        
        Assumes that self is periodic. '''
        
        assert self.is_periodic()
        
        polygonalisation = self.invariant_polygonalisation()
        
        conjugator = polygonalisation.shorten()
        short = conjugator(polygonalisation)
        
        h = conjugator * self * conjugator.inverse()
        
        # Some short names.
        h_order = h.order()
        triangulation = short.triangulation
        components = short.triangulation.components()
        surface = triangulation.surface()
        
        OrientedArc = namedtuple('OrientedArc', ['arc', 'hc', 'boundary'])
        # Build the oriented arcs.
        oriented = dict()
        for edge in triangulation.edges:
            if short(edge) < 0:
                arc = triangulation.edge_arc(edge)
                if len(arc.vertices()) == 1:
                    [v] = arc.vertices()
                    v_edges = curver.kernel.utilities.cyclic_slice(v, edge, ~edge)
                    boundary = triangulation.curve_from_cut_sequence(v_edges[1:])
                else:  # two vertices:
                    boundary = arc.boundary()
                oriented[edge] = OrientedArc(arc, triangulation.edge_homology(edge), boundary)
        oriented_arcs = list(oriented.values())
        
        # Some useful maps.
        # How h permutes the oriented arcs.
        h_oriented = dict((oriented_arc, OrientedArc(h(oriented_arc.arc), h(oriented_arc.hc), h(oriented_arc.boundary))) for oriented_arc in oriented_arcs)
        # The component of triangulation each oriented_arc lives in.
        component_lookup = dict((oriented[edge], component) for component in components for edge in component if edge in oriented)
        # The length of the orbit of each component of triangulation under the action of h.
        component_orbit_length = dict()
        for oriented_arc in oriented_arcs:
            component = component_lookup[oriented_arc]
            if component not in component_orbit_length:  # Save repeating calculations.
                image = oriented_arc
                for i in range(1, h_order+1):
                    image = h_oriented[image]
                    if component_lookup[image] == component:
                        component_orbit_length[component] = i
                        break
        # The (orbifold) Euler characteristic of each components quotient.
        euler_characteristic = dict((component, Fraction((2 - 2*surface[component][0] - surface[component][1]) * component_orbit_length[component], h_order)) for component in components)
        
        # Compute the polygons cut out by short.
        # Remember to walk around their boundary in the correct direction.
        polygons = []
        used = set()
        for edge in triangulation.edges:
            if short(edge) < 0 and edge not in used:
                polygon_edges = [edge]
                while True:
                    # Correct direction requires taking the last oriented arc out of the vertex each time.
                    polygon_edges.append([edgy for edgy in curver.kernel.utilities.cyclic_slice(triangulation.vertex_lookup[~polygon_edges[-1]], ~polygon_edges[-1]) if short(edgy) < 0][-1])
                    if polygon_edges[-1] == polygon_edges[0]:  # if back where we started.
                        polygon_edges = polygon_edges[:-1]  # Remeber to discard the last edge as it duplicates the first.
                        break
                
                used = used.union(polygon_edges)  # Mark everything as used.
                polygons.append(polygon_edges)
        
        ConePoint = namedtuple('ConePoint', ['punctured', 'rotation', 'preimages'])
        # There are three places to look for cone points:
        # 1) at vertices,
        # 2) at the midpoints of edges, and
        # 3) at the centres of polygons (note that these are not once-punctured polygons since we started with an invariant polygonalisation).
        candidates = [(True, [oriented[edge] for edge in vertex if short(edge) < 0]) for vertex in triangulation.vertices] + \
            [(False, [oriented[edge], oriented[~edge]]) for edge in triangulation.positive_edges if short(edge) < 0] + \
            [(False, [oriented[edge] for edge in polygon]) for polygon in polygons]
        
        cone_points = defaultdict(list)
        for punctured, oriented_arcs in candidates:
            oriented_arc = image = oriented_arcs[0]
            for i in range(1, h_order+1):
                image = h_oriented[image]
                if image in oriented_arcs:
                    rotation = Fraction(oriented_arcs.index(image), len(oriented_arcs))
                    if punctured or rotation > 0:  # Real cone point.
                        cone_points[component_lookup[oriented_arc]].append(ConePoint(punctured, rotation, i))  # Put it in the correct component.
                    break
        
        # Remember to make all the data canonical by sorting.
        Orbifold = namedtuple('Orbifold', ['euler_characteristic', 'preimages', 'cone_points'])
        signature = sorted(Orbifold(euler_characteristic[component], component_orbit_length[component], sorted(cone_points[component])) for component in components)
        
        # Group.
        signature = [key for key, group in groupby(signature) for _ in range(len(list(group)) // key.preimages)]
        signature = [Orbifold(orbifold.euler_characteristic, orbifold.preimages, [key for key, group in groupby(orbifold.cone_points) for _ in range(len(list(group)) // key.preimages)]) for orbifold in signature]
        return signature

    def is_conjugate_to(self, other):
        ''' Return whether this mapping class is conjugate to other.
        
        It would also be straightforward to check whether self^i ~~ other^j for some i, j.
        
        Currently assumes that at least one mapping class is periodic. '''
        
        assert isinstance(other, curver.kernel.MappingClass)
        
        if self.source_triangulation.surface() != other.source_triangulation.surface():  # Defined on different surfaces.
            return False
        
        if not self.vertex_permutation().is_conjugate_to(other.vertex_permutation()):  # Induce non-conjugate permutations of the vertices.
            return False
        
        if self.is_periodic() != other.is_periodic():
            return False
        
        if self.is_periodic():
            if self.order() != other.order():  # Conjugacy invariant.
                return False
            return self.quotient_orbifold_signature() == other.quotient_orbifold_signature()  # Total conjugacy invariant.
        else:
            raise ValueError('is_conjugate_to is currently only implemented when one of the mapping classes is periodic. Consider using flipper.')

def create_encoding(source_triangulation, sequence):
    ''' Return the encoding defined by sequence starting at source_triangulation.
    
    This is only really here to help with pickling. Users should use
    source_triangulation.encode(sequence) directly. '''
    
    assert isinstance(source_triangulation, curver.kernel.Triangulation)
    
    return source_triangulation.encode(sequence)

