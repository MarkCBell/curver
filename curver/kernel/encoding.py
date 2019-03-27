
''' A module for representing and manipulating maps between Triangulations. '''

from fractions import Fraction
import operator
import numpy as np

import curver
from curver.kernel.decorators import memoize

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
        
        is_lamination = isinstance(other, curver.kernel.Lamination)
        is_homology = isinstance(other, curver.kernel.HomologyClass)
        if not is_lamination and not is_homology: raise TypeError('Unknown type %s' % other)
        
        for item in reversed(self):
            if is_lamination:
                other = item.apply_lamination(other)
            elif is_homology:
                other = item.apply_homology(other)
        
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
    
    @memoize
    def self_image(self):
        ''' Return the image of self.source_triangulation under self. '''
        return self(self.source_triangulation.as_lamination())
    
    @memoize
    def intersection_matrix(self):
        ''' Return the matrix M = {signed_intersection(self(e_i), e'_j)}_{ij}.
        Here e_i and e'_j are the edges of self.source_triangulation and self.target_triangulation respectively.
        
        Except when on S_{1,1}, this uniquely determines self. '''
        
        return np.array([list(self(arc)) for arc in self.source_triangulation.edge_arcs()], dtype=object)
    
    @memoize
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
            
            return self.self_image() == other.self_image() and np.array_equal(self.homology_matrix(), other.homology_matrix())  # We only really need this for S_{1,1}.
        else:
            return NotImplemented
    
    @memoize
    def __hash__(self):
        return hash(self.self_image())
    
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
        
        potential_closers = [isom.encode() for isom in self.source_triangulation.isometries_to(conjugator.target_triangulation)]
        
        # We used to test:
        #   if potential_closer.inverse() * conjugator * self == identity
        # However since conjugator * self preserves self.source_triangulation.as_lamination() it is periodic, and so we can find the correct closer simply
        # by finding the one that induces the same action on homology since the only periodic mapping class in the Torelli group is the identity.
        
        homology_images = [conjugator(self(hc)) for hc in self.source_triangulation.edge_homologies()]
        [closer] = [potential_closer for potential_closer in potential_closers if all(potential_closer(hc) == hci for hc, hci in zip(self.source_triangulation.edge_homologies(), homology_images))]
        
        return conjugator.inverse() * closer
    
    def flip_mapping(self):
        ''' Return a Mapping equal to self that only uses EdgeFlips and Isometries. '''
        
        return self.__class__([item for move in self for item in move.flip_mapping()])

class MappingClass(Mapping):
    ''' A Mapping from a Triangulation to itself.
    
    That is, one where self.source_triangulation == self.target_triangulation. '''
    def __call__(self, other, power=1):
        h = self if power >= 0 else self.inverse()
        power = abs(power)
        for _ in range(power):
            other = super(MappingClass, h).__call__(other)
        return other
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
        return np.array_equal(homology_matrix, np.identity(homology_matrix.shape[0], dtype=object))
    
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
        originals = [np.identity(homology_matrix.shape[0], dtype=object), self.source_triangulation.as_lamination()]
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
            genus = lambda orbifold: (2 - orbifold.euler_characteristic - sum(1 - (0 if cone_point.punctured else Fraction(1, cone_point.order)) for cone_point in orbifold.cone_points)) // 2
            return not all(len(orbifold.cone_points) == 3 and genus(orbifold) == 0 for orbifold in self.subgroup().quotient_orbifold_signature())
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
        geodesic = C.geodesic(c, self(c, power=C.M))
        m = geodesic[len(geodesic)//2]  # midpoint
        
        numerator = C.distance(m, self(m, power=C.M))
        denominator = C.M
        return Fraction(numerator, denominator).limit_denominator(C.D)
    
    def positive_asymptotic_translation_length(self):
        ''' Return whether the asymptotic translation length of this mapping class on the curve complex is positive.
        
        This uses Remark 4.7 of [BellWebb16]_ which is based on [GadreTsai11]_ and so is more efficient than doing::
        
            self.asymptotic_translation_length() > 0 '''
        
        C = curver.kernel.CurveGraph(self.source_triangulation)
        c = self.source_triangulation.edge_arc(0).boundary()  # A "short" curve.
        geodesic = C.geodesic(c, self(c, power=C.M2))
        m = geodesic[len(geodesic)//2]  # midpoint
        
        return C.distance(m, self(m, power=C.M2)) > 4
    
    @memoize
    def subgroup(self):
        ''' Return the FiniteSubgroup generated by this mapping class.
        
        This mapping class must be periodic. '''
        
        order = self.order()
        if order == 0:  # self is not periodic.
            raise ValueError('MappingClass is not periodic.')
        
        return curver.kernel.FiniteSubgroup({i: self**i for i in range(order)}, [0 if order == 1 else 1])

    def is_conjugate_to(self, other):
        ''' Return whether this mapping class is conjugate to other.
        
        It would also be straightforward to check whether self^i ~~ other^j for some i, j.
        
        In the periodic case we use the quotient orbifold and its covering map, the covering map is recorded via the `preimage` and `holonomy` fields.
        This is a total conjugacy invariant for periodic mapping classes by Theorem 9 of [Mosher07]_.
        
        Currently, at least one mapping class must be is periodic. '''
        
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
            return self.subgroup().quotient_orbifold_signature() == other.subgroup().quotient_orbifold_signature()  # Total conjugacy invariant.
        else:
            raise ValueError('is_conjugate_to is currently only implemented when one of the mapping classes is periodic. Consider using flipper.')
    
    def extract_twisting_multicurve(self):
        ''' Return a MultiCurve c such that c.encode_twist() == self.
        
        This raises a ValueError if no such MultiCurve exists.'''
        
        lamination = self.source_triangulation.as_lamination()
        for _ in range(self.zeta):
            image = self(lamination)
            try:
                multicurve = self.source_triangulation([x - y for x, y in zip(image, lamination)])
                if isinstance(multicurve, curver.kernel.MultiCurve):
                    weighted_multicurve = self.source_triangulation.disjoint_sum([(multiplicity // lamination.intersection(component)) * component for component, multiplicity in multicurve.components().items()])
                    if not weighted_multicurve.is_empty() and self == weighted_multicurve.encode_twist():
                        return weighted_multicurve
            except AssertionError:
                pass
            lamination = image
        
        raise ValueError('Mapping Class is not a twist about a MultiCurve.')

def create_encoding(source_triangulation, sequence):
    ''' Return the encoding defined by sequence starting at source_triangulation.
    
    This is only really here to help with pickling. Users should use
    source_triangulation.encode(sequence) directly. '''
    
    assert isinstance(source_triangulation, curver.kernel.Triangulation)
    
    return source_triangulation.encode(sequence)

