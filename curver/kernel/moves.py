
''' A module for representing basic ways of changing triangulations.
These moves can also track how laminations and homology classes move through those changes. '''

import curver

class Move(object):
    ''' A basic move from one triangulation to another. '''
    def __init__(self, source_triangulation, target_triangulation):
        assert isinstance(source_triangulation, curver.kernel.Triangulation)
        assert isinstance(target_triangulation, curver.kernel.Triangulation)
        
        self.source_triangulation = source_triangulation
        self.target_triangulation = target_triangulation
        self.zeta = self.source_triangulation.zeta
    def __repr__(self):
        return str(self)
    def __invert__(self):
        return self.inverse()
    def __call__(self, other):
        if isinstance(other, curver.kernel.Lamination):
            return self.apply_lamination(other)
        elif isinstance(other, curver.kernel.HomologyClass):
            return self.apply_homology(other)
        else:
            raise TypeError('Unknown type %s' % other)
    
    def encode(self):
        ''' Return the Encoding induced by this move. '''
        
        return curver.kernel.Encoding([self])
    
    def package(self):
        ''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
        
        return self
    
    def inverse(self):  # pylint: disable=no-self-use
        ''' Return the inverse of this move. '''
        
        return NotImplemented
    def apply_lamination(self, lamination):  # pylint: disable=no-self-use,unused-argument
        ''' Return the lamination obtained by mapping the given lamination through this move. '''
        
        return NotImplemented
    def apply_homology(self, homology_class):  # pylint: disable=no-self-use,unused-argument
        ''' Return the homology class obtained by mapping the given homology class through this move. '''
        
        return NotImplemented

class FlipGraphMove(Move):
    ''' A Move between two triangulations in the same flip graph. '''
    def encode(self):
        if self.source_triangulation != self.target_triangulation:
            return curver.kernel.Mapping([self])
        else:
            return curver.kernel.MappingClass([self])
    def flip_mapping(self):  # pylint: disable=no-self-use
        ''' Return a Mapping equal to self.encoding() but that only uses EdgeFlips and Isometries. '''
        
        return NotImplemented

class Isometry(FlipGraphMove):
    ''' This represents an isometry from one Triangulation to another.
    
    Triangulations can create the isometries between themselves and this
    is the standard way users are expected to create these. '''
    def __init__(self, source_triangulation, target_triangulation, label_map):
        ''' This represents an isometry from source_triangulation to target_triangulation.
        
        It is given by a map taking each edge label of source_triangulation to a label of target_triangulation.
        
        Assumes that this map is defined on all labels. '''
        
        super(Isometry, self).__init__(source_triangulation, target_triangulation)
        
        assert isinstance(label_map, dict)
        self.label_map = dict(label_map)
        
        # Quick sanity check.
        assert all(i in self.label_map for i in self.source_triangulation.labels)
        
        self.index_map = dict((i, curver.kernel.norm(self.label_map[i])) for i in self.source_triangulation.indices)
        # Store the inverses too while we're at it.
        self.inverse_label_map = dict((self.label_map[label], label) for label in self.source_triangulation.labels)
        self.inverse_index_map = dict((index, curver.kernel.norm(self.inverse_label_map[index])) for index in self.source_triangulation.indices)
    
    def __str__(self):
        return 'Isometry ' + str([curver.kernel.Edge(self.label_map[index]) for index in self.source_triangulation.indices])
    def __reduce__(self):
        return (self.__class__, (self.source_triangulation, self.target_triangulation, self.label_map))
    def package(self):
        if not all(self.label_map[i] == i for i in self.source_triangulation.indices):  # If self is not the identity isometry.
            return {i: self.label_map[i] for i in self.source_triangulation.labels}
        else:
            return None
    def __eq__(self, other):
        if isinstance(other, Isometry):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation and self.label_map == other.label_map
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    
    def apply_lamination(self, lamination):
        geometric = [lamination(self.inverse_index_map[index]) for index in self.source_triangulation.indices]
        return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
    
    def apply_homology(self, homology_class):
        algebraic = [homology_class(self.inverse_label_map[index]) for index in self.source_triangulation.indices]
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def inverse(self):
        
        return Isometry(self.target_triangulation, self.source_triangulation, self.inverse_label_map)
    
    def flip_mapping(self):
        return self.encode()

class EdgeFlip(FlipGraphMove):
    ''' Represents the change to a curve caused by flipping an edge. '''
    def __init__(self, source_triangulation, target_triangulation, edge):
        super(EdgeFlip, self).__init__(source_triangulation, target_triangulation)
        
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        self.edge = edge
        assert self.source_triangulation.is_flippable(self.edge)
        
        self.square = self.source_triangulation.square(self.edge)
    
    def __str__(self):
        return 'Flip %s' % self.edge
    def __reduce__(self):
        return (self.__class__, (self.source_triangulation, self.target_triangulation, self.edge))
    def package(self):
        return self.edge.label
    def __eq__(self, other):
        if isinstance(other, EdgeFlip):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation and self.edge == other.edge
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    
    def apply_lamination(self, lamination):
        L = lamination  # Shorter name.
        a, b, c, d, e = self.square
        ai, bi, ci, di, ei = [max(L(edge), 0) for edge in self.square]
        
        # Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
        geometric = list(L.geometric)
        
        if L(e) >= ai + bi and ai >= di and bi >= ci:  # CASE: A(ab)
            geometric[e.index] = ai + bi - L(e)
        elif L(e) >= ci + di and di >= ai and ci >= bi:  # CASE: A(cd)
            geometric[e.index] = ci + di - L(e)
        elif L(e) <= 0 and ai >= bi and di >= ci:  # CASE: D(ad)
            geometric[e.index] = ai + di - L(e)
        elif L(e) <= 0 and bi >= ai and ci >= di:  # CASE: D(bc)
            geometric[e.index] = bi + ci - L(e)
        elif L(e) >= 0 and ai >= bi + L(e) and di >= ci + L(e):  # CASE: N(ad)
            geometric[e.index] = ai + di - 2*L(e)
        elif L(e) >= 0 and bi >= ai + L(e) and ci >= di + L(e):  # CASE: N(bc)
            geometric[e.index] = bi + ci - 2*L(e)
        elif ai + bi >= L(e) and bi + L(e) >= 2*ci + ai and ai + L(e) >= 2*di + bi:  # CASE: N(ab)
            geometric[e.index] = (ai + bi - L(e)) // 2
        elif ci + di >= L(e) and di + L(e) >= 2*ai + ci and ci + L(e) >= 2*bi + di:  # CASE: N(cd)
            geometric[e.index] = (ci + di - L(e)) // 2
        else:
            geometric[e.index] = max(ai + ci, bi + di) - L(e)
        
        return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
    
    def apply_homology(self, homology_class):
        a, b, c, d, e = self.square
        
        algebraic = list(homology_class)
        # Move the homology on e onto a & b.
        algebraic[a.index] -= a.sign() * homology_class(e)
        algebraic[b.index] -= b.sign() * homology_class(e)
        algebraic[e.index] = 0
        
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def inverse(self):
        return EdgeFlip(self.target_triangulation, self.source_triangulation, ~self.edge)
    
    def flip_mapping(self):
        return self.encode()

