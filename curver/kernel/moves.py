
''' A module for representing basic ways of changing triangulations.
These moves can also track how laminations and homology classes move through those changes. '''

from abc import ABC, abstractmethod
import numpy as np

import curver

class Move(ABC):
    ''' A basic move from one triangulation to another. '''
    def __init__(self, source_triangulation, target_triangulation):
        assert isinstance(source_triangulation, curver.kernel.Triangulation)
        assert isinstance(target_triangulation, curver.kernel.Triangulation)
        
        self.source_triangulation = source_triangulation
        self.target_triangulation = target_triangulation
        self.zeta = self.source_triangulation.zeta
        self._inverse = None
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
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation
        else:
            return NotImplemented
    
    def encode(self):
        ''' Return the Encoding induced by this move. '''
        
        return curver.kernel.Encoding([self]).promote()
    
    def package(self):
        ''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
    
    def inverse(self):
        ''' Return the inverse of this move. '''
        
        return self._inverse
    
    @abstractmethod
    def package(self):
        ''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
    
    @abstractmethod
    def apply_lamination(self, lamination):
        ''' Return the lamination obtained by mapping the given lamination through this move. '''
    
    @abstractmethod
    def apply_homology(self, homology_class):
        ''' Return the homology class obtained by mapping the given homology class through this move. '''

class FlipGraphMove(Move):
    ''' A Move between two triangulations in the same flip graph. '''
    @abstractmethod
    def flip_mapping(self):
        ''' Return a Mapping equal to self.encoding() but that only uses EdgeFlips and Isometries. '''
        
        return NotImplemented
    
    @abstractmethod
    def pl_action(self, multicurve):
        ''' Return the PartialLinearFunction that this FlipGraphMove applies to the given multicurve. '''
        
        return NotImplemented

class Isometry(FlipGraphMove):
    ''' This represents an isometry from one Triangulation to another.
    
    Triangulations can create the isometries between themselves and this
    is the standard way users are expected to create these. '''
    def __init__(self, source_triangulation, target_triangulation, label_map):
        ''' This represents an isometry from source_triangulation to target_triangulation.
        
        It is given by a map taking each edge label of source_triangulation to a label of target_triangulation.
        
        This map must be defined on all labels. '''
        
        super().__init__(source_triangulation, target_triangulation)
        
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
    def package(self):
        return None if self.is_identity() else self.label_map
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.label_map == other.label_map
    
    def apply_lamination(self, lamination):
        geometric = [lamination(self.inverse_index_map[index]) for index in self.target_triangulation.indices]
        return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
    
    def apply_homology(self, homology_class):
        algebraic = [homology_class(self.inverse_label_map[index]) for index in self.target_triangulation.indices]
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def flip_mapping(self):
        return self.encode()
    
    def is_identity(self):
        ''' Return whether this isometry is the identity. '''
        
        return all(key == value for key, value in self.label_map.items())
    
    def pl_action(self, multicurve):
        action = np.array([[1 if j == self.index_map[i] else 0 for i in range(self.zeta)] for j in range(self.zeta)], dtype=object)
        condition = np.array([[0] * self.zeta], dtype=object)
        return curver.kernel.PartialLinearFunction(action, condition)

class EdgeFlip(FlipGraphMove):
    ''' Represents the change to a curve caused by flipping an edge. '''
    def __init__(self, source_triangulation, target_triangulation, edge):
        super().__init__(source_triangulation, target_triangulation)
        
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        self.edge = edge
        assert self.source_triangulation.is_flippable(self.edge)
        
        self.square = self.source_triangulation.square(self.edge)
    
    def __str__(self):
        return 'Flip %s' % self.edge
    def package(self):
        return self.edge.label
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.edge == other.edge
    
    def apply_lamination(self, lamination):
        ''' See Lemma 5.1.3 of [Bell15]_ for details of the cases involved in performing a flip. '''
        
        ei = lamination(self.edge)
        ai0, bi0, ci0, di0, ei0 = [max(lamination(edge), 0) for edge in self.square]
        
        # Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
        geometric = list(lamination.geometric)
        
        if ei >= ai0 + bi0 and ai0 >= di0 and bi0 >= ci0:  # CASE: A(ab)
            geometric[self.edge.index] = ai0 + bi0 - ei
        elif ei >= ci0 + di0 and di0 >= ai0 and ci0 >= bi0:  # CASE: A(cd)
            geometric[self.edge.index] = ci0 + di0 - ei
        elif ei <= 0 and ai0 >= bi0 and di0 >= ci0:  # CASE: D(ad)
            geometric[self.edge.index] = ai0 + di0 - ei
        elif ei <= 0 and bi0 >= ai0 and ci0 >= di0:  # CASE: D(bc)
            geometric[self.edge.index] = bi0 + ci0 - ei
        elif ei >= 0 and ai0 >= bi0 + ei and di0 >= ci0 + ei:  # CASE: N(ad)
            geometric[self.edge.index] = ai0 + di0 - 2*ei
        elif ei >= 0 and bi0 >= ai0 + ei and ci0 >= di0 + ei:  # CASE: N(bc)
            geometric[self.edge.index] = bi0 + ci0 - 2*ei
        elif ai0 + bi0 >= ei and bi0 + ei >= 2*ci0 + ai0 and ai0 + ei >= 2*di0 + bi0:  # CASE: N(ab)
            geometric[self.edge.index] = curver.kernel.utilities.half(ai0 + bi0 - ei)
        elif ci0 + di0 >= ei and di0 + ei >= 2*ai0 + ci0 and ci0 + ei >= 2*bi0 + di0:  # CASE: N(cd)
            geometric[self.edge.index] = curver.kernel.utilities.half(ci0 + di0 - ei)
        else:
            geometric[self.edge.index] = max(ai0 + ci0, bi0 + di0) - ei
        
        return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
    
    def apply_homology(self, homology_class):
        a, b, c, d, e = self.square
        
        algebraic = list(homology_class)
        # Move the homology on e onto a & b.
        algebraic[a.index] -= a.sign() * homology_class(e)
        algebraic[b.index] -= b.sign() * homology_class(e)
        algebraic[e.index] = 0
        
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def flip_mapping(self):
        return self.encode()
    
    def pl_action(self, multicurve):
        identity = np.identity(self.zeta, dtype=object)
        
        def V(x):
            ''' Return the vector of length self.zeta which has a 1 at x. '''
            return np.array([1 if i == x else 0 for i in range(self.zeta)], dtype=object)
        
        def E(x, y):
            ''' Return the self.zeta x self.zeta matrix that has a 1 at (x, y). '''
            return np.array([V(x) if j == y else [0] * self.zeta for j in range(self.zeta)], dtype=object)
        
        ai, bi, ci, di, ei = [edge.index for edge in self.square]
        ai0, bi0, ci0, di0, ei0 = [max(multicurve(edge), 0) for edge in self.square]
        if ai0 + ci0 - bi0 - di0 >= 0:
            action = identity + E(ai, ei) + E(ci, ei) - 2*E(ei, ei)
            condition = np.array([V(ai) + V(ci) - V(bi) - V(di)])
        else:
            action = identity + E(bi, ei) + E(di, ei) - 2*E(ei, ei)
            condition = np.array([V(bi) + V(di) - V(ai) - V(ci)])
        
        return curver.kernel.PartialLinearFunction(action, condition)

class MultiEdgeFlip(FlipGraphMove):
    ''' Represents the change to a curve caused by flipping multiple edges. '''
    def __init__(self, source_triangulation, target_triangulation, edges):
        super().__init__(source_triangulation, target_triangulation)
        
        self.edges = set(curver.kernel.Edge(edge) if isinstance(edge, curver.IntegerType) else edge for edge in edges)  # If given any integers.
        self.squares = dict((edge, self.source_triangulation.square(edge)) for edge in self.edges)
        
        support = set(self.source_triangulation.triangle_lookup[e] for edge in edges for e in [edge, ~edge])
        assert len(support) == 2 * len(edges)  # Check disjoint support.
        # Disjoint support implies flippable.
    
    def __str__(self):
        return 'Flips %s' % self.edges
    def package(self):
        return set(edge.label for edge in self.edges)
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.edges == other.edges
    
    def apply_lamination(self, lamination):
        ''' See Lemma 5.1.3 of [Bell15]_ for details of the cases involved in performing a flip. '''
        
        # Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
        geometric = list(lamination.geometric)
        
        for edge in self.edges:
            ei = lamination(edge)
            ai0, bi0, ci0, di0, ei0 = [max(lamination(e), 0) for e in self.squares[edge]]
            if ei >= ai0 + bi0 and ai0 >= di0 and bi0 >= ci0:  # CASE: A(ab)
                geometric[edge.index] = ai0 + bi0 - ei
            elif ei >= ci0 + di0 and di0 >= ai0 and ci0 >= bi0:  # CASE: A(cd)
                geometric[edge.index] = ci0 + di0 - ei
            elif ei <= 0 and ai0 >= bi0 and di0 >= ci0:  # CASE: D(ad)
                geometric[edge.index] = ai0 + di0 - ei
            elif ei <= 0 and bi0 >= ai0 and ci0 >= di0:  # CASE: D(bc)
                geometric[edge.index] = bi0 + ci0 - ei
            elif ei >= 0 and ai0 >= bi0 + ei and di0 >= ci0 + ei:  # CASE: N(ad)
                geometric[edge.index] = ai0 + di0 - 2*ei
            elif ei >= 0 and bi0 >= ai0 + ei and ci0 >= di0 + ei:  # CASE: N(bc)
                geometric[edge.index] = bi0 + ci0 - 2*ei
            elif ai0 + bi0 >= ei and bi0 + ei >= 2*ci0 + ai0 and ai0 + ei >= 2*di0 + bi0:  # CASE: N(ab)
                geometric[edge.index] = curver.kernel.utilities.half(ai0 + bi0 - ei)
            elif ci0 + di0 >= ei and di0 + ei >= 2*ai0 + ci0 and ci0 + ei >= 2*bi0 + di0:  # CASE: N(cd)
                geometric[edge.index] = curver.kernel.utilities.half(ci0 + di0 - ei)
            else:
                geometric[edge.index] = max(ai0 + ci0, bi0 + di0) - ei
        
        return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
    
    def apply_homology(self, homology_class):
        algebraic = list(homology_class)
        
        for edge in self.edges:
            a, b, c, d, e = self.squares[edge]
            
            # Move the homology on e onto a & b.
            algebraic[a.index] -= a.sign() * homology_class(e)
            algebraic[b.index] -= b.sign() * homology_class(e)
            algebraic[e.index] = 0
        
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def flip_mapping(self):
        return self.source_triangulation.encode([edge.label for edge in self.edges])
    
    def pl_action(self, multicurve):
        action = np.identity(self.zeta, dtype=object)
        conditions = []
        
        def V(x):
            ''' Return the vector of length self.zeta which has a 1 at x. '''
            return np.array([1 if i == x else 0 for i in range(self.zeta)], dtype=object)
        
        def E(x, y):
            ''' Return the self.zeta x self.zeta matrix that has a 1 at (x, y). '''
            return np.array([V(x) if j == y else [0] * self.zeta for j in range(self.zeta)], dtype=object)
        
        for edge in self.edges:
            ai, bi, ci, di, ei = [e.index for e in self.squares[edge]]
            ai0, bi0, ci0, di0, ei0 = [max(multicurve(e), 0) for e in self.squares[edge]]
            if ai0 + ci0 - bi0 - di0 >= 0:
                action = action + E(ai, ei) + E(ci, ei) - 2*E(ei, ei)
                conditions.append(V(ai) + V(ci) - V(bi) - V(di))
            else:
                action = action + E(bi, ei) + E(di, ei) - 2*E(ei, ei)
                conditions.append(V(bi) + V(di) - V(ai) - V(ci))
        
        return curver.kernel.PartialLinearFunction(action, np.array(conditions))

class PartialIsometry(Move):
    ''' This represents an isometry from one Triangulation to another.
    
    Triangulations can create the isometries between themselves and this
    is the standard way users are expected to create these. '''
    def __init__(self, source_triangulation, target_triangulation, label_map):
        ''' This represents a partial isometry from source_triangulation to target_triangulation.
        
        It is given by a map taking each edge label of source_triangulation to a label of target_triangulation.
        
        This map must be defined on all labels. '''
        
        super(PartialIsometry, self).__init__(source_triangulation, target_triangulation)
        
        assert isinstance(label_map, dict)
        self.label_map = dict(label_map)
        
        self.index_map = dict((curver.kernel.norm(key), curver.kernel.norm(value)) for key, value in self.label_map.items())
        self.inverse_label_map = dict((value, key) for key, value in self.label_map.items())
        self.inverse_index_map = dict((curver.kernel.norm(value), curver.kernel.norm(key)) for key, value in self.label_map.items())
    
    def __str__(self):
        return 'PartialIsometry ' + str([curver.kernel.Edge(self.label_map[index]) if index in self.label_map else '-' for index in self.source_triangulation.indices])
    def package(self):
        return sorted(self.label_map)
    def __eq__(self, other):
        eq = super(PartialIsometry, self).__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.label_map == other.label_map
    
    def apply_lamination(self, lamination):
        geometric = [lamination(self.inverse_index_map[index]) if index in self.inverse_index_map else 0 for index in self.target_triangulation.indices]
        return self.target_triangulation(geometric)  # Have to promote.
    
    def apply_homology(self, homology_class):
        algebraic = [homology_class(self.inverse_label_map[index]) if index in self.inverse_label_map else 0 for index in self.target_triangulation.indices]
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def pl_action(self, multicurve):
        action = np.array([[1 if i in self.index_map and j == self.index_map[i] else 0 for i in range(self.source_triangulation.zeta)] for j in range(self.target_triangulation.zeta)], dtype=object)
        condition = np.array([[0] * self.zeta], dtype=object)
        return curver.kernel.PartialLinearFunction(action, condition)

