
''' A module for building and manipulating splitting sequences. '''

from collections import Counter, defaultdict
import numpy as np

import curver

class SplittingSequence(object):
    ''' This represents a splitting sequence.
    
    This is an encoding which splits a lamination open along its large branches,
    ensuring that the lamination is a bipod or empty in every triangle. '''
    def __init__(self, lamination, puncture, preperiodic, periodic):
        self.puncture = puncture
        self.preperiodic = preperiodic
        self.periodic = periodic
        
        self.lamination = lamination
        self.punctured_lamination = self.puncture(self.lamination)
        self.periodic_lamination = self.preperiodic(self.punctured_lamination)
        
        # Save some useful triangulations.
        self.triangulation = self.puncture.source_triangulation
        self.punctured_triangulation = self.preperiodic.source_triangulation  # = self.puncture.target_triangulation.
        self.periodic_triangulation = self.periodic.source_triangulation  # = self.preperiodic.target_triangulation.
        
        # The boundary of lamination.
        self.periodic_boundary = self.periodic_triangulation([2 if weight > 0 else 0 for weight in self.periodic_lamination])
        self.punctured_boundary = self.preperiodic.inverse()(self.periodic_boundary)
        
        assert self.periodic(self.periodic_boundary) == self.periodic_boundary  # TODO: Make a unittest from this.
    
    @classmethod
    def from_lamination(cls, starting_lamination, mapping_class):
        ''' Return a splitting sequence from a projectively invariant starting_lamination.
        
        This is the encoding obtained by flipping edges to repeatedly split
        the branches of the corresponding train track with maximal weight
        until you reach a projectively periodic sequence (with the required
        dilatation if given).
        
        Each entry of self.geometric must be an Integer or a RealAlgebraic (over
        the same RealNumberField). '''
        
        # In this method we use Lamination.projective_hash to store the laminations
        # we encounter efficiently and so avoid a quadratic algorithm. This currently
        # only ever uses the default precision HASH_DENOMINATOR. At some point this
        # should change dynamically depending on the algebraic numbers involved in
        # this lamination.
        
        def projective_hash(L, precision=30):
            w = L.weight()
            PL = [x * 10**precision // w for x in L]
            
            # We'll try to preserve as much of the structure as possible to try to reduce hash collisions.
            # In this version we'll store the sorted, cyclically ordered, triangles.
            triples = [tuple(PL[edge.index] for edge in triangle) for triangle in L.triangulation]
            return tuple(sorted([min(triple[i:] + triple[:i] for i in range(len(triple))) for triple in triples]))
        
        def projectivise(L):
            w = L.weight()
            return L * (1 / w)
        
        assert projectivise(mapping_class(starting_lamination)) == projectivise(starting_lamination)
        assert all(weight >= 0 for weight in starting_lamination)
        assert all(starting_lamination.dual_weight(edge) >= 0 for edge in starting_lamination.triangulation.edges)
        
        puncture = starting_lamination.triangulation.encode_pachner_1_3()
        lamination = puncture(starting_lamination)
        
        encodings = [lamination.triangulation.id_encoding()]
        laminations = dict()  # i |--> L_i.
        seen = defaultdict(list)  # hash |--> [i]  # This is a dict taking the hash of each lamination to the index where we saw it.
        
        while True:
            # Save lamination.
            laminations[len(encodings)] = lamination
            seen[projective_hash(lamination)].append(len(encodings))
            
            # Remove all of the obvious boundary.
            non_peripheral_boundary = lamination.triangulation([2 if lamination(index) > 0 else 0 for index in lamination.triangulation.indices]).non_peripheral()
            if non_peripheral_boundary:  # is not empty.
                move = non_peripheral_boundary.crush()
                encodings.append(move)
                lamination = move(lamination)
                
                # Reset dictionaries ...
                laminations = dict()
                seen = defaultdict(list)
                # ... and re-save lamination.
                laminations[len(encodings)] = lamination
                seen[projective_hash(lamination)].append(len(encodings))
            
            # Split all of the maximal branches
            max_weight = max(lamination)
            max_edges = [edge for edge in lamination.triangulation.positive_edges if lamination(edge) == max_weight]
            for edge in max_edges:
                move = lamination.triangulation.encode_flip(edge)
                encodings.append(move)
                lamination = move(lamination)
            
            # Check if lamination now (projectively) matches a lamination we've already seen.
            for index in seen.get(projective_hash(lamination), []):
                weight = lamination.weight()
                old_lamination = laminations[index]
                old_weight = old_lamination.weight()
                
                for isometry in lamination.triangulation.isometries_to(old_lamination.triangulation):
                    isom_e = isometry.encode()
                    if isom_e(lamination) * old_weight == weight * old_lamination:  # A projective isometry.
                        preperiodic = curver.kernel.Encoding([move for item in reversed(encodings[:index]) for move in item]).promote()
                        open_periodic = curver.kernel.Encoding([move for item in reversed(encodings[index:]) for move in item]).promote()
                        periodic = isom_e * open_periodic
                        # We really should only return for the correct isometry.
                        # This should be determined by mapping_class.homology_matrix() and periodic.homology_matrix().
                        # if np.array_equal((preperiodic * mapping_class).homology_matrix(), (periodic * preperiodic).homology_matrix()):  # if isometry gives correct map.
                        return cls(starting_lamination, puncture, preperiodic, periodic)
        
        raise RuntimeError('Unreachable code.')
    
    def essential_punctured_boundary(self):
        ''' Return the MultiCurve consisiting of the components of punctured_boundary that are essential in self.triangulation.
        
        We do this by looking for components of self.punctured_boundary that do not bound a disk or punctured disk.
        We test for this by crushing along a candidate curve and checking whether all the components that correspond to it either:
        have genus or have more than one real vertex. '''
        
        def is_essential(curve):
            ''' Return whether the given curve is essential in the original surface. '''
            
            if curve.is_peripheral():  # x is obviously peripheral or null-homotopic in the original surface.
                return False
            
            crush = curve.crush()
            lift = crush.inverse()
            T = crush.target_triangulation
            T_components = set(c.containing_components().pop() for c in T([2] * T.zeta).components() if lift(c) == curve)  # The components of T that contain the curve vertices.
            assert len(T_components) <= 2
            
            S = T.surface()
            real_punctures = crush(self.puncture(self.triangulation([2] * self.triangulation.zeta)))  # Find the punctures of T that are real.
            num_real_punctures = Counter(component.containing_components().pop() for component in real_punctures.components())
            
            return not any(S[component].g == 0 and num_real_punctures[component] <= 1 for component in T_components)
        
        return self.punctured_triangulation.disjoint_sum([curve for curve in self.punctured_boundary.components() if is_essential(curve)])

