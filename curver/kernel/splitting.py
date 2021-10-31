
''' A module for building and manipulating splitting sequences. '''

from collections import defaultdict

import curver

class SplittingSequence:  # pylint: disable=too-few-public-methods
    ''' This represents a splitting sequence.
    
    This is an encoding which splits a lamination open along its large branches,
    ensuring that the lamination is a bipod or empty in every triangle. '''
    def __init__(self, puncture, preperiodic, periodic, lamination, punctured_lamination=None, periodic_lamination=None):  # pylint: disable=too-many-arguments
        self.puncture = puncture
        self.preperiodic = preperiodic
        self.periodic = periodic
        
        self.lamination = lamination
        self.punctured_lamination = punctured_lamination if punctured_lamination is not None else self.puncture(self.lamination)
        self.periodic_lamination = periodic_lamination if periodic_lamination is not None else self.preperiodic(self.punctured_lamination)
        
        # Save some useful triangulations.
        self.triangulation = self.lamination.triangulation
        self.punctured_triangulation = self.punctured_lamination.triangulation
        self.periodic_triangulation = self.periodic_lamination.triangulation
        
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
        # we encounter efficiently and so avoid a quadratic algorithm.
        
        def projective_hash(L):
            PL = projectivise(L)
            
            # We'll try to preserve as much of the structure as possible to try to reduce hash collisions.
            # In this version we'll store the sorted, cyclically ordered, triangles.
            triples = [tuple(PL[edge.index] for edge in triangle) for triangle in L.triangulation]
            return tuple(sorted([min(triple[i:] + triple[:i] for i in range(len(triple))) for triple in triples]))
        
        def projectivise(L):
            weight = L.weight()
            return [entry / weight for entry in L]
        
        assert projectivise(mapping_class(starting_lamination)) == projectivise(starting_lamination)
        assert all(starting_lamination.dual_weight(edge) >= 0 for edge in starting_lamination.triangulation.edges)
        
        while True:
            if any(weight == 0 for weight in starting_lamination):
                raise ValueError('Lamination is not filling')
            
            puncture = starting_lamination.encode_puncture_to_train_track()
            lamination = punctured_lamination = puncture(starting_lamination)
            
            encodings = [lamination.triangulation.id_encoding()]
            laminations = dict()  # i |--> L_i.
            seen = defaultdict(list)  # hash |--> [i]  # This is a dict taking the hash of each lamination to the index where we saw it.
            
            while True:
                # Save lamination.
                laminations[len(encodings)] = lamination
                seen[projective_hash(lamination)].append(len(encodings))
                
                assert all(sum(1 if lamination.dual_weight(edge) > 0 else 0 for edge in triangle) in (0, 1, 2) for triangle in lamination.triangulation)
                
                # Write down the multiarc which is obviously disjoint from lamination.
                disjoint_multiarc = lamination.triangulation([-1 if lamination(index) == 0 else 0 for index in lamination.triangulation.indices])
                if disjoint_multiarc:
                    preperiodic = curver.kernel.Encoding([move for item in reversed(encodings) for move in item])  # Don't bother promoting.
                    original_multiarc = preperiodic.inverse()(disjoint_multiarc)  # The disjoint arcs back on punctured_lamination.
                    
                    punctured_edges = [edge for edge in punctured_lamination.triangulation.edges if original_multiarc.dual_weight(edge) < 0]
                    assert all(edge.index < starting_lamination.zeta for edge in punctured_edges)
                    
                    edge = punctured_edges[0]
                    while True:
                        corner = starting_lamination.triangulation.corner_lookup[~edge]
                        right = starting_lamination.right_weight(edge)
                        left = starting_lamination.left_weight(~edge)
                        
                        # We could remember all of the moves that we have applied to the starting_lamination.
                        move = starting_lamination.triangulation.encode_flip(edge)
                        starting_lamination = move(starting_lamination)
                        
                        if left < right:
                            edge = corner[1]
                        elif left > right:
                            edge = corner[2]
                        else:  # if left == right:
                            break
                    
                    break
                
                assert all(sum(1 if lamination.dual_weight(edge) > 0 else 0 for edge in triangle) in (0, 2) for triangle in lamination.triangulation)
                
                # Split all of the maximal branches
                move = lamination.triangulation.encode_multiflip(curver.kernel.utilities.maxes(lamination.triangulation.positive_edges, key=lamination))
                encodings.append(move)
                lamination = move(lamination)
                
                # Check if lamination now (projectively) matches a lamination we've already seen.
                for index in seen.get(projective_hash(lamination), []):
                    old_lamination = laminations[index]
                    
                    scaled_old_lamination = old_lamination * lamination.weight()
                    scaled_lamination = lamination * old_lamination.weight()
                    
                    for isometry in scaled_lamination.isometries_to(scaled_old_lamination):
                        isom_e = isometry.encode()
                        preperiodic = curver.kernel.Encoding([move for item in reversed(encodings[:index]) for move in item]).promote()
                        open_periodic = curver.kernel.Encoding([move for item in reversed(encodings[index:]) for move in item]).promote()
                        periodic = isom_e * open_periodic
                        # We really should only return for the correct isometry.
                        # This should be determined by mapping_class.homology_matrix() and periodic.homology_matrix().
                        # if np.array_equal((preperiodic * mapping_class).homology_matrix(), (periodic * preperiodic).homology_matrix()):  # if isometry gives correct map.
                        return cls(puncture, preperiodic, periodic, starting_lamination, punctured_lamination, old_lamination)
        
        raise RuntimeError('Unreachable code')

