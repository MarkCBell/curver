
''' A module for building and manipulating splitting sequences. '''

from collections import defaultdict

import curver

class SplittingSequence:  # pylint: disable=too-few-public-methods
    ''' This represents a splitting sequence.
    
    This is an encoding which splits a lamination open along its large branches,
    ensuring that the lamination is a bipod or empty in every triangle. '''
    def __init__(self, lamination, puncture, refine, preperiodic, periodic):  # pylint: disable=too-many-arguments
        self.lamination = lamination
        self.puncture = puncture
        self.refine = refine
        self.preperiodic = preperiodic
        self.periodic = periodic
        self.image = self.preperiodic(self.refine(self.puncture(self.lamination)))
    
    @classmethod
    def from_lamination(cls, lamination, mapping_class):
        ''' Return a splitting sequence from a projectively invariant lamination.
        
        This is the encoding obtained by flipping edges to repeatedly split
        the branches of the corresponding train track with maximal weight
        until you reach a projectively periodic sequence.
        
        Each entry of self.geometric must be an Integer or a RealAlgebraic (over
        the same RealNumberField).
        
        Raises a ValueError describing a disjoint curve if the lamination is not filling. '''
        
        def projective_hash(L):
            ''' Return a hash key for L that is scale and isometry invariant. '''
            weight = L.weight()
            PL = [entry / weight for entry in L]
            
            # We'll try to preserve as much of the structure as possible to try to reduce hash collisions.
            # In this version we'll store the sorted, cyclically ordered, triangles.
            triples = [tuple(PL[edge.index] for edge in triangle) for triangle in L.triangulation]
            return tuple(sorted([min(triple[i:] + triple[:i] for i in range(len(triple))) for triple in triples]))
        
        assert all(lamination.dual_weight(edge) >= 0 for edge in lamination.triangulation.edges)  # No arcs.
        image = mapping_class(lamination)
        assert (image * lamination.weight()).is_isometric_to(lamination * image.weight())
        
        if not all(lamination):
            null_index = next(edge for edge in lamination.triangulation.positive_edges if lamination(edge) == 0)
            curve = lamination.triangulation.edge_curve(null_index)
            assert curve.intersection(lamination) == 0
            assert mapping_class(curve) == curve
            raise ValueError(f'Lamination is not filling, it is disjoint from {curve}')
        
        starting_lamination = lamination  # Remember where we started.
        pairs = []  # List of pairs {side, side} whose corridors are connected.
        while True:
            lamination = starting_lamination  # Reset.
            classes = curver.kernel.UnionFind(lamination.triangulation)
            for side, sidy in pairs:
                try:
                    classes.merge(lamination.triangulation.triangle_lookup[side], lamination.triangulation.triangle_lookup[sidy])
                except ValueError:
                    # We tried to merge two tripods in the same class, so they didn't form a polygon.
                    # TODO: 3) Determine an invariant curve.
                    raise ValueError('Lamination is not filling') from None
            
            seeds = []
            for group in classes:
                num_tripods = sum(1 for triangle in group if all(lamination.dual_weight(side) > 0 for side in triangle))
                num_bipods = len(group) - num_tripods
                if num_bipods == 0:
                    seeds.append(group[0])
                elif num_bipods == 1:
                    continue
                else:  # num_bipods > 1:
                    raise ValueError('Lamination is not filling')
            
            puncture = lamination.triangulation.encode_pachner_1_3(seeds)
            lamination = puncture(lamination)  # Start again.
            
            refine = lamination.triangulation.id_encoding()
            to_open = {side for pair in pairs for side in pair}
            while to_open:
                edge = next(side for side in to_open if lamination.left_weight(side) > 0 and lamination.right_weight(side) > 0 and lamination.dual_weight(side) == 0)
                
                while True:
                    a, b, c, d, _ = lamination.triangulation.square(edge)
                    ad, bd, cd, dd, ed = [lamination.dual_weight(side) for side in lamination.triangulation.square(edge)]
                    assert ad > 0
                    assert bd > 0
                    assert ed == 0
                    
                    move = lamination.triangulation.encode_flip(+edge)
                    refine = move * refine
                    lamination = move(lamination)
                    
                    mapping = dict()
                    if ad > dd:  # Horizontal chord.
                        mapping[edge] = d
                        mapping[d] = ~edge
                        mapping[~edge] = b
                        mapping[b] = edge
                    elif ad < dd:  # Vertical chord.
                        mapping[edge] = c
                        mapping[c] = edge
                        mapping[~edge] = a
                        mapping[a] = ~edge
                    else:  # ad == dd  # No chord.
                        to_open.remove(edge)
                        to_open.remove(~edge)
                        break
                    
                    to_open = set(mapping.get(side, side) for side in to_open)
                    edge = mapping[edge]
            
            assert all(sum(1 for side in triangle if lamination.dual_weight(side) > 0) == 2 for triangle in lamination.triangulation)
            
            encoding = lamination.triangulation.id_encoding()
            laminations = [lamination]
            proj_indices = defaultdict(list)  # hash |--> [i]  # This is a dict taking the hash of each lamination to the index where we saw it.
            
            # Maximal split.
            while True:
                # Split all of the maximal branches
                move = lamination.triangulation.encode_multiflip(curver.kernel.utilities.maxes(lamination.triangulation.positive_edges, key=lamination))
                encoding = move * encoding
                lamination = move(lamination)
                
                if not all(lamination):  # We did a central split and so have found a refinement.
                    break
                
                # Check if lamination now (projectively) matches a lamination we've already seen.
                for index in proj_indices.get(projective_hash(lamination), []):
                    old_lamination = laminations[index]
                    
                    for isometry in (lamination * old_lamination.weight()).isometries_to(old_lamination * lamination.weight()):
                        isom_e = isometry.encode()
                        preperiodic = encoding[-index:]
                        periodic = isom_e * encoding[:-index]  # Don't forget to close up
                        # We really should only return for the correct isometry.
                        # This should be determined by mapping_class.homology_matrix() and periodic.homology_matrix().
                        # if np.array_equal((preperiodic * mapping_class).homology_matrix(), (periodic * preperiodic).homology_matrix()):  # if isometry gives correct map.
                        return cls(starting_lamination, puncture, refine, preperiodic, periodic)
                
                proj_indices[projective_hash(lamination)].append(len(laminations))
                laminations.append(lamination)
            
            # We've encountered a refinement.
            
            # Work backwards to find more pairs to start with.
            pairs = []
            for move in encoding * refine:
                if isinstance(move, curver.kernel.Isometry):  # Skip over the id_encoding.
                    continue
                
                assert isinstance(move, (curver.kernel.EdgeFlip, curver.kernel.MultiEdgeFlip))
                mapping = dict()
                new_pairs = []
                for edge in [move.edge] if isinstance(move, curver.kernel.EdgeFlip) else move.edges:
                    assert edge.sign() == +1
                    
                    a, b, c, d, e = move.target_triangulation.square(edge)
                    ad, _, _, dd, _ = [lamination.dual_weight(side) for side in lamination.triangulation.square(edge)]
                    if ad > dd:  # Horizontal chord.
                        mapping[edge] = d
                        mapping[d] = edge
                        mapping[b] = ~edge
                        mapping[~edge] = b
                    elif ad < dd:  # Vertical chord.
                        mapping[edge] = c
                        mapping[c] = ~edge
                        mapping[~edge] = a
                        mapping[a] = edge
                    else:  # ad == dd  # No chord.
                        new_pairs.append({e, ~e})
                
                lamination = move.inverse()(lamination)
                pairs = [set(mapping.get(x, x) for x in pair) for pair in pairs] + new_pairs
        
        raise RuntimeError('Unreachable code')

