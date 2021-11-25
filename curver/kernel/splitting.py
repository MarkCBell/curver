
''' A module for building and manipulating splitting sequences. '''

import curver

# These rely on knowing exactly how edge flips work.
H_MAPPING = [0, 4, 2, 5, 3, 1]
V_MAPPING = [5, 1, 4, 3, 2, 0]

INV_H_MAPPING = [H_MAPPING.index(i) for i in range(6)]
INV_V_MAPPING = [V_MAPPING.index(i) for i in range(6)]

class SplittingSequence:  # pylint: disable=too-few-public-methods
    ''' This represents a splitting sequence.
    
    This is an encoding which splits a lamination open along its large branches,
    ensuring that the lamination is a bipod or empty in every triangle. '''
    def __init__(self, mapping_class):
        ''' Return the splitting sequence of a pseudo-Anosov mapping class.
        
        This is the encoding obtained by flipping edges to repeatedly split
        the branches of the corresponding train track with maximal weight
        until you reach a projectively periodic sequence.
        
        Raises a ValueError if mapping_class is not pseudo-Anosov.
        Sometimes this ValueError describes an invariant multicurve. '''
        
        dilatation, lamination = mapping_class.projectively_invariant_lamination()  # May raise a ValueError if self is not pA.
        assert all(lamination.dual_weight(edge) >= 0 for edge in lamination.triangulation.edges)  # No arcs.
        
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
            
            # Determine where we need to puncture.
            seeds = []
            for group in classes:
                num_tripods = sum(1 for triangle in group if all(lamination.dual_weight(side) > 0 for side in triangle))
                num_bipods = len(group) - num_tripods
                if num_bipods == 0:
                    seeds.append(group[0])
                elif num_bipods == 1:
                    continue
                else:  # num_bipods > 1:
                    # TODO: 3) Determine an invariant curve.
                    raise ValueError('Lamination is not filling')
            
            puncture = lamination.triangulation.encode_pachner_1_3(seeds)
            lamination = puncture(lamination)  # Start again.
            
            # Open up to a train track.
            refine = lamination.triangulation.id_encoding()
            to_open = {side for pair in pairs for side in pair}
            while to_open:
                edge = next(side for side in to_open if lamination.left_weight(side) > 0 and lamination.right_weight(side) > 0 and lamination.dual_weight(side) == 0)
                
                while True:
                    tetra = lamination.triangulation.square(edge) + [~edge]
                    ad, bd, _, dd, ed, _ = [lamination.dual_weight(side) for side in tetra]
                    assert ad > 0 and bd > 0 and ed == 0
                    
                    move = lamination.triangulation.encode_multiflip([+edge])
                    refine = move * refine
                    lamination = move(lamination)
                    
                    if ad == dd:  # No chord.
                        assert ~edge in to_open
                        to_open.remove(edge)
                        to_open.remove(~edge)
                        break
                    
                    mapping = dict()
                    for i, j in enumerate(H_MAPPING if ad > dd else V_MAPPING):
                        mapping[tetra[i]] = tetra[j]
                    
                    to_open = set(mapping.get(side, side) for side in to_open)
                    edge = mapping[edge]
            
            assert all(sum(1 for side in triangle if lamination.dual_weight(side) > 0) == 2 for triangle in lamination.triangulation)
            
            # Compute the preperiodic map.
            # We use Floyd's tortoise and hare algorithm to detect a (projective) cycle.
            preperiodic = lamination.triangulation.id_encoding()
            hare = lamination  # Lamination is our tortoise.
            while True:
                # Split all of the maximal branches
                move = lamination.triangulation.encode_multiflip(curver.kernel.utilities.maxes(lamination.triangulation.positive_edges, key=lamination))
                preperiodic = move * preperiodic
                lamination = move(lamination)
                
                if not all(lamination):  # We did a central split and so have found a refinement.
                    break
                
                for _ in range(2):
                    if all(hare):  # Careful, the move may only be well defined when there isn't a refinement.
                        move = hare.triangulation.encode_multiflip(curver.kernel.utilities.maxes(hare.triangulation.positive_edges, key=hare))
                        hare = move(hare)
                
                if lamination.is_projectively_isometric_to(hare):  # Tortoise in now inside the loop.
                    break
        
            if not all(lamination):  # We've encountered a refinement.
                # Work backwards to find more pairs to start with.
                pairs = []
                for move in preperiodic * refine:
                    if isinstance(move, curver.kernel.Isometry):  # Skip over the id_encoding.
                        assert move.is_identity()
                        continue
                    
                    assert isinstance(move, curver.kernel.MultiEdgeFlip)
                    
                    lamination = move.inverse()(lamination)  # Pull back lamination.
                    
                    # Determine the mapping to pull back the pairs.
                    mapping = dict()
                    new_pairs = []
                    for edge in move.edges:
                        assert edge.sign() == +1
                        
                        tetra = lamination.triangulation.square(edge) + [~edge]
                        ad, _, _, dd, _, _ = [lamination.dual_weight(side) for side in tetra]
                        if ad == dd:  # No chord.
                            new_pairs.append({edge, ~edge})
                            continue
                        
                        for i, j in enumerate(INV_H_MAPPING if ad > dd else INV_V_MAPPING):
                            mapping[tetra[i]] = tetra[j]
                    
                    pairs = [set(mapping.get(x, x) for x in pair) for pair in pairs] + new_pairs
                
                continue  # Start again.
            
            # Compute the periodic map.
            cycle_start = lamination  # Remember where we started.
            periodic = lamination.triangulation.id_encoding()
            while True:
                move = lamination.triangulation.encode_multiflip(curver.kernel.utilities.maxes(lamination.triangulation.positive_edges, key=lamination))
                periodic = move * periodic
                lamination = move(lamination)
                
                if lamination.weight() * dilatation > cycle_start.weight():  # pylint: disable=no-else-continue
                    continue  # Haven't gone far enough around periodic.
                elif lamination.weight() * dilatation == cycle_start.weight():
                    for isometry in (lamination * cycle_start.weight()).isometries_to(cycle_start * lamination.weight()):
                        closed = (isometry.encode() * periodic).inverse()  # Don't forget to close up
                        # Currently we can't use **.homology_matrix() since that isn't defined for encodings.
                        if all(closed(preperiodic(refine(puncture(hc)))) == preperiodic(refine(puncture(mapping_class(hc)))) for hc in starting_lamination.triangulation.edge_homologies()):
                            self.puncture = puncture
                            self.refine = refine
                            self.preperiodic = preperiodic
                            self.periodic = closed
                            self.lamination = starting_lamination
                            self.stable_lamination = self.preperiodic(self.refine(self.puncture(self.lamination)))
                            self.punctures = self.preperiodic(self.refine(self.puncture(self.lamination.triangulation.peripheral_multicurve())))
                            return
                else:  # lamination.weight() * dilatation < cycle_start.weight()
                    raise RuntimeError(f'No cycle with dilatation {dilatation}')
        
        raise RuntimeError('Unreachable code')

