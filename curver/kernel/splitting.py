
''' A module for building and manipulating splitting sequences. '''

from collections import defaultdict

import curver

class SplittingSequence:  # pylint: disable=too-few-public-methods
    ''' This represents a splitting sequence.
    
    This is an encoding which splits a lamination open along its large branches,
    ensuring that the lamination is a bipod or empty in every triangle. '''
    def __init__(self, lamination, refine, puncture, preperiodic, periodic):  # pylint: disable=too-many-arguments
        self.lamination = lamination
        self.refine = refine
        self.puncture = puncture
        self.preperiodic = preperiodic
        self.periodic = periodic
        self.image = self.preperiodic(self.puncture(self.refine(self.lamination)))
    
    @classmethod
    def from_lamination(cls, lamination, mapping_class):
        ''' Return a splitting sequence from a projectively invariant lamination.
        
        This is the encoding obtained by flipping edges to repeatedly split
        the branches of the corresponding train track with maximal weight
        until you reach a projectively periodic sequence (with the required
        dilatation if given).
        
        Each entry of self.geometric must be an Integer or a RealAlgebraic (over
        the same RealNumberField). '''
        
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
        
        starting_lamination = lamination  # Remember for later.
        refine = lamination.triangulation.id_encoding()
        while True:
            triangulation = lamination.triangulation  # Shorthand.
            if not all(lamination):  # Lamination is not filling.
                null_index = next(index for index in triangulation.indices if lamination(index) == 0)
                curve = triangulation.edge_curve(null_index)
                disjoint_curve = refine.inverse()(curve)  # Pull back.
                assert disjoint_curve and mapping_class(disjoint_curve) == disjoint_curve
                raise ValueError(f'Lamination is not filling, it is disjoint from {disjoint_curve}')
            
            refined_lamination = lamination  # Remember this point in case we later need to refine further.
            
            # Compute the puncture map.
            tripods = set(triangle for triangle in triangulation if all(lamination.dual_weight(side) > 0 for side in triangle))
            internal_sides = set(
                side
                for triangle in tripods
                for side in triangle
                if triangulation.triangle_lookup[~side] in tripods
                and lamination.left_weight(side) == lamination.right_weight(~side)  # Central split.
                and lamination.right_weight(side) == lamination.left_weight(~side)
                )
            boundary = set(
                side
                for triangle in tripods
                for side in triangle
                if side not in internal_sides
                )
            
            # Grab one triangle from each collection to puncture.
            classes = curver.kernel.UnionFind(tripods)
            try:
                for side in internal_sides:
                    if side.sign() == +1:
                        classes.merge(triangulation.triangle_lookup[side], triangulation.triangle_lookup[~side])
            except ValueError:
                # We tried to merge two tripods in the same class, so they didn't form a polygon.
                # TODO: 3) Determine an invariant curve.
                raise ValueError('Lamination is not filling') from None
            seeds = [cls[0] for cls in classes]
            
            # Compute the flips to cone over the collections via a depth first search.
            flips = []
            stack = [side.index for triangle in seeds for side in triangle]
            while stack:
                side = stack.pop()
                if side in internal_sides:
                    flips.append(side)
                    stack.extend(triangulation.corner_lookup[~side][1:])  # Put the other two sides on the stack.
            
            puncture = triangulation.encode_pachner_1_3(seeds)
            puncture = puncture.target_triangulation.encode(flips[::-1]) * puncture
            
            # puncture = lamination.encode_puncture_to_train_track()  # Could raise a ValueError.
            lamination = puncture(lamination)
            
            encodings = [lamination.triangulation.id_encoding()]
            laminations = dict()  # i |--> L_i.
            seen = defaultdict(list)  # hash |--> [i]  # This is a dict taking the hash of each lamination to the index where we saw it.
            
            # Maximal split.
            while True:
                # Save lamination.
                laminations[len(encodings)] = lamination
                seen[projective_hash(lamination)].append(len(encodings))
                
                # Split all of the maximal branches
                move = lamination.triangulation.encode_multiflip(curver.kernel.utilities.maxes(lamination.triangulation.positive_edges, key=lamination))
                encodings.append(move)
                lamination = move(lamination)
                
                if not all(lamination):  # We did a central split and so have found a refinement.
                    break
                
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
                        return cls(starting_lamination, refine, puncture, preperiodic, periodic)
            
            # We've encountered a refinement.
            # One of the zero weight edges will guide a refinement.
            # TODO: 4) Can we use all of the zero weight edges simultaneously?
            null_index = next(index for index in lamination.triangulation.indices if lamination(index) == 0)
            arc = lamination.triangulation.edge_arc(null_index)
            # Pull back to the punctured triangulation.
            for encoding in reversed(encodings):
                arc = encoding.inverse()(arc)
            
            end, _ = [label for label in arc.triangulation.labels if arc.dual_weight(label) < 0]
            
            # Jump back to before doing the puncture.
            lamination = refined_lamination
            
            while True:
                _, _, c, d, _ = lamination.triangulation.square(end)
                ad, bd, cd, dd, _ = [lamination.dual_weight(edgy) for edgy in lamination.triangulation.square(end)]
                assert ad > 0 and bd > 0, f'No cusp on edge {end}'
                
                move = lamination.triangulation.encode_flip(end)
                refine = move * refine
                lamination = move(lamination)
                
                # Follow the cusp.
                if ad < dd:  # Right split.
                    assert bd > cd
                    end = c
                elif ad > dd:  # Left split.
                    assert bd < cd
                    end = d
                else:  # ad == dd  # Central split.
                    assert bd == cd
                    break
        
        raise RuntimeError('Unreachable code')

