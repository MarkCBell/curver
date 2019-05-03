
import heapq
from collections import defaultdict

import curver
from curver.kernel.moves import Move

def puncture():
    pass

class SplittingSequence(object):
    def __init__(self, x):
        self.x = x
    
    @classmethod
    def from_lamination(cls, lamination, mapping_class):
        ''' Return a splitting sequence from a projectively invariant lamination.
        
        This is the encoding obtained by flipping edges to repeatedly split
        the branches of the corresponding train track with maximal weight
        until you reach a projectively periodic sequence (with the required
        dilatation if given).
        
        Assumes (and checks) that this lamination is filling.
        
        Each entry of self.geometric must be an Integer or a RealAlgebraic (over
        the same RealNumberField). '''
        
        # In this method we use Lamination.projective_hash to store the laminations
        # we encounter efficiently and so avoid a quadratic algorithm. This currently
        # only ever uses the default precision HASH_DENOMINATOR. At some point this
        # should change dynamically depending on the algebraic numbers involved in
        # this lamination.
        
        def projective_hash(L, precision=30):
            w = self.weight()
            PL = [x * 10**precision // w for x in L]
            
            # We'll try to preserve as much of the structure as possible to try to reduce hash collisions.
            # In this version we'll store the sorted, cyclically ordered, triangles.
            triples = [tuple(PL[edge.index] for edge in triangle) for triangle in L.triangulation]
            return tuple(sorted([min(triple[i:] + triple[:i] for i in range(len(triple))) for triple in triples]))
        
        def projectivise(L):
            w = L.weight()
            return L * (1 / w)
        
        assert projectivise(mapping_class(lamination)) == projectivise(lamination)
        
        # This is a dict taking the hash of each lamination to the index where we saw it.
        encodings = []
        laminations = dict()  # i |--> L_i.
        seen = defaultdict(list)  # hash |--> [i]
        
        def store_move(move, lamination):
            encodings.append(move)
            return move(lamination)
        
        def store_lamination(lamination):
            laminations[len(encodings)] = lamination
            seen[projective_hash(lamination)].append(len(encodings))
        
        # Remove all of the obvious places.
        for edge_index in range(lamination.zeta):
            if lamination(edge_index) == 0:
                # TODO: Choose crush curve best.
                lamination = store_move(lamination.triangulation.edge_curve(edge_index).crush(), lamination)
        
        # Puncture all the triangles where the lamination is a tripod.
        lamination = store_move(puncture_tripods(lamination), lamination)
        
        store_lamination(lamination)
        
        while True:
            # Get the indices of the largest weight edges since we will be changing lamination as we iterate over the edges.
            max_weight = max(lamination)
            max_indices = [index for index, weight in enumerate(lamination) if weight == max_weight]
            
            # Flip all of these edges.
            for flip_index in max_indices:
                lamination = store_move(lamination.triangulation.encode_flip(flip_index), lamination)
                # Check if we have created any edges of weight 0. Of course it is enough to just check the flip_index.
                if lamination(flip_index) == 0:
                    lamination = store_move(lamination.triangulation.edge_curve(flip_index).crush(), lamination)
            
            # Check if lamination now (projectively) matches a lamination we've already seen.
            for index in seen.get(projective_hash(lamination), []):
                weight = lamination.weight()
                old_lamination = laminations[index]
                old_weight = old_lamination.weight()
                
                for isometry in lamination.triangulation.isometries_to(old_lamination.triangulation):
                    if isometry.encode()(lamination) * old_weight == weight * old_lamination:  # A projective isometry.
                        encoding = flipper.kernel.Encoding([move for item in reversed(encodings) for move in item])
                        return cls(encoding, isometry, index, old_lamination)
            
            store_lamination(lamination)
        
        raise RuntimeError('Unreachable code.')

