
from collections import defaultdict

import curver

class SplittingSequence(object):
    def __init__(self, preperiodic, periodic, lamination):
        self.preperiodic = preperiodic
        self.periodic = periodic
        self.lamination = lamination
        
        self.boundary = self.lamination.triangulation([2 if weight > 0 else 0 for weight in self.lamination])
    
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
            w = L.weight()
            PL = [x * 10**precision // w for x in L]
            
            # We'll try to preserve as much of the structure as possible to try to reduce hash collisions.
            # In this version we'll store the sorted, cyclically ordered, triangles.
            triples = [tuple(PL[edge.index] for edge in triangle) for triangle in L.triangulation]
            return tuple(sorted([min(triple[i:] + triple[:i] for i in range(len(triple))) for triple in triples]))
        
        def projectivise(L):
            w = L.weight()
            return L * (1 / w)
        
        assert projectivise(mapping_class(lamination)) == projectivise(lamination)
        assert all(weight >= 0 for weight in lamination)
        assert all(lamination.dual_weight(edge) >= 0 for edge in lamination.triangulation.edges)
        
        # Puncture all the triangles where the lamination is a tripod.
        zeta = lamination.zeta
        geometric = list(lamination)
        
        triangles = []
        for triangle in lamination.triangulation:
            a, b, c = triangle.edges
            if all(lamination.dual_weight(edge) > 0 for edge in triangle):  # Is tripod.
                s, t, u = curver.kernel.Edge(zeta), curver.kernel.Edge(zeta+1), curver.kernel.Edge(zeta+2)  # New edges.
                triangles.extend([curver.kernel.Triangle([a, ~u, t]), curver.kernel.Triangle([b, ~s, u]), curver.kernel.Triangle([c, ~t, s])])
                geometric.extend([lamination.dual_weight(a), lamination.dual_weight(b), lamination.dual_weight(c)])
                
                zeta += 3
            else:
                triangles.append(curver.kernel.Triangle([a, b, c]))
        
        lamination = curver.kernel.Triangulation(triangles)(geometric)
        
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
                
                # Re-save lamination.
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
                        return cls(preperiodic, periodic, old_lamination)
        
        raise RuntimeError('Unreachable code.')

