
''' A module for representing more advanced ways of changing triangulations. '''

import numpy as np

import curver
from curver.kernel.moves import Move  # Special import needed for subclassing.

class Crush(Move):
    ''' This represents the effect of crushing along a curve. '''
    def __init__(self, source_triangulation, target_triangulation, curve):
        super().__init__(source_triangulation, target_triangulation)
        
        assert isinstance(curve, curver.kernel.Curve)
        assert not curve.is_peripheral() and curve.is_short()
        assert curve.triangulation == self.source_triangulation
        
        self.curve = curve
    
    def __str__(self):
        return 'Crush ' + str(self.curve)
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.curve == other.curve
    
    def apply_lamination(self, lamination):
        geometric = list(lamination)
        
        # Get some edges.
        a = self.curve.parallel()
        _, b, e = self.source_triangulation.corner_lookup[a]
        
        v = self.curve.triangulation.vertex_lookup[a]  # = self.triangulation.vertex_lookup[~a].
        v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
        around_v = curver.kernel.utilities.minimal((lamination.left_weight(edgy) for edgy in v_edges), lower_bound=0)
        out_v = sum(max(-lamination.left_weight(edge), 0) for edge in v_edges) + sum(max(-lamination(edge), 0) for edge in v_edges[1:])
        # around_v > 0 ==> out_v == 0; out_v > 0 ==> around_v == 0.
        twisting = curver.kernel.utilities.minimal((lamination.left_weight(edgy) - around_v for edgy in v_edges[1:-1]), lower_bound=0)
        
        # We could have initially removed the twisting via the fact that:
        # twisting == abs(self.curve.slope(lamination) * lamination(a))
        
        # Computing around_v and twisting can be done more efficiently.
        
        # We work by manipulating the side weights around v.
        sides = dict((edge, lamination.left_weight(edge) - (self.curve.left_weight(edge)*twisting + around_v if edge in v_edges and lamination.left_weight(edge) >= 0 else 0)) for edge in self.source_triangulation.edges)
        parallels = dict((edge.index, max(-lamination(edge), 0)) for edge in v_edges)
        
        # TODO: 4) Add comments explaining what is going on in the next two blocks and how the different tightening cases work.
        
        # Tighten to the left.
        drop = max(sides[a], 0) + max(-sides[b], 0)
        for edge in v_edges[1:-1]:
            x, y, z = lamination.triangulation.corner_lookup[edge]
            if sides[x] >= 0 and sides[y] >= 0 and sides[z] >= 0:
                if drop <= sides[x]:
                    sides[x] = sides[x] - drop
                else:  # sides[x] < drop.
                    sides[x], sides[y], drop = sides[x] - drop, sides[y] + sides[x] - drop, sides[x]
            elif sides[x] < 0:
                sides[x], sides[y], drop = sides[x] - drop, sides[y] - drop, 0
            elif sides[y] < 0:
                sides[x] = sides[x] - drop
            else:  # sides[z] < 0.
                if drop <= sides[x]:
                    sides[x] = sides[x] - drop
                elif sides[x] < drop <= sides[x] - sides[z]:
                    parallels[z.index] = parallels[z.index] + (drop - sides[x])
                    sides[x], sides[z], drop = 0, sides[z] + (drop - sides[x]), sides[x]
                else:  # sides[x] - sides[z] < drop:
                    parallels[z.index] = parallels[z.index] - sides[z]
                    sides[x], sides[y], sides[z], drop = sides[x] - sides[z] - drop, sides[y] - (drop - sides[x] + sides[z]), 0, sides[x]
            
            if drop == 0: break  # Stop early.
        
        # Tighten to the right.
        drop = max(-sides[a], 0) + max(sides[b], 0)
        for edge in reversed(v_edges[1:-1]):
            x, y, z = lamination.triangulation.corner_lookup[edge]
            if sides[x] >= 0 and sides[y] >= 0 and sides[z] >= 0:
                if drop <= sides[x]:
                    sides[x] = sides[x] - drop
                else:  # sides[x] < drop.
                    sides[x], sides[z], drop = sides[x] - drop, sides[z] + sides[x] - drop, sides[x]
            elif sides[x] < 0:
                sides[x], sides[z], drop = sides[x] - drop, sides[z] - drop, 0
            elif sides[y] < 0:
                if drop <= sides[x]:
                    sides[x] = sides[x] - drop
                elif sides[x] < drop <= sides[x] - sides[y]:
                    parallels[x.index] = parallels[x.index] + (drop - sides[x])
                    sides[x], sides[y], drop = 0, sides[y] + (drop - sides[x]), sides[x]
                else:  # sides[x] - sides[y] < drop:
                    parallels[x.index] = parallels[x.index] - sides[y]
                    sides[x], sides[y], sides[z], drop = sides[x] - sides[y] - drop, 0, sides[z] - (drop - sides[x] + sides[y]), sides[x]
            else:  # sides[z] < 0.
                sides[x] = sides[x] - drop
            
            if drop == 0: break  # Stop early.
        
        # Now rebuild the intersection.
        for edge in v_edges:
            if edge not in (a, b, e, ~b, ~e):
                x, y, z = lamination.triangulation.corner_lookup[edge]
                if parallels[edge.index] > 0:
                    geometric[edge.index] = -parallels[x.index]
                else:
                    geometric[edge.index] = max(sides[x], 0) + max(sides[y], 0) + max(-sides[z], 0)
                    
                    # Sanity check:
                    x2, y2, z2 = lamination.triangulation.corner_lookup[~edge]
                    assert geometric[edge.index] == max(sides[x2], 0) + max(sides[y2], 0) + max(-sides[z2], 0)
        
        # We have to rebuild the ~e edge separately since it now pairs with ~b.
        x, y, z = lamination.triangulation.corner_lookup[~e]
        if parallels[e.index] + parallels[b.index] + max(-sides[e], 0) > 0:
            geometric[e.index] = -(parallels[e.index] + parallels[b.index] + max(-sides[e], 0))
        else:
            geometric[e.index] = max(sides[x], 0) + max(sides[y], 0) + max(-sides[z], 0)
            
            # Sanity check:
            x2, y2, z2 = lamination.triangulation.corner_lookup[~b]
            assert geometric[e.index] == max(sides[x2], 0) + max(sides[y2], 0) + max(-sides[z2], 0)
        
        # And finally the b edge, which is now paired with e.
        # Since around_v > 0 ==> out_v == 0 & out_v > 0 ==> around_v == 0, this is equivalent to: around_v if around_v > 0 else -out_v
        geometric[b.index] = around_v - out_v
        
        return self.target_triangulation(geometric)  # Have to promote.
    
    def apply_homology(self, homology_class):
        return NotImplemented  # I don't think we ever need this.
    
    def package(self):
        return (self.curve.parallel(), 0)

class LinearTransformation(Move):
    ''' This represents a linear transformation between two triangulations. '''
    def __init__(self, source_triangulation, target_triangulation, matrix):
        super().__init__(source_triangulation, target_triangulation)
        assert matrix.shape == (target_triangulation.zeta, source_triangulation.zeta)
        
        self.matrix = matrix
    
    def __str__(self):
        return 'LT to %s' % self.target_triangulation
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return np.array_equal(self.matrix, other.matrix)
    
    def package(self):
        return (self.target_triangulation.sig(), self.matrix.tolist())
    
    def apply_lamination(self, lamination):
        return self.target_triangulation(self.matrix.dot(lamination.geometric).tolist())
    
    def apply_homology(self, homology_class):
        return NotImplemented  # I don't think we ever need this.

class Lift(LinearTransformation):
    ''' This represents the inverse of crushing along a curve. '''
    def __init__(self, source_triangulation, target_triangulation, matrix):
        super().__init__(source_triangulation, target_triangulation, matrix)
        
        # We need to use super again since we have not found the vertices needed so that we can call self yet.
        apply_lamination = super().apply_lamination
        self.vertices = [vertex for vertex in self.source_triangulation.vertices if apply_lamination(self.source_triangulation.curve_from_cut_sequence(vertex)).is_peripheral()]
    
    def __str__(self):
        return 'Lift to %s' % self.target_triangulation
    
    def apply_lamination(self, lamination):
        assert all(lamination(edge) >= 0 and lamination.left_weight(edge) >= 0 for vertex in self.vertices for edge in vertex)
        
        return super().apply_lamination(lamination)

