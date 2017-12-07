
''' A module for representing more advanced ways of changing triangulations. '''

import curver
from curver.kernel.moves import Move  # Special import needed for subclassing.

class Crush(Move):
    ''' This represents the effect of crushing along a curve. '''
    def __init__(self, source_triangulation, target_triangulation, curve, matrix):
        super(Crush, self).__init__(source_triangulation, target_triangulation)
        
        assert(isinstance(curve, curver.kernel.Curve))
        assert(curve.is_short() and not curve.is_peripheral())
        assert(curve.triangulation == self.source_triangulation)
        
        self.curve = curve
        self.matrix = matrix
    
    def __str__(self):
        return 'Crush ' + str(self.curve)
    def __reduce__(self):
        return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve, self.matrix))
    def __eq__(self, other):
        if isinstance(other, Crush):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation and self.curve == other.curve and self.matrix == other.matrix
        else:
            return NotImplemented
    def __ne__(self, other):
        return not (self == other)
    
    def apply_lamination(self, lamination):
        geometric = list(lamination)
        
        if self.curve.is_isolating():
            # Get some edges.
            a = self.curve.parallel()
            v = self.curve.triangulation.vertex_lookup[a.label]  # = self.triangulation.vertex_lookup[~a.label].
            _, b, e = self.source_triangulation.corner_lookup[a.label]
            
            v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
            around_v = min(max(lamination.side_weight(edge), 0) for edge in v_edges)
            out_v = sum(max(-lamination.side_weight(edge), 0) for edge in v_edges) + sum(max(-lamination(edge), 0) for edge in v_edges[1:])
            # around_v > 0 ==> out_v == 0; out_v > 0 ==> around_v == 0.
            
            # TODO: 4) Add comments explaining what is going on here and how the different cases work.
            # We work by manipulating the dual weights around v.
            twisting = min(max(lamination.side_weight(edge) - around_v, 0) for edge in v_edges[1:-1])
            sides = dict()
            for edge in v_edges:
                if lamination.side_weight(edge) < 0:
                    sides[edge] = lamination.side_weight(edge)
                else:  # lamination.side_weight(edge) >= 0:
                    if edge == a or edge == b:
                        sides[edge] = lamination.side_weight(edge) - around_v
                    elif edge == e:
                        sides[edge] = lamination.side_weight(edge) - 2*twisting - around_v
                    else:
                        sides[edge] = lamination.side_weight(edge) - twisting - around_v
            parallels = dict((edge.index, max(-lamination(edge), 0)) for edge in v_edges)
            
            # Tighten to the left.
            drop = max(sides[a], 0) + max(-sides[b], 0)
            for edge in v_edges[1:-1]:
                x, y, z = lamination.triangulation.corner_lookup[edge.label]
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
            for edge in v_edges[-2:0:-1]:
                x, y, z = lamination.triangulation.corner_lookup[edge.label]
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
                    x, y, z = lamination.triangulation.corner_lookup[edge.label]
                    if parallels[edge.index] > 0:
                        geometric[edge.index] = -parallels[x.index]
                    else:
                        geometric[edge.index] = max(sides[x], 0) + max(sides[y], 0) + max(-sides[z], 0)
                        
                        # Sanity check:
                        x2, y2, z2 = lamination.triangulation.corner_lookup[~edge.label]
                        assert(geometric[edge.index] == max(sides[x2], 0) + max(sides[y2], 0) + max(-sides[z2], 0))
            
            # We have to rebuild the ~e edge separately since it now pairs with ~b.
            x, y, z = lamination.triangulation.corner_lookup[~e.label]
            if parallels[e.index] + parallels[b.index] > 0:
                geometric[e.index] = -parallels[e.index] - parallels[b.index]
            else:
                geometric[e.index] = max(sides[x], 0) + max(sides[y], 0) + max(-sides[z], 0)
                
                # Sanity check:
                x2, y2, z2 = lamination.triangulation.corner_lookup[~b.label]
                assert(geometric[e.index] == max(sides[x2], 0) + max(sides[y2], 0) + max(-sides[z2], 0))
            
            # And finally the b edge, which is now paired with e.
            geometric[b.index] = around_v - out_v  # Same trick as below.
        else:  # self.curve is non-isolating.
            # Get some edges.
            a = self.curve.parallel()
            _, b, e = self.source_triangulation.corner_lookup[a.label]
            _, c, d = self.source_triangulation.corner_lookup[~e.label]
            
            u = self.source_triangulation.vertex_lookup[a.label]  # = self.triangulation.vertex_lookup[~a.label].
            u_edges = curver.kernel.utilities.cyclic_slice(u, a, ~a)
            around_u = min(max(lamination.side_weight(edge), 0) for edge in u_edges)  # The number of components that go through a, around u and then back through ~a.
            out_u = sum(max(-lamination.side_weight(edge), 0) for edge in u_edges) + sum(max(-lamination(edge), 0) for edge in u_edges[1:])  # The number of arcs that come out of u (from a around to ~a).
            # Since around_u > 0 ==> out_u == 0 & out_u > 0 ==> around_u == 0, this is equivalent to around_u if around_u > 0 else -out_u
            geometric[b.index] = around_u - out_u
            
            v = self.source_triangulation.vertex_lookup[c.label]  # = self.triangulation.vertex_lookup[~c.label].
            v_edges = curver.kernel.utilities.cyclic_slice(v, c, ~c)
            around_v = min(max(lamination.side_weight(edge), 0) for edge in v_edges)
            out_v = sum(max(-lamination.side_weight(edge), 0) for edge in v_edges) + sum(max(-lamination(edge), 0) for edge in v_edges[1:])
            geometric[e.index] = around_v - out_v  # Same trick.
        
        return self.target_triangulation(geometric)  # Have to promote.
    
    def apply_homology(self, homology_class):
        return NotImplemented  # I don't think we ever need this.
    
    def inverse(self):
        return Lift(self.target_triangulation, self.source_triangulation, self.curve, self.matrix)

class Lift(Move):
    ''' This represents the inverse of crushing along a curve. '''
    def __init__(self, source_triangulation, target_triangulation, curve, matrix):
        super(Lift, self).__init__(source_triangulation, target_triangulation)
        
        assert(isinstance(curve, curver.kernel.Curve))
        assert(curve.triangulation == self.target_triangulation)
        
        self.curve = curve
        self.matrix = matrix
    
    def __str__(self):
        return 'Lift ' + str(self.curve)
    def __reduce__(self):
        return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve, self.matrix))
    def __eq__(self, other):
        if isinstance(other, Lift):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation and self.curve == other.curve and self.matrix == other.matrix
        else:
            return NotImplemented
    def __ne__(self, other):
        return not (self == other)
    
    def apply_lamination(self, lamination):
        # Really should check that the dual weights around a vertex are all non-negative.
        geometric = curver.kernel.matrix_vector(self.matrix, lamination.geometric)
        return lamination.__class__(self.target_triangulation, geometric)  # Avoid promote since the lift has to be the same type as the given lamination.
    
    def apply_homology(self, homology_class):
        return NotImplemented  # I don't think we ever need this.
    
    def inverse(self):
        return Crush(self.target_triangulation, self.source_triangulation, self.curve, self.matrix)

