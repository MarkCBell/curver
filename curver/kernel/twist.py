
''' A module for representing more advanced ways of changing triangulations. '''

import numpy as np

import curver
from curver.kernel.moves import FlipGraphMove  # Special import needed for subclassing.
from curver.kernel.decorators import ensure

class Twist(FlipGraphMove):
    ''' This represents the effect of twisting a short curve.
    
    This format allows us to efficiently perform powers of twists. '''
    def __init__(self, curve, power):
        super().__init__(curve.triangulation, curve.triangulation)
        
        assert isinstance(curve, curver.kernel.Curve)
        assert curve.is_short() and not curve.is_peripheral()
        assert power != 0
        
        self.curve = curve
        self.power = power
        
        a = self.curve.parallel()
        # Theorem: The right number of flips to do is:
        #  - weight - self.curve.dual_weight(parallel) in the non-isolating case
        #  - 3*num_tripods is in the isolating case.
        # Proof: TODO.
        num_flips = self.curve.weight() - self.curve.dual_weight(a)
        
        twist = self.curve.triangulation.id_encoding()
        for _ in range(num_flips):
            twist = twist.target_triangulation.encode_flip(twist.target_triangulation.corner_lookup[a][2]) * twist
        twist = twist.target_triangulation.find_isometry(twist.source_triangulation, {a.label: a.label}).encode() * twist
        
        self.encoding = twist
        self.power_sign = +1 if self.power > 0 else -1
        self.signed_encoding = self.encoding if self.power > 0 else self.encoding.inverse()
    
    def __str__(self):
        return f'Twist^{self.power}_{self.curve}'
    def package(self):
        return (self.curve.parallel().label, self.power)
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.curve == other.curve and self.power == other.power
    
    def apply_lamination(self, lamination):
        # Take care of an easy case for speed.
        if abs(self.power) == 1:
            return self.signed_encoding(lamination)
        
        intersection = self.curve.intersection(lamination)
        if intersection == 0:  # Disjoint twists have no effect.
            return lamination
        
        # Naive way would be to do:
        # return self.encoding(lamination, power=self.power)
        # which is roughly equivalent to:
        # for i in range(self.power):
        #     lamination = self.encoding(lamination)
        # return lamination
        # But we can be cleverer and perform this calculation in O(log(self.power)) instead.
        
        power = self.power
        slope = self.curve.slope(lamination)
        slope_sign = +1 if slope > 0 else -1
        steps = min(abs(power), int(abs(slope)))  # This is how far we can definitely move without reaching the dangerous region.
        
        if power * slope_sign < 0:  # We are heading towards the dangerous region.
            lamination = lamination.__class__(self.target_triangulation, [w - steps * intersection * c for w, c in zip(lamination, self.curve)])  # Avoids promote.
            power = power + slope_sign * steps
        
        # We have to go slowly through the dangerous region, but we cross it in at most three twists.
        for _ in range(3):
            if power:
                lamination = self.signed_encoding(lamination)
                power = power - self.power_sign
        
        if power:  # Since we are now moving away from the dangerous region, we can accelerate.
            lamination = lamination.__class__(self.target_triangulation, [w + abs(power) * intersection * c for w, c in zip(lamination, self.curve)])  # Avoids promote.
        
        return lamination
    
    def apply_homology(self, homology_class):
        a = self.curve.parallel()
        
        v = self.source_triangulation.vertex_lookup[a]  # = self.source_triangulation.vertex_lookup[~a].
        v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
        
        algebraic = list(homology_class)
        algebraic[a.index] += a.sign() * self.power * sum(homology_class(edge) for edge in v_edges[1:])
        
        return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
    
    def inverse(self):
        return Twist(self.curve, -self.power)
    
    def flip_mapping(self):
        return self.encoding**self.power
    
    @ensure(lambda data: data.result(data.multicurve.geometric) == data.self(data.multicurve).geometric)
    def pl_action(self, multicurve):
        # Take care of an easy case for speed.
        if abs(self.power) == 1:
            return self.signed_encoding.pl_action(multicurve)
        
        # Helper functions for building matrices.
        def V(edge):
            return np.array([1 if i == edge.index else 0 for i in range(self.zeta)], dtype=object)
        
        def C2(edge):
            corner = self.source_triangulation.corner_lookup[edge]
            return V(corner[0]) + V(corner[2]) - V(corner[1])
        
        power = self.power
        
        # Slope calculation:
        # Get some edges.
        a = self.curve.parallel()
        v = self.curve.triangulation.vertex_lookup[a]
        
        v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)
        around = curver.kernel.utilities.maximin([0], (multicurve.left_weight(edgy) for edgy in v_edges))
        around_edge = next(edge for edge in v_edges if multicurve.left_weight(edge) == around)  # The edge that realises around.
        
        twisting = curver.kernel.utilities.maximin([0], (multicurve.left_weight(edgy) - around for edgy in v_edges[1:-1]))
        twisting_edge = next(edge for edge in v_edges[1:-1] if multicurve.left_weight(edge) - around == twisting)  # The edge that realises twisting.
        
        slope_sign = -1 if multicurve.left_weight(a) - around > 0 else +1
        intersection = multicurve(a) - 2 * around  # = self.curve.intersection(multicurve)
        
        # Condition matrices which restricts to multicurves with the same around_edge, twisting_edge and slope sign respectively.
        around_condition = np.array([C2(edge) - C2(around_edge) for edge in v_edges])
        twisting_condition = np.array([C2(edge) - C2(twisting_edge) for edge in v_edges[1:-1]])
        slope_sign_condition = np.array([slope_sign * (C2(around_edge) - C2(a))])
        
        numerator, denominator = twisting, intersection
        if denominator == 0:  # Disjoint twists have no effect.
            return curver.kernel.PartialLinearFunction(
                np.identity(self.zeta, dtype=object),
                np.concatenate([around_condition, np.array([C2(around_edge) - V(a), V(a) - C2(around_edge)])])
                )
        
        floor_abs_slope = numerator // denominator
        steps = min(abs(power), floor_abs_slope)
        
        # A condition matrix which restricts to multicurves with the same floor_abs_slope.
        # Note we use an extra factor of two to avoid fractions.
        floor_abs_slope_condition = np.array([
            (C2(twisting_edge) - C2(around_edge)) - 2 * floor_abs_slope * (V(a) - C2(around_edge)),  # 2 * twisting - 2*floor_abs_slope * intersection.
            2 * (floor_abs_slope+1) * (V(a) - C2(around_edge)) - (C2(twisting_edge) - C2(around_edge))  # 2 * (floor_abs_slope+1) * intersection - 2 * twisting.
            ])
        
        # Start with all of these constraints in the PL function.
        F = curver.kernel.PartialLinearFunction(
            np.identity(self.zeta, dtype=object),
            np.concatenate([around_condition, twisting_condition, slope_sign_condition, floor_abs_slope_condition])
            )
        
        if power * slope_sign < 0:
            F = curver.kernel.PartialLinearFunction(
                np.array([V(edge) - steps * (V(a) - C2(around_edge)) * self.curve(edge) for edge in self.source_triangulation.positive_edges]),
                np.array([[0] * self.zeta], dtype=object)
                ) * F
            multicurve = multicurve.__class__(self.target_triangulation, [w - steps * intersection * c for w, c in zip(multicurve, self.curve)])  # Avoids promote.
            power = power + slope_sign * steps
        
        for _ in range(3):
            if power:
                F = self.signed_encoding.pl_action(multicurve) * F
                multicurve = self.signed_encoding(multicurve)
                power = power - self.power_sign
        
        if power:
            # We now have to recalculate around.
            around = curver.kernel.utilities.maximin([0], (multicurve.left_weight(edgy) for edgy in v_edges))
            around_edge = next(edge for edge in v_edges if multicurve.left_weight(edge) == around)  # The edge that realises around.
            around_condition = np.array([C2(edge) - C2(around_edge) for edge in v_edges])
            
            F = curver.kernel.PartialLinearFunction(
                np.array([V(edge) + abs(power) * (V(a) - C2(around_edge)) * self.curve(edge) for edge in self.source_triangulation.positive_edges]),
                around_condition,
                ) * F
            multicurve = multicurve.__class__(self.target_triangulation, [w + abs(power) * intersection * c for w, c in zip(multicurve, self.curve)])  # Avoids promote.
        
        return F

class HalfTwist(FlipGraphMove):
    ''' This represents the effect of half-twisting a short arc.
    
    This format allows us to efficiently perform powers of twists. '''
    def __init__(self, arc, power):
        super().__init__(arc.triangulation, arc.triangulation)
        
        assert isinstance(arc, curver.kernel.Arc)
        assert arc.is_short()
        assert power != 0
        
        self.arc = arc
        self.power = power
        
        conjugator = arc.triangulation.id_encoding()
        # We need to get to a really good configuration, one where self.arc is not just short
        # but where valence(self.arc.initial_vertex) == 1.
        
        edge = self.arc.parallel()
        # Reverse the orientation if the valence of the other end is less.
        # This reduces the number of flips needed to reach a really good configuration.
        if len(arc.triangulation.vertex_lookup[edge]) > len(arc.triangulation.vertex_lookup[~edge]):
            edge = ~edge
        
        # Since self.arc is short it is an edge of the triangulation so we just keep moving
        # edges away from this edge's initial vertex to get to a really good triangulation.
        while len(conjugator.target_triangulation.vertex_lookup[edge]) > 1:  # valence(initial vertex) > 1.
            flip = conjugator.target_triangulation.encode_flip(conjugator.target_triangulation.corner_lookup[edge][2])
            conjugator = flip * conjugator
        
        # We can now perform the half twist. To do this we move all the edges back across to the other vertex.
        # Again, we keep moving edges away from this edge's terminal vertex.
        # TODO: 4) Prove this always works.
        # NOTE: William Worden checked that this works for genus <= 20.
        half_twist = conjugator.target_triangulation.id_encoding()  # valence(terminal vertex) > 1.
        while len(half_twist.target_triangulation.vertex_lookup[~edge]) > 1:
            flip = half_twist.target_triangulation.encode_flip(half_twist.target_triangulation.corner_lookup[~edge][2])
            half_twist = flip * half_twist
        
        # No close up to complete the half twist. Use the isometry that inverts this edge.
        half_twist = half_twist.target_triangulation.find_isometry(half_twist.source_triangulation, {edge.label: ~edge.label}).encode() * half_twist
        
        self.encoding = half_twist.conjugate_by(conjugator)
        
        # We handle large powers by replacing (T^1/2_self)^2 with T_boundary, which includes acceleration.
        # We handle small powers separately to increase performance.
        if abs(self.power) <= 1:
            self.encoding_power = self.encoding**self.power
        elif self.power % 2 == 0:
            self.encoding_power = self.arc.boundary().encode_twist(self.power // 2)
        else:  # self.power % 2 == 1:  # Division rounds down so, regardless of power, we need an extra right half-twist.
            self.encoding_power = self.arc.boundary().encode_twist(self.power // 2) * self.encoding
    
    def __str__(self):
        return f'HalfTwist^{self.power}_{self.arc}'
    def package(self):
        return (self.arc.parallel().label, self.power)
    def __eq__(self, other):
        eq = super().__eq__(other)
        if eq in [NotImplemented, False]:
            return eq
        
        return self.arc == other.arc and self.power == other.power
    
    def apply_lamination(self, lamination):
        return self.encoding_power(lamination)
    
    def apply_homology(self, homology_class):
        return self.encoding_power(homology_class)
    
    def inverse(self):
        return HalfTwist(self.arc, -self.power)
    
    def flip_mapping(self):
        return self.encoding**self.power
    
    def pl_action(self, multicurve):
        return self.encoding_power.pl_action(multicurve)
