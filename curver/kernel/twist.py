
''' A module for representing more advanced ways of changing triangulations. '''

import curver
from curver.kernel.moves import FlipGraphMove  # Special import needed for subclassing.

class Twist(FlipGraphMove):
    ''' This represents the effect of twisting a short curve.
    
    This format allows us to efficiently perform powers of twists. '''
    def __init__(self, curve, power):
        super(Twist, self).__init__(curve.triangulation, curve.triangulation)
        
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
    
    def __str__(self):
        return 'Twist^%d_%s ' % (self.power, self.curve)
    def __reduce__(self):
        return (self.__class__, (self.curve, self.power))
    def package(self):
        return (self.curve.parallel().label, self.power)
    def __eq__(self, other):
        if isinstance(other, Twist):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation and self.curve == other.curve and self.power == other.power
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    
    def apply_lamination(self, lamination):
        intersection = self.curve.intersection(lamination)
        if intersection == 0:  # Disjoint twists have no effect.
            return lamination
        
        # Naive way would be to do:
        # for i in range(self.power):
        #     lamination = self.encoding(lamination)
        # return lamination
        # But we can be cleverer and perform this calculation in O(log(self.power)) instead.
        
        power = self.power
        while power:  # This loop will run at most four times.
            slope = self.curve.slope(lamination)
            if power > 0:  # Right twist (increases slope).
                if 1 < slope:  # pylint: disable=misplaced-comparison-constant
                    geometric = [w + power * intersection * c for w, c in zip(lamination, self.curve)]
                    power_applied = power
                elif -1 <= slope <= 1:  # Dangerous region.
                    geometric = self.encoding(lamination).geometric
                    power_applied = 1
                else:  # if slope < -1:
                    steps = min(power, -(slope.numerator // slope.denominator) - 1)  # ceil(-slope) - 1.
                    assert steps >= 0  # Sanity.
                    geometric = [w - steps * intersection * c for w, c in zip(lamination, self.curve)]
                    power_applied = steps
            else:  # power < 0:  # Left twist (decreases slope).
                if slope < -1:
                    geometric = [w + -power * intersection * c for w, c in zip(lamination, self.curve)]
                    power_applied = power
                elif -1 <= slope <= 1:  # Dangerous region.
                    geometric = self.encoding.inverse()(lamination).geometric
                    power_applied = -1
                else:  # 1 < slope:
                    steps = min(-power, -(-slope.numerator // slope.denominator) - 1)  # ceil(slope) - 1.
                    assert steps >= 0  # Sanity.
                    geometric = [w - steps * intersection * c for w, c in zip(lamination, self.curve)]
                    power_applied = -steps
            new_lamination = lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
            # assert new_lamination == (self.encoding**power_applied)(lamination)
            # assert -1 <= slope <= 1 or self.curve.slope(new_lamination) == slope + power_applied
            lamination = new_lamination
            power = power - power_applied
        
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

class HalfTwist(FlipGraphMove):
    ''' This represents the effect of half-twisting a short arc.
    
    This format allows us to efficiently perform powers of twists. '''
    def __init__(self, arc, power):
        super(HalfTwist, self).__init__(arc.triangulation, arc.triangulation)
        
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
        
        self.encoding = conjugator.inverse() * half_twist * conjugator
        
        # We handle large powers by replacing (T^1/2_self)^2 with T_boundary, which includes acceleration.
        # We handle small powers separately to increase performance.
        if abs(self.power) <= 1:
            self.encoding_power = self.encoding**self.power
        elif self.power % 2 == 0:
            self.encoding_power = self.arc.boundary().encode_twist(self.power // 2)
        else:  # self.power % 2 == 1:  # Division rounds down so, regardless of power, we need an extra right half-twist.
            self.encoding_power = self.arc.boundary().encode_twist(self.power // 2) * self.encoding
    
    def __str__(self):
        return 'HalfTwist^%d_%s ' % (self.power, self.arc)
    def __reduce__(self):
        return (self.__class__, (self.arc, self.power))
    def package(self):
        return (self.arc.parallel().label, self.power)
    def __eq__(self, other):
        if isinstance(other, HalfTwist):
            return self.source_triangulation == other.source_triangulation and self.target_triangulation == other.target_triangulation and self.arc == other.arc and self.power == other.power
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    
    def apply_lamination(self, lamination):
        return self.encoding_power(lamination)
    
    def apply_homology(self, homology_class):
        return self.encoding_power(homology_class)
    
    def inverse(self):
        return HalfTwist(self.arc, -self.power)
    
    def flip_mapping(self):
        return self.encoding**self.power
