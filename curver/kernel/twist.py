
''' A module for representing more advanced ways of changing triangulations.

Provides: Twist, HalfTwist. '''

import curver
from curver.kernel.moves import Move  # Special import needed for subclassing.

class Twist(Move):
	''' This represents the effect of twisting a short curve.
	
	This format allows us to efficiently perform powers of twists. '''
	def __init__(self, curve, power):
		super(Twist, self).__init__(curve.triangulation, curve.triangulation)
		
		assert(isinstance(curve, curver.kernel.Curve))
		assert(curve.is_short())
		
		self.curve = curve
		self.power = power
		
		a = self.curve.parallel()
		if self.curve.is_isolating():
			# Theorem: 3*num_tripods is the right number of flips to do.
			# Proof: TODO.
			num_tripods = len([triangle for triangle in self.curve.triangulation if sum(self.curve(edge) for edge in triangle) == 6])
			num_flips = 3*num_tripods
		else:  # curve is non-isolating.
			num_flips = 1
		
		twist = self.curve.triangulation.id_encoding()
		for _ in range(num_flips):
			twist = twist.target_triangulation.encode_flip(twist.target_triangulation.corner_lookup[a.label][2]) * twist
		twist = twist.target_triangulation.find_isometry(twist.source_triangulation, {a.label: a.label}).encode() * twist
		
		self.encoding = twist
	
	def __str__(self):
		return 'Twist^%d %s ' % (self.power, self.curve)
	def __reduce__(self):
		return (self.__class__, (self.curve, self.power))
	
	def apply_lamination(self, lamination):
		intersection = self.curve.intersection(lamination)
		if intersection == 0:  # Disjoint twists have no effect.
			return lamination
		
		power = self.power
		while power:
			slope = self.curve.slope(lamination)
			if power > 0:  # Left twist. (Reduces slope)
				if slope < -2:
					geometric = [w + power * intersection * c for w, c in zip(lamination, self.curve)]
					power = 0
				elif -2 <= slope <= 2:  # Dangerous region.
					geometric = self.encoding(lamination).geometric
					power = power - 1
				else:  # if 2 < slope:
					steps = min(power, int(slope-1))
					geometric = [w - steps * intersection * c for w, c in zip(lamination, self.curve)]
					power = power - steps
			else:  # power < 0:  # Right twist.
				if 2 < slope:
					geometric = [w + -power * intersection * c for w, c in zip(lamination, self.curve)]
					power = 0
				elif -2 <= slope <= 2:  # Dangerous region.
					geometric = self.encoding.inverse()(lamination).geometric
					power = power + 1
				else:  # slope < -2:
					steps = min(-power, int(-1-slope))
					geometric = [w - steps * intersection * c for w, c in zip(lamination, self.curve)]
					power = power + steps
			lamination = lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
		
		return lamination
	
	def apply_homology(self, homology_class):
		# I don't think we even need this case.
		if self.curve.is_isolating(): # Isolating ==> separating so no effect on homology.
			return homology_class
		
		a = self.curve.parallel()
		_, b, e = self.source_triangulation.corner_lookup[a.label]
		
		algebraic = list(homology_class)
		algebraic[a.index] -= a.sign() * self.power * (homology_class(e) - homology_class(b))
		
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		return Twist(self.curve, -self.power)

class HalfTwist(Move):
	''' This represents the effect of half-twisting a short arc.
	
	This format allows us to efficiently perform powers of twists. '''
	def __init__(self, arc, power):
		super(HalfTwist, self).__init__(arc.triangulation, arc.triangulation)
		
		assert(isinstance(arc, curver.kernel.Arc))
		assert(arc.is_short())
		
		self.arc = arc
		self.power = power
		
		edge = self.arc.parallel()
		
		conjugator = arc.triangulation.id_encoding()
		# We need to get to a really good configuration, one where self.curve is not just short
		# but where valence(self.curve.initial_vertex) == 1.
		
		# Since self.curve is short it is an edge of the triangulation so we just keep moving
		# edges away from this edge's initial vertex to get to a really good triangulation.
		while len(conjugator.target_triangulation.vertex_lookup[edge.label]) > 1:  # valence(initial vertex) > 1.
			flip = conjugator.target_triangulation.encode_flip(conjugator.target_triangulation.corner_lookup[edge.label][2])
			conjugator = flip * conjugator
		
		# We can now perform the half twist. To do this we move all the edges back across to the other vertex.
		# Again, we keep moving edges away from this edge's terminal vertex.
		# TODO: 4) Prove this always works.
		# NOTE: Will Worden checked that this works for genus <= 20.
		half_twist = conjugator.target_triangulation.id_encoding()  # valence(terminal vertex) > 1.
		while len(half_twist.target_triangulation.vertex_lookup[~edge.label]) > 1:
			flip = half_twist.target_triangulation.encode_flip(half_twist.target_triangulation.corner_lookup[~edge.label][2])
			half_twist = flip * half_twist
		
		# No close up to complete the half twist. Use the isometry that inverts this edge.
		half_twist = half_twist.target_triangulation.find_isometry(half_twist.source_triangulation, {edge.label: ~edge.label}).encode() * half_twist
		
		self.encoding = conjugator.inverse() * half_twist * conjugator
		
		# We handle large powers by replacing (T^1/2_self)^2 with T_boundary, which includes acceleration.
		if self.power % 2 == 0:
			self.encoding_power = self.arc.boundary().encode_twist(self.power // 2)
		else:  # self.power % 2 == 1:  # Division rounds down so, regardless of power, we need an extra left half-twist.
			self.encoding_power = self.arc.boundary().encode_twist(self.power // 2) * self.encoding
	
	def __str__(self):
		return 'HalfTwist^%d %s ' % (self.power, self.arc)
	def __reduce__(self):
		return (self.__class__, (self.arc, self.power))
	
	def apply_lamination(self, lamination):
		return self.encoding_power(lamination)
	
	def apply_homology(self, homology_class):
		return self.encoding_power(homology_class)
	
	def inverse(self):
		return HalfTwist(self.arc, -self.power)
