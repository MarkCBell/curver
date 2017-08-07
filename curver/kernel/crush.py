
import curver
from curver.kernel.moves import Move  # Special import needed for subclassing.

class Crush(Move):
	def __init__(self, source_triangulation, target_triangulation, curve, matrix):
		super(Crush, self).__init__(source_triangulation, target_triangulation)
		
		assert(isinstance(curve, curver.kernel.Curve))
		assert(curve.is_short())
		assert(curve.triangulation == self.source_triangulation)
		
		self.curve = curve
		self.matrix = matrix
	
	def __str__(self):
		return 'Crush ' + str(self.curve)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve, self.matrix))
	
	def apply_lamination(self, lamination):
		
		if self.curve.weight() == 2:
			e1, e2 = [index for index in self.source_triangulation.indices if self.curve(index) > 0]
			# We might need to swap these edge indices so we have a good frame of reference.
			if self.source_triangulation.corner_lookup[e1].indices[2] != e2: e1, e2 = e2, e1
			
			a, b, c, d, e = self.source_triangulation.square(e1)
			u, v, w, x, y, z = lamination.dual_square(e)
			
			geometric = list(lamination)
			if lamination(b) < 0 or lamination(e) < 0:
				new_b = new_e = lamination(b) + lamination(e)
			else:
				new_b = min(v, w, y)
				new_e = min(u, x, z)
			geometric[b.index] = new_b
			geometric[e.index] = new_e
		else:
			# TODO: 4) Implement LP to find intersection for general configuration.
			raise curver.AssumptionError('Currently can only crush along non-isolating curves.')
		
		return curver.kernel.Lamination(self.target_triangulation, geometric).promote()  # Have to promote.
	
	def apply_homology(self, homology_class):
		return NotImplemented  # I don't think we ever need this.
	
	def inverse(self):
		return Lift(self.target_triangulation, self.source_triangulation, self.curve, self.matrix)

class Lift(Move):
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
	
	def apply_lamination(self, lamination):
		# Really should check that the dual weights around a vertex are all non-negative.
		geometric = [sum(x * y for x, y in zip(row, lamination)) for row in self.matrix]  # Dot product.
		return curver.kernel.Lamination(self.target_triangulation, geometric).promote()  # Have to promote.
	
	def apply_homology(self, homology_class):
		return NotImplemented  # I don't think we ever need this.
	
	def inverse(self):
		return Crush(self.target_triangulation, self.source_triangulation, self.curve, self.matrix)

