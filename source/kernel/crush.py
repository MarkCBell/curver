
import curver
from curver.kernel.moves import Move  # Special import needed for subclassing.

class Crush(Move):
	def __init__(self, source_triangulation, target_triangulation, curve):
		super(Crush, self).__init__(source_triangulation, target_triangulaion)
		
		assert(isinstance(curve, curver.Curve))
		assert(curve.is_short())
		assert(curve.triangulation == self.source_triangulation)
		
		self.curve = curve
	
	def __str__(self):
		return 'Crush ' + str(self.curve)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return NotImplemented
	
	def apply_lamination(self, lamination):
		
		if self.curve.weight() == 2:
			triangulation = short.triangulation
			e1, e2 = [index for index in triangulation.indices if short(index) > 0]
			# We might need to swap these edge indices so we have a good frame of reference.
			if triangulation.corner_lookup[e1].indices[2] != e2: e1, e2 = e2, e1
			
			a, b, c, d, e = triangulation.square(e1)
			u, v, w, x, y, z = lamination.dual_square(e)
			
			geometric = list(lamination)
			if lamination(b) < 0 or lamination(e) < 0:
				new_b = new_e = lamination(b) + lamination(e)
			else:
				new_b = min(v, w, y)
				new_e = min(u, x, z))
			geometric[b.index] = new_b
			geometric[e.index] = new_e
		else:
			# TODO: 4) Implement LP to find intersection for general configuration.
			raise curver.AssumptionError('Currently can only crush along non-isolating curves.')
		
		return curver.kenel.Lamination(self.target_triangulation, geometric).promote()  # Have to promote.
	
	def apply_homology(self, homology_class):
		return NotImplemented  # I don't think we ever need this.
		algebraic = []
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		return Lift(self.target_triangulation, self.source_triangulation, self.curve)

class Lift(Move):
	def __init__(self, source_triangulation, target_triangulation, curve):
		super(Crush, self).__init__(source_triangulation, target_triangulaion)
		
		assert(isinstance(curve, curver.Curve))
		assert(curve.triangulation == self.target_triangulation)
		
		self.curve = curve
	
	
	def __str__(self):
		return 'Crush ' + str(self.curve)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return NotImplemented
	
	def apply_lamination(self, lamination):
		geometric = [lamination(self.inverse_index_map[index]) for index in self.source_triangulation.indices]
		return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
	
	def apply_homology(self, homology_class):
		return NotImplemented  # I don't think we ever need this.
		algebraic = []
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		return Crush(self.target_triangulation, self.source_triangulation, self.curve)

