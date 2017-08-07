
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
		geometric = [lamination(self.inverse_index_map[index]) for index in self.source_triangulation.indices]
		return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
	
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

