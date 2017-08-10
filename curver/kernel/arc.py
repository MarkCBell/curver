
import curver
from curver.kernel.lamination import Shortenable  # Special import needed for subclassing.

class MultiArc(Shortenable):
	''' A Lamination in which every component is an Arc. '''
	def is_multicurve(self):
		return False
	def is_multiarc(self):
		return True
	
	def is_short(self):
		return self.weight() == 0
	
	def shorten_score(self, edge):
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		# Low score == bad.
		if not self.triangulation.is_flippable(edge):
			return 0
		if self.dual_weight(edge) < 0:
			return 1
		
		return 0
	
	def boundary(self):
		''' Return the multicurve which is the boundary of a regular neighbourhood of this multiarc. '''
		
		short, conjugator = self.shorten()
		# short is a subset of the edges of the triangulation it is defined on.
		# So its geometric vector is non-positive.
		
		# TODO: 2) Make linear instead of quadratic by using a queue / stack.
		geometric = [0 if weight < 0 else 2 for weight in short]
		# Tighten.
		changed = True
		while changed:
			changed = False
			for triangle in short.triangulation:
				if sum(geometric[index] for index in triangle.indices) == 2:
					for index in triangle.indices:
						geometric[index] = 0
					changed = True
		
		boundary = short.triangulation.lamination(geometric)
		
		return conjugator.inverse()(boundary)
	
	def is_filling(self):
		short, _ = self.shorten()
		
		avoid = set(index for index in short.triangulation.indices if short(index) < 0)  # All of the edges used.
		dual_tree = short.triangulation.dual_tree(avoid=avoid)
		
		return all(dual_tree[index] or index in avoid for index in short.triangulation.indices)
	
	def explore_ball(self, radius):
		''' Extend this MultiArc to a triangulation and return all triangulations within the ball of the given radius of that one.
		
		Note that this is only well-defined if self is filling. '''
		
		short, conjugator = self.shorten()
		
		triangulations = set()
		for encoding in short.triangulation.all_encodings(radius):
			T = encoding.target_triangulation.as_lamination()
			triangulations.add(conjugator.inverse()(encoding.inverse()(T)))
		
		return triangulations

class Arc(MultiArc):
	''' A MultiArc with a single component. '''
	def mcomponents(self):
		return [(self, 1)]
	
	def parallel(self):
		''' Return an edge that this curve is parallel to. '''
		assert(self.is_short())
		
		return min([edge for edge in self.triangulation.edges if self(edge) < 0], key=lambda e: e.label)
	
	def encode_halftwist(self, k=1):
		''' Return an Encoding of a left half twist about a regular neighbourhood of this arc, raised to the power k.
		
		Assumes (and checks) that this arc connects between distinct vertices. '''
		# Will Worden checked that this works for genus <= 20.
		
		# Some easy cases:
		if k == 0: return self.triangulation.id_encoding()
		if k < 0: return self.encode_halftwist(-k).inverse()
		
		# We need to get to a really good configuration, one where self
		# is an edge of the triangulation whose valence(initial vertex) == 1.
		#
		# We achieve this in two steps. First conjugate to make self an edge of some triangulation.
		short, conjugator = self.shorten()
		arc = short.parallel()
		
		# Check where it connects.
		if short.triangulation.vertex_lookup[arc.label] == short.triangulation.vertex_lookup[~arc.label]:
			raise curver.AssumptionError('Arc connects a vertex to itself.')
		
		# Now keep moving edges away from this edge's initial vertex to get to a really good triangulation.
		while len(short.triangulation.vertex_lookup[arc.label]) > 1:  # valence(initial vertex) > 1.
			flip = short.triangulation.encode_flip(short.triangulation.corner_lookup[arc.label][2])
			conjugator = flip * conjugator
			short = flip(short)
		
		# We can now perform the half twist. To do this we move all the edges back across to the other vertex.
		# Again, we keep moving edges away from this edge's terminal vertex.
		half_twist = short.triangulation.id_encoding()  # valence(terminal vertex) > 1.
		while len(short.triangulation.vertex_lookup[~arc.label]) > 1:
			flip = short.triangulation.encode_flip(short.triangulation.corner_lookup[~arc.label][2])
			half_twist = flip * half_twist
			short = flip(short)
		
		# No close up to complete the half twist. This means finding the correct isometry back to the
		# really good triangulation. We want the isometry to be the identity on all other components
		# and on the one containing this arc to invert this arc.
		half_twist = short.triangulation.find_isometry(half_twist.source_triangulation, {arc.label: ~arc.label}).encode() * half_twist
		
		# We handle large powers by replacing (T^1/2_self)^2 with T_boundary which includes acceleration.
		return self.boundary().encode_twist(k // 2) * conjugator.inverse() * half_twist * conjugator
	
	def intersection(self, other):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(other, Lamination))
		assert(other.triangulation == self.triangulation)
		
		# short = [0, 0, ..., 0, -1, 0, ..., 0]
		short, conjugator = self.shorten()
		short_other = conjugator(other)
		
		return sum(b for a, b in zip(short, short_other) if a == -1 if b >= 0)

