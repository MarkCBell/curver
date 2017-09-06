
''' A module for representing (multi)arcs on triangulations.

Provides: MultiArc, Arc. '''

import curver
from curver.kernel.lamination import Shortenable  # Special import needed for subclassing.

class MultiArc(Shortenable):
	''' A Lamination in which every component is an Arc. '''
	def is_multicurve(self):
		return False
	def is_multiarc(self):
		return True
	
	def shorten_strategy(self, edge):
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
		
		geometric = [0 if weight < 0 else 2 for weight in short]
		# Tighten by retracting from any triangle where geometric meets only one side.
		to_fix = [triangle for triangle in short.triangulation if sum(geometric[index] for index in triangle.indices) == 2]  # A stack of work.
		while to_fix:
			triangle = to_fix.pop()
			if sum(geometric[index] for index in triangle.indices) == 2:
				for edge in triangle:
					geometric[edge.index] = 0
					to_fix.append(short.triangulation.triangle_lookup[~edge.label])
		
		boundary = short.triangulation.lamination(geometric)
		
		return conjugator.inverse()(boundary)
	
	def is_polygonalisation(self):
		''' Return if this MultiArc is a polygonalisation, that is, if it cuts the surface into polygons. '''
		short, _ = self.shorten()
		
		avoid = set(index for index in short.triangulation.indices if short(index) < 0)  # All of the edges used.
		dual_tree = short.triangulation.dual_tree(avoid=avoid)
		
		return all(dual_tree[index] or index in avoid for index in short.triangulation.indices)
	
	def explore_ball(self, radius):
		''' Extend this MultiArc to a triangulation and return all triangulations within the ball of the given radius of that one.
		
		Runs in exp(radius) time.
		Note that this is only well-defined if this multiarc is filling. '''
		
		assert(self.is_filling())
		
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
		''' Return an edge that this arc is parallel to. '''
		assert(self.is_short())
		
		return min([edge for edge in self.triangulation.edges if self(edge) < 0], key=lambda e: e.label)
	
	def encode_halftwist(self, power=1):
		''' Return an Encoding of a right half twist about a regular neighbourhood of this arc, raised to the given power.
		
		Assumes (and checks) that this arc connects between distinct vertices. '''
		
		short, conjugator = self.shorten()
		
		edge = short.parallel()
		
		# Check where it connects.
		if short.triangulation.vertex_lookup[edge.label] == short.triangulation.vertex_lookup[~edge.label]:
			raise curver.AssumptionError('Arc connects a vertex to itself.')
		
		return conjugator.inverse() * curver.kernel.HalfTwist(short, power).encode() * conjugator
	
	def intersection(self, lamination):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(lamination, curver.kernel.Lamination))
		assert(lamination.triangulation == self.triangulation)
		
		short, conjugator = self.shorten()
		short_lamination = conjugator(lamination)
		
		arc = short.parallel()
		
		return max(short_lamination(arc), 0)

