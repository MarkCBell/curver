
''' A module for representing (multi)curves on triangulations.

Provides: MultiCurve, Curve. '''

import curver
from curver.kernel.lamination import Lamination, Shortenable  # Special import needed for subclassing.

class MultiCurve(Lamination):
	''' A Lamination in which every component is a Curve. '''
	def is_multicurve(self):
		return True
	def is_multiarc(self):
		return False
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn (multi)twist about the components of this multicurve, raised to the power k. '''
		
		h = self.triangulation.id_encoding()
		for curve, multiplicity in self.mcomponents():
			h = curve.encode_twist(k * multiplicity) * h
		
		return h
	
	def boundary_union(self, other):
		''' Return \\partial N(self \\cup other). '''
		assert(isinstance(other, Lamination))
		
		crush = self.crush()
		lift = crush.inverse()
		other_prime = crush(other)
		m_prime = other_prime.boundary()
		return lift(m_prime)  # = m.
	
	def fills_with(self, other):
		''' Return whether self \\cup other fills. '''
		assert(isinstance(other, Lamination))
		
		return self.boundary_union(other).is_empty()
	
	def crush(self):
		''' Return the crush map associated to this MultiCurve. '''
		
		g = self.triangulation.id_encoding()
		for curve in self.components():
			h = g(curve).crush()  # Map forward under crushes first.
			g = h * g
		
		return g

class Curve(MultiCurve, Shortenable):
	''' A MultiCurve with a single component. '''
	def mcomponents(self):
		return [(self, 1)]
	
	def is_short(self):
		# Theorem: A curve is short iff either:
		#  - it meets T exactly twice, or
		#  - it meets every edge of T either 0 or 2 times and has one corridor [BellWebb16a].
		is_corridor = lambda triangle: sum(self(edge) for edge in triangle) == 4
		return self.weight() == 2 or (all(weight in [0, 2] for weight in self) and len([triangle for triangle in self.triangulation if is_corridor(triangle)]) == 1)
	
	def shorten_score(self, edge):
		# This relies on the following
		# Lemma: If there are no bipods then self(edge) \in [0,2] for each edge.
		#
		# This follows from the fact that if this is a tripod / monopod in every triangle then
		# connectedness implies that there can't be any part which can't see a vertex.
		
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		# Low score == bad.
		if not self.triangulation.is_flippable(edge):
			return 0
		
		a, b, c, d, e = self.triangulation.square(edge)
		ia, ib, ic, id, ie = [self(edgey) for edgey in self.triangulation.square(edge)]
		da, db, dc, dd, de = [self.dual_weight(edgey) for edgey in self.triangulation.square(edge)]
		if ie == 0:
			return 0
		if max(ia+ ic, ib + id) == ie:  # Drops to zero.
			return 3
		if de > 0:
			return 0
		
		if da > 0 and db > 0:
			return 2
		
		return 1
	
	def parallel(self):
		''' Return an edge that this curve is parallel to. '''
		assert(self.is_short())
		
		return min([edge for edge in self.triangulation.edges if self(edge) == 0 and self.dual_weight(edge) > 0], key=lambda e: e.label)  # Take the minimum of two.
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve, raised to the power k. '''
		
		# Some easy cases:
		if k == 0: return self.triangulation.id_encoding()
		if k < 0: return self.encode_twist(-k).inverse()
		
		short, conjugator = self.shorten()
		
		# TODO: 2) Make polynomial-time by taking advantage of spiralling.
		
		a = short.parallel()
		if short.weight() == 2:  # curve is non-isolating.
			num_flips = 1
		else:  # curve is isolating.
			# Theorem: 3*num_tripods is the right number of flips to do.
			# Proof: TODO.
			num_tripods = len([triangle for triangle in short.triangulation if sum(short(edge) for edge in triangle) == 6])
			num_flips = 3*num_tripods
		
		twist = short.triangulation.id_encoding()
		for _ in range(num_flips):
			twist = twist.target_triangulation.encode_flip(twist.target_triangulation.corner_lookup[a.label][2]) * twist
		twist = twist.target_triangulation.find_isometry(twist.source_triangulation, {a.label: a.label}).encode() * twist
		
		return conjugator.inverse() * twist**k * conjugator
	
	def intersection(self, lamination):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(lamination, curver.kernel.Lamination))
		assert(lamination.triangulation == self.triangulation)
		
		short, conjugator = self.shorten()
		short_lamination = conjugator(lamination)
		
		a = short.parallel()
		v = short.triangulation.vertex_lookup[a.label]  # = self.triangulation.vertex_lookup[~a.label].
		edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
		
		around = min(short_lamination.side_weight(edge) for edge in edges)
		if around > 0:  # All side_weights and edge weights are non-negative.
			return short_lamination(a) - 2 * min(short_lamination.side_weight(edge) for edge in edges)
		else:
			return short_lamination(a) - sum(min(short_lamination.side_weight(edge), 0) for edge in edges) - sum(min(short_lamination(edge), 0) for edge in edges[1:])
	
	def crush(self):
		''' Return the crush map associated to this Curve. '''
		
		short, conjugator = self.shorten()
		
		# Use the following for reference:
		#             #<----------#                #  #-----------#  #
		#            /|     a    ^|               /|  |     a    /  /|
		#           / |         / |              / |  |         /  / |
		#          /  |        /  |             /  |  |        /  /  |
		#         /   |       /   |            /   |  |       /  /   |
		#        /    |b    e/    |   ===>>   /    |  |b   ~b/  /    |
		#       /   ~b|     /~e   |          /    e|  |     /  /~e   |
		#      /      |    /      |         /      |  |    /  /      |
		#     /       |   /       |        /       |  |   /  /       |
		#    /        |  /        |       /        |  |  /  /        |
		#   /         | /         |      /         |  | /  /         |
		#  /          V/          |     /          |  |/  /          |
		# #-----------#-----------#    #-----------#  #  #-----------#
		# Where a is parallel to short.
		
		a = short.parallel()
		a, b, e = short.triangulation.corner_lookup[a.label]
		
		# Build the new triangulation.
		edge_map = dict((edge, curver.kernel.Edge(edge.label)) for edge in short.triangulation.edges)
		# Remap some edges.
		edge_map[e] = curver.kernel.Edge(~b.label)
		edge_map[~b] = curver.kernel.Edge(e.label)
		
		new_triangulation = curver.kernel.Triangulation([curver.kernel.Triangle([edge_map[edgey] for edgey in triangle]) for triangle in short.triangulation])
		
		# Build the lifting matrix back.
		v = short.triangulation.vertex_lookup[a.label]  # = short.triangulation.vertex_lookup[~a.label].
		indices = [edge.index for edge in curver.kernel.utilities.cyclic_slice(v, a, ~a)[1:]]  # The indices that appear walking around v from a to ~a. Note need to exclude the initial a.
		matrix = [[1 if i == j else indices.count(j) if i == b.index else 1 if (i == e.index and j == b.index) else 0 for i in range(self.zeta)] for j in range(self.zeta)]
		
		crush = curver.kernel.Crush(short.triangulation, new_triangulation, short, matrix).encode()
		return crush * conjugator

