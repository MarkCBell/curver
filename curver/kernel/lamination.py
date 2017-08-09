
''' A module for representing laminations on Triangulations.

Provides one class: Curve. '''

import networkx
from itertools import permutations

import curver
from curver.kernel.utilities import memoize  # Special import needed for decorating.

INFTY = float('inf')

def dual_weight(a, b, c):
	a, b, c = max(a, 0), max(b, 0), max(c, 0)  # Correct for negatives.
	correction = min(a + b - c, b + c - a, c + a - b, 0)
	return (b + c - a + correction) // 2

class Lamination(object):
	''' This represents a lamination on a triangulation.
	
	Users can use Triangulation.lamination(). '''
	def __init__(self, triangulation, geometric):
		assert(isinstance(triangulation, curver.kernel.Triangulation))
		
		self.triangulation = triangulation
		self.zeta = self.triangulation.zeta
		self.geometric = geometric
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(self.geometric)
	def __iter__(self):
		return iter(self.geometric)
	def __len__(self):
		return sum(multiplicity for _, multiplicity in self.mcomponents())  # Total number of components.
	def __call__(self, edge):
		''' Return the geometric measure assigned to item. '''
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]
		
		return self.geometric[edge.index]
	def __eq__(self, other):
		return self.triangulation == other.triangulation and self.geometric == other.geometric
	def __ne__(self, other):
		return not (self == other)
	def __hash__(self):
		return hash(tuple(self.geometric))
	def __add__(self, other):
		# Haken sum.
		if isinstance(other, Lamination):
			if other.triangulation != self.triangulation:
				raise ValueError('Laminations must be on the same triangulation to add them.')
			
			geometric = [x + y for x, y in zip(self.geometric, other.geometric)]
			return self.triangulation.lamination(geometric)  # Have to promote.
		else:
			return NotImplemented
	def __radd__(self, other):
		return self + other  # Commutative.
	def __mul__(self, other):
		geometric = [other * x for x in self]
		# TODO: 3) Could save components if they have already been computed.
		return self.__class__(self.triangulation, geometric)  # Preserve promotion.
	def __rmul__(self, other):
		return self * other  # Commutative.
	
	def weight(self):
		''' Return the geometric intersection of this lamination with its underlying triangulation. '''
		
		return sum(max(weight, 0) for weight in self)
	
	def dual_weight(self, edge):
		''' Return the number of component of this lamination dual to the given edge.
		
		Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		corner = self.triangulation.corner_lookup[edge.label]
		weights = [self(edge) for edge in corner]
		return dual_weight(*weights)
	
	def side_weight(self, edge):
		''' Return the number of component of this lamination dual to the given edge.
		
		Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		corner = self.triangulation.corner_lookup[edge.label]
		return self.dual_weight(corner[1])
	
	def dual_square(self, edge):
		# Remove?
		a, b, c, d, e = self.triangulation.square(edge)
		return [self.dual_weight(edgey) for edgey in [a, b, c, d, e, ~e]]
	
	def is_empty(self):
		''' Return if this lamination has no components. '''
		
		return not any(self)  # len(self) == 0
	
	def is_multicurve(self):
		return not self.is_empty() and all(isinstance(component, Curve) for component in self.components())
	
	def is_curve(self):
		return self.is_multicurve() and len(self) == 1
	
	def is_multiarc(self):
		return not self.is_empty() and all(isinstance(component, Arc) for component in self.components())
	
	def is_arc(self):
		return self.is_multiarc() and len(self) == 1
	
	def promote(self):
		if self.is_multicurve():
			if self.is_curve():
				self.__class__ = Curve
			else:
				self.__class__ = MultiCurve
		elif self.is_multiarc():
			if self.is_arc():
				self.__class__ = Arc
			else:
				self.__class__ = MultiArc
		return self
	
	def remove_peripheral(self):
		''' Return a new lamination with any peripheral components removed.
		
		Most functions will assume that any lamination they are given does not have any peripheral components. '''
		
		peripherals = [0] * self.zeta
		for vertex in self.triangulation.vertices:
			peripheral = max(min(self.side_weight(edge) for edge in vertex), 0)
			for edge in vertex:
				peripherals[edge.index] += peripheral
		weights = [weight - peripheral for weight, peripheral in zip(self, peripherals)]  # Remove the peripheral components.
		
		return Lamination(self.triangulation, weights)
	
	def skeleton(self):
		''' Return the lamination obtained by collapsing parallel components. '''
		
		return sum([component for component in self.components()], self.triangulation.empty_lamination())
	
	def peek_component(self):
		''' Return one component of this Lamination. '''
		
		return self.components().pop()
	
	def intersection(self, lamination):
		''' Return the geometric intersection number between this lamination and the given one. '''
		
		assert(isinstance(lamination, Lamination))
		
		return sum(m1 * m2 * max(c1.intersection(c2), 0) for c1, m1 in self.mcomponents() for c2, m2 in lamination.mcomponents())
	
	def no_common_component(self, lamination):
		''' Return that self does not share any components with the given Lamination. '''
		assert(isinstance(lamination, Lamination))
		
		components = self.components()
		return not any(component in components for component, _ in other.components())
	
	def train_track(self):
		# Discards all arcs parallel to edges.
		# Now in each triangle lamination looks like one of the following types:
		# 0) Empty    # 1) One arc  # 2) Two arcs  # 3) Three arcs
		#     /\      #     /\      #     /\       #     /\
		#    /  \     #    /  \     #    /  \      #    /--\
		#   /    \    #   /\   \    #   /\  /\     #   /\  /\
		#  /      \   #  /  |   \   #  /  ||  \    #  /  ||  \
		#  --------   #  --------   #  --------    #  --------
		#
		# 0a), 1a) or 2a); which are the same as 0), 1) and 2) but with an extra arc. For example, 2a):
		#     /|\
		#    / | \
		#   /\ | /\
		#  /  |||  \
		#  ---------
		# These cases are determined by the fact that (exactly) one of their dual weights is negative.
		
		# We introduce new edges to subdivide a triangle (p, q, r) as follows:
		#            /^\
		#           / | \
		#          /  |  \
		#         /   |s(i)
		#        /    |    \
		#     r /    / \    \ q
		#      /   /     \   \
		#     /  /t(j) u(k)\  \
		#    /</             \>\
		#   /-------------------\
		#             p
		# So that afterwards every complementary region can reach a vertex.
		
		geometric = list(self.geometric)
		triangles = []
		num_subdivided = 0  # Number of subdivided triangles.
		for triangle in self.triangulation:
			dual_weights = [self.dual_weight(label) for label in triangle.labels]
			if all(weight > 0 for weight in dual_weights):  # Type 3).
				p, q, r = [curver.kernel.Edge(label) for label in triangle.labels]
				i, j, k = range(self.zeta + 3*num_subdivided, self.zeta + 3*num_subdivided + 3)
				s, t, u = [curver.kernel.Edge(i), curver.kernel.Edge(j), curver.kernel.Edge(k)]
				triangles.append(curver.kernel.Triangle([p, ~u, t]))
				triangles.append(curver.kernel.Triangle([q, ~s, u]))
				triangles.append(curver.kernel.Triangle([r, ~t, s]))
				num_subdivided += 1
				
				geometric.extend(dual_weights)  # Record intersections with new edges.
			else:
				p, q, r = [curver.kernel.Edge(label) for label in triangle.labels]
				triangles.append(curver.kernel.Triangle([p, q, r]))
		
		T = curver.kernel.Triangulation(triangles)
		return curver.kernel.TrainTrack(T, geometric)
	
	def components(self):
		''' Return the set of Arcs and Curves that appear within self. '''
		
		return set(component for component, _ in self.mcomponents())
	
	@memoize
	def mcomponents(self):
		''' Return a set of pairs (component, multiplicity). '''
		
		components = []
		for component, multiplicity in self.train_track().mcomponents():
			# Project an Arc or Curve on T back to self.triangulation.
			if isinstance(component, Curve):
				components.append((Curve(self.triangulation, component.geometric[:self.zeta]), multiplicity))
			elif isinstance(component, Arc):
				components.append((Arc(self.triangulation, component.geometric[:self.zeta]), multiplicity))
			else:
				raise ValueError('')
		
		return components
	
	def sublaminations(self):
		''' Return all sublaminations that appear within self. '''
		components = self.components()
		return [sum(sub, self.triangulation.empty_lamination()) for i in range(len(components)) for sub in permutations(components, i)]
	
	def multiarc(self):
		''' Return the maximal MultiArc contained within this lamination. '''
		
		empty = self.triangulation.empty_lamination()
		return sum([multiplicity * component for component, multiplicity in self.mcomponents() if isinstance(component, Arc)], empty)
	
	def multicurve(self):
		''' Return the maximal MultiCurve contained within this lamination. '''
		
		empty = self.triangulation.empty_lamination()
		return sum([multiplicity * component for component, multiplicity in self.mcomponents() if isinstance(component, Curve)], empty)
	
	def boundary(self):
		''' Return the boundary of a regular neighbourhood of this lamination. '''
		
		multiarc = self.multiarc()  # Might be empty.
		multicurve = self.multicurve()  # Might be empty.
		return multicurve + (multiarc if multiarc.is_empty() else multiarc.boundary)

class MultiCurve(Lamination):
	''' A Lamination in which every component is a Curve. '''
	def is_multicurve(self):
		return True
	def is_multiarc(self):
		return False
	
	def encode_twist(self, k=1):
		h = self.triangulation.id_encoding()
		for curve, multiplicity in self.mcomponents():
			h = curve.encode_twist(k * multiplicity) * h
		
		return h
	
	def boundary_union(self, other):
		''' Return \\partial N(self \\cup other). '''
		assert(isinstance(other, Lamination))
		
		crush = self.crush()
		lift = crush.inverse()
		other_prime = crush(other).skeleton()
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

class Curve(MultiCurve):
	''' A MultiCurve with a single component. '''
	def mcomponents(self):
		return [(self, 1)]
	
	def score(self, edge):
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
			return 10
		if de > 0:
			return 0
		
		if da > 0 and db > 0:
			return 2
		
		return 1
	
	def shorten(self):
		''' Return an encoding which maps this curve to a curve to a short one, one with as little weight as possible, together with its image. '''
		
		# Repeatedly flip to reduce the weight of this curve as much as possible.
		# TODO: 3) Make polynomial-time by taking advantage of spiralling.
		
		# This relies on the following
		# Lemma: If there are no bipods then self(edge) \in [0,2] for each edge.
		#
		# This follows from the fact that if this is a tripod / monopod in every triangle then
		# connectedness implies that there can't be any part which can't see a vertex.
		
		curve = self
		conjugator = curve.triangulation.id_encoding()
		
		extra = []
		while not curve.is_short():
			edge = max(extra + curve.triangulation.edges, key=curve.score)
			# This edge is always flippable.
			
			move = curve.triangulation.encode_flip(edge)
			conjugator = move * conjugator
			curve = move(curve)
			
			# TODO: 3) Accelerate!!
			a, b, c, d, e = curve.triangulation.square(edge)
			extra = [a, d]
		
		return curve, conjugator
	
	def is_short(self):
		# Theorem: A curve is short iff either:
		#  - it meets T exactly twice, or
		#  - it meets every edge of T either 0 or 2 times and has one corridor [BellWebb16a].
		is_corridor = lambda triangle: sum(self(edge) for edge in triangle) == 4
		return self.weight() == 2 or (all(weight in [0, 2] for weight in self) and len([triangle for triangle in self.triangulation if is_corridor(triangle)]) == 1)
	
	def parallel(self):
		
		assert(self.is_short())
		
		return min([edge for edge in self.triangulation.edges if self(edge) == 0 and self.dual_weight(edge) > 0], key=lambda e: e.label)  # Take the minimum of two.
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve, raised to the power k.
		
		Currently, this must be a non-isolating curve. '''
		
		if k == 0: return self.triangulation.id_encoding()
		
		short, conjugator = self.shorten()
		
		if short.weight() == 2:  # curve is non-isolating.
			# Get some edges.
			a = short.parallel()
			_, b, e = short.triangulation.corner_lookup[a.label]
			_, c, d = short.triangulation.corner_lookup[~e.label]
			
			# We now have:
			# #<----------#
			# |     a    ^^
			# |         / |
			# |---->------|
			# |       /   |
			# |b    e/   d|
			# |     /     |
			# |    /      |
			# |   /       |
			# |  /        |
			# | /         |
			# V/    c     |
			# #---------->#
			# where b == ~d.
			
			twist = short.triangulation.encode([{i: i for i in short.triangulation.indices if i not in [e.index, b.index]}, e.label])
			
			# TODO: 4) Once Spiral is working we can do:
			# twist_k = triangulation.encode([(e.label, k)])
			# return conjugator.inverse() * twist_k * conjugator
		else:  # curve is isolating.
			a = short.parallel()
			v = short.triangulation.vertex_lookup[a.label]
			b = short.triangulation.corner_lookup[a.label][1]
			
			# Build the image of the most twisted up arc.
			geometric = [2 * weight for weight in short]
			for edge in curver.kernel.utilities.cyclic_slice(v, ~b, b) + (b,):
				geometric[edge.index] -= 1
			arc = Arc(short.triangulation, geometric)
			
			twist = short.triangulation.id_encoding()
			while not arc.is_short():
				flip = arc.triangulation.encode_flip(arc.triangulation.corner_lookup[a.label][2])
				twist = flip * twist
				arc = flip(arc)
			
			# Relabel back.
			[this_component] = [component for component in twist.target_triangulation.components() if a in component]
			label_map = dict(
				[(edge.label, edge.label) for edge in twist.target_triangulation.edges if edge not in this_component] + \
				[(a.label, a.label)]
				)
			twist = twist.target_triangulation.find_isometry(twist.source_triangulation, label_map).encode() * twist
		
		return conjugator.inverse() * twist**k * conjugator
	
	def intersection(self, other):
		''' Return the geometric intersection between self and the given other.
		
		Currently assumes (and checks) that self is a non-isolating curve. '''
		
		assert(isinstance(other, (Arc, Curve)))
		if isinstance(other, Arc):
			return other.intersection(self)
		
		assert(other.triangulation == self.triangulation)
		
		short, conjugator = self.shorten()
		short_other = conjugator(other)
		
		a = short.parallel()
		v = self.triangulation.vertex_lookup[a.label]  # = self.triangulation.vertex_lookup[~a.label].
		edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
		
		# This assumes that other is a curve.
		return short_other(a) - 2 * min(short_other.side_weight(edge) for edge in edges)
	
	def crush(self):
		''' Return the crush map associated to this Curve.
		
		Currently assumes that this a non-isolating curve. '''
		
		short, conjugator = self.shorten()
		
		if short.weight() == 2:  # curve is non-isolating.
			triangulation = short.triangulation
			# Grab the indices of the two edges we meet.
			e1, e2 = [edge_index for edge_index in short.triangulation.indices if short(edge_index) > 0]
			
			# We might need to swap these edge indices so we have a good frame of reference.
			if short.triangulation.corner_lookup[e1].indices[2] != e2: e1, e2 = e2, e1
			
			a, b, c, d, e = triangulation.square(e1)
			
			# Use the following for reference:
			# #<----------#     #-----------#  #
			# |     a    ^^     |     a    /  /|
			# |         / |     |         /  / |
			# |        /  |     |        /  /  |
			# |       /   |     |       /  /   |
			# |b    e/   d| --> |b   ~b/  /~e e|
			# |     /     |     |     /  /     |
			# |    /      |     |    /  /      |
			# |   /       |     |   /  /       |
			# |  /        |     |  /  /        |
			# | /         |     | /  /         |
			# V/    c     |     |/  /    c     |
			# #---------->#     #  #-----------#
			
			edge_map = dict((edge, curver.kernel.Edge(edge.label)) for edge in triangulation.edges)
			
			# Most triangles don't change.
			triangles = [curver.kernel.Triangle([edge_map[edgey] for edgey in triangle]) for triangle in triangulation if e not in triangle and ~e not in triangle]
			
			triangle_A2 = curver.kernel.Triangle([edge_map[a], edge_map[b], edge_map[~b]])
			triangle_B2 = curver.kernel.Triangle([edge_map[c], edge_map[e], edge_map[~e]])
			new_triangulation = curver.kernel.Triangulation(triangles + [triangle_A2, triangle_B2])
			
			matrix = [[1 if i == j or (i == b.index and j == e.index) or (i == e.index and j == b.index) else 0 for i in range(self.zeta)] for j in range(self.zeta)]
		else:  # curve is isolating.
			raise curver.AssumptionError('Curve is isolating.')  # TODO: 4) Handle isolating case.
		
		crush = curver.kernel.Crush(triangulation, new_triangulation, short, matrix).encode()
		return crush * conjugator


class MultiArc(Lamination):
	''' A Lamination in which every component is an Arc. '''
	def is_multicurve(self):
		return False
	def is_multiarc(self):
		return True
	
	def score(self, edge):
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		# Low score == bad.
		if not self.triangulation.is_flippable(edge):
			return 0
		if self.dual_weight(edge) < 0:
			return 1
		
		return 0
	
	def shorten(self):
		''' Return an encoding which maps this arc to one with as little weight as possible together with its image. '''
		
		# TODO: 3) Make polynomial-time by taking advantage of spiralling.
		
		arc = self
		conjugator = arc.triangulation.id_encoding()
		extra = []
		while not arc.is_short():
			edge = max(extra + arc.triangulation.edges, key=arc.score)
			# This edge is always flippable.
			
			move = arc.triangulation.encode_flip(edge)
			conjugator = move * conjugator
			curve = move(arc)
			
			# TODO: 3) Accelerate!!
			a, b, c, d, e = arc.triangulation.square(edge)
			extra = [a, d]
		
		return arc, conjugator
	
	def is_short(self):
		return self.weight() == 0
	
	def boundary(self):
		''' Return the multicurve which is the boundary of a regular neighbourhood of this multiarc. '''
		
		short, conjugator = self.shorten()
		# short is a subset of the edges of the triangulation it is defined on.
		# So its geometric vector is non-positive.
		
		geometric = [0 if weight < 0 else 2 for weight in short]
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
		assert(self.is_short())
		
		return min([edge for edge in self.triangulation.edges if self(edge) < 0], key=lambda e: e.label)  # Take the minimum of two.
	
	def encode_halftwist(self, k=1):
		''' Return an Encoding of a left half twist about a regular neighbourhood of this arc, raised to the power k.
		
		Assumes (and checks) that this arc connects between distinct vertices. '''
		# Will Worden checked that this works for genus <= 20.
		
		boundary = self.boundary()
		if len(boundary) != 1:
			raise curver.AssumptionError('Arc connects a vertex to itself.')
		
		# TODO: 4) Once we can twist about any curve we can go back to:
		# if k % 2 == 0:  # k is even so use a Dehn twist about the boundary.
		#	return boundary.encode_twist(k // 2)
		
		# We need to get to a really good configuration, one where self
		# is an edge of the triangulation whose valence(initial vertex) == 1.
		#
		# We achieve this in two steps. First conjugate to make self an edge of some triangulation.
		short, conjugator = self.shorten()
		arc = short.parallel()
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
		# and on this component (the one containing this arc) to invert this arc.
		[this_component] = [component for component in short.triangulation.components() if arc in component]
		label_map = dict(
			[(edge.label, edge.label) for edge in short.triangulation.edges if edge not in this_component] + \
			[(arc.label, ~arc.label)]
			)
		half_twist = short.triangulation.find_isometry(half_twist.source_triangulation, label_map).encode() * half_twist
		
		# Until all twists are working we must do:
		return conjugator.inverse() * half_twist**k * conjugator
		
		# TODO: 4) Afterwards we can go back to:
		# We accelerate large powers by replacing (T^1/2_self)^2 with T_self which includes acceleration.
		if abs(k) == 1:
			return conjugator.inverse() * half_twist**k * conjugator
		else:  # k is odd so we need to add in an additional half twist.
			# Note: k // 2 always rounds down, so even if k < 0 the additional half twist we need to do is positive.
			return boundary.encode_twist(k // 2) * conjugator.inverse() * half_twist * conjugator
	
	def intersection(self, other):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(other, Lamination))
		assert(other.triangulation == self.triangulation)
		
		# short = [0, 0, ..., 0, -1, 0, ..., 0]
		short, conjugator = self.shorten()
		short_other = conjugator(other)
		
		return sum(b for a, b in zip(short, short_other) if a == -1 if b >= 0)

