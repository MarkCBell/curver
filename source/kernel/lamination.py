
''' A module for representing laminations on Triangulations.

Provides one class: Curve. '''

import curver

INFTY = float('inf')

def sign(number):
	return 1 if number > 0 else 0 if number  == 0 else -1  # if number < 0.

class Lamination(object):
	''' This represents a lamination on an triangulation.
	
	It is given by a list of its geometric intersection numbers.
	Users should use Triangulation.lamination() to create laminations. '''
	def __init__(self, triangulation, geometric):
		assert(isinstance(triangulation, curver.kernel.Triangulation))
		assert(isinstance(geometric, (list, tuple)))
		# We should check that geometric satisfies reasonable relations.
		
		self.triangulation = triangulation
		self.zeta = self.triangulation.zeta
		self.geometric = list(geometric)
		assert(len(self.geometric) == self.zeta)
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(self.geometric)
	def __iter__(self):
		return iter(self.geometric)
	
	def __call__(self, item):
		''' Return the geometric measure assigned to item. '''
		if isinstance(item, curver.kernel.Edge):
			return self.geometric[item.index]
		else:
			return self.geometric[curver.kernel.norm(item)]
	
	def __len__(self):
		return self.zeta
	
	def __eq__(self, other):
		return self.triangulation == other.triangulation and \
			all(v == w for v, w in zip(self.geometric, other.geometric))
	def __ne__(self, other):
		return not (self == other)
	
	def __hash__(self):
		# This should be done better.
		return hash(tuple(self.geometric))
	
	def __add__(self, other):
		if isinstance(other, Lamination):
			if other.triangulation != self.triangulation:
				raise ValueError('Laminations must be on the same triangulation to add them.')
			
			# Haken sum.
			geometric = [x + y for x, y in zip(self, other)]
			return Lamination(self.triangulation, geometric)
		else:
			if other == 0:  # So we can use sum.
				return self
			else:
				return NotImplemented
	def __radd__(self, other):
		return self + other
	def __sub__(self, other):
		if isinstance(other, Lamination):
			if other.triangulation != self.triangulation:
				raise ValueError('Laminations must be on the same triangulation to subtract them.')
			
			geometric = [x - y for x, y in zip(self, other)]
			return Lamination(self.triangulation, geometric)
		else:
			return NotImplemented
	def __mul__(self, other):
		geometric = [other * x for x in self]
		return Lamination(self.triangulation, geometric)
	def __rmul__(self, other):
		return self * other
	
	def is_empty(self):
		''' Return if this lamination is equal to the empty lamination. '''
		
		return not any(self)
	
	def isometries_to(self, other):
		''' Return a list of isometries taking this lamination to other. '''
		
		assert(isinstance(other, Lamination))
		
		return [isom for isom in self.triangulation.isometries_to(other.triangulation) if other== isom.encode()(self)]
	
	def self_isometries(self):
		''' Return a list of isometries taking this lamination to itself. '''
		
		return self.isometries_to(self)
	
	def weight(self):
		''' Return the geometric intersection of this lamination with its underlying triangulation. '''
		
		return sum(abs(x) for x in self.geometric)
	
	def conjugate_short(self):
		''' Return an encoding which maps this lamination to a lamination with as little weight as possible. '''
		
		# Repeatedly flip to reduce the weight of this lamination as much as possible.
		# Needs to be made polynomial-time by taking advantage of spiralling.
		
		lamination = self
		conjugation = lamination.triangulation.id_encoding()
		
		def weight_change(lamination, edge_index):
			''' Return how much the weight would change by if this flip was done. '''
			
			if lamination(edge_index) == 0 or not lamination.triangulation.is_flippable(edge_index): return INFTY
			a, b, c, d = lamination.triangulation.square_about_edge(edge_index)
			return max(lamination(a) + lamination(c), lamination(b) + lamination(d)) - 2 * lamination(edge_index)
		
		weight = lamination.weight()
		old_weight = weight + 1
		old_old_weight = old_weight + 1
		possible_edges = lamination.triangulation.indices
		drops = sorted([(weight_change(lamination, i), i) for i in possible_edges])
		# If we ever fail to make progress more than once then the lamination is as short as it's going to get.
		while weight < old_old_weight:
			# Find the edge which decreases our weight the most.
			# If none exist then it doesn't matter which edge we flip, so long as it meets the lamination.
			_, edge_index = drops[0]
			
			forwards = lamination.triangulation.encode_flip(edge_index)
			conjugation = forwards * conjugation
			lamination = forwards(lamination)
			weight, old_weight, old_old_weight = lamination.weight(), weight, old_weight
			
			# Update new neighbours.
			changed_indices = set([edge.index for edge in lamination.triangulation.square_about_edge(edge_index)] + [edge_index])
			drops = sorted([(d, i) if i not in changed_indices else (weight_change(lamination, i), i) for (d, i) in drops])  # This should be lightning fast since the list was basically sorted already.
			# If performance really becomes an issue then we could look at using heapq.
		
		return lamination, conjugation
	
	def dual_weights(self, triangle):
		''' Return the number of component of this lamination dual to each side of the given triangle.
		
		Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
		
		weights = [self(index) for index in triangle.indices]
		correction = min(weights[(i+1)%3] + weights[(i+2)%3] - weights[i] for i in range(3))
		if correction >= 0:  # No terminal arcs.
			return [(weights[(i+1)%3] + weights[(i+2)%3] - weights[i]) // 2 for i in range(3)]
		else:  # Terminal arc.
			return [(weights[(i+1)%3] + weights[(i+2)%3] - weights[i] + correction) // 2 for i in range(3)]
	
	def arc_side(self, triangle):
		dual_weights = self.dual_weights(triangle)
		arc = min(dual_weights)
		if arc >= 0:  # No arc.
			return None
		else:
			return dual_weights.index(arc)
	
	def components(self):
		''' Return a dictionary mapping the Curves and Arcs that appear in this lamination to their multiplicities. '''
		
		# Probably should cache.
		
		return NotImplemented
		
		comp = dict()
		
		lamination = self
		
		# Remove all the obvious arcs. This reduces the number of cases we have to consider later.
		for index in lamination.triangulation.indices:
			if lamination(index) < 0:
				geometric = [0 if i != index else -1 for i in range(self.zeta)]
				component, multiplicity = Arc(lamination.triangulation, geometric), abs(lamination(index))
				comp[arc] = multiplicity
				lamination = lamination - multiplicity * arc
		
		# Remove all the obvious curves. This reduces the number of places we have to look for new curves later.
		for index in lamination.triangulation.indices:
			if lamination.triangulation.is_flippable(index):
				a, b, c, d = lamination.triangulation.square_about_edge(index)
				if b == ~d:
					multiplicity = ??
					if multiplicity > 0:
						curve = Curve(lamination.triangulation, geometric)
						comp[project(encoding.inverse()(arc))] = multiplicity
						lamination = lamination - multiplicity * curve
		
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
		
		# We will subdivide the triangles so that afterwards each triangle contains at most one switch.
		# It is only necessary to subdivide 3) and 2a) but it's easier to do them all.
		
		
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
		# WLOG: If there is an terminal normal arc in this triangle then it goes through side p.
		
		geometric = lamination.geometric + [None] * (2*zeta)
		triangles = []
		for count, triangle in enumerate(lamination.triangulation):
			dual_weights = self.dual_weights(triangle)
			arc_side = self.arc_side(triangle)
			if arc_side is not None:  # Rotate the triangle so that arc_side is side #0.
				triangle = lamination.triangulation.corner_lookup[triangle.labels[arc_side]]
				dual_weights = dual_weights[arc_side:] + dual_weights[:arc_side]
				arc_side = 0
			
			i, j, k = self.zeta + 3*count, self.zeta + 3*count + 1, self.zeta + 3*count + 2
			p, q, r = [curver.kernel.Edge(label) for label in triangle.labels]
			s, t, u = [curver.kernel.Edge(i), curver.kernel.Edge(j), curver.kernel.Edge(k)]
			triangles.append(curver.kernel.Triangle([p, ~u, t]))
			triangles.append(curver.kernel.Triangle([q, ~s, u]))
			triangles.append(curver.kernel.Triangle([r, ~t, s]))
			
			# Record intersections with new edges.
			if arc_side is None:
				geometric[i], geometric[j], geometric[k] = dual_weights
			else:
				pass
		
		T = curver.kernel.Triangulation(triangles)
		lamination = Lamination(T, geometric)
		
		def project(lamination):
			''' Project a good lamination on T back to one on self.triangulation. '''
			assert(lamination.triangulation == T)
			return Lamination(self.triangulation, lamination.geometric[:self.zeta])
		
		to_flip = None
		encoding = T.id_encoding()
		while not lamination.is_empty():
			if to_flip is None:
				???
			
			move = lamination.triangulation.encode_flip(to_flip)
			encoding = encoding * move
			lamination = move(lamination)
			
			# The only place where an obvious arc could appear is along the new edge we have just introduced.
			if lamination(to_flip) < 0:
				geometric = [0 if i != to_flip else -1 for i in range(lamination.zeta)]
				component, multiplicity = Arc(lamination.triangulation, geometric), abs(lamination(to_flip))
				comp[project(encoding.inverse()(arc))] = multiplicity
				lamination = lamination - multiplicity * arc
			
			# The only places where a short curve could appear is across an edge adjacent to the one we just flipped.
			for edge in lamination.triangulation.square_about_edge(to_flip):
				if lamination.triangulation.is_flippable(edge.label):
					a, b, c, d = lamination.triangulation.square_about_edge(edge.label)
					if b == ~d:
						multiplicity = ??
						if multiplicity > 0:
							curve = Curve(lamination.triangulation, geometric)
							comp[project(encoding.inverse()(arc))] = multiplicity
							lamination = lamination - multiplicity * curve
			
			# Accelerate!!
			
			to_flip = ??
		
		return comp
	
	def skeleton(self):
		''' Return the lamination obtained by collapsing parallel components. '''
		
		return sum(self.components())
	
	def is_multicurve(self):
		return all(isinstance(component, Curve) for component in self.components())
	
	def is_curve(self):
		return self.is_multicurve() and sum(self.components().values()) == 1
	
	def is_multiarc(self):
		return all(isinstance(component, Arc) for component in self.components())
	
	def is_arc(self):
		return self.is_multiarc() and sum(self.components().values()) == 1
	
	def promote(self):
		if self.is_multicurve():
			if self.is_curve():
				self.__class__ = Curve
			else:
				self.__class__ = MultiCurve
		if self.is_multiarc():
			if self.is_arc():
				self.__class__ = Arc
			else:
				self.__class__ = MultiArc
		return self
	
	def intersection(self, lamination):
		''' Return the geometric intersection number between this lamination and the given one. '''
		
		return sum(multiplicity * component.intersection(lamination) for component, multiplicity in self.components().items())


class MultiCurve(Lamination):
	def is_multicurve(self):
		return True
	def is_multiarc(self):
		return False

class Curve(MultiCurve):
	def is_curve(self):
		return True
	def components(self):
		return {self: 1}
	
	def is_isolating(self):
		''' Return if this curve is an isolating curve.
		
		That is, is there a component of S - self that does not contain a puncture. '''
		
		# This is based off of self.encode_twist(). See the documentation there as to why this works.
		
		short, _ = self.conjugate_short()
		
		return short.weight() == 2
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve raised to the power k.
		
		Currently, this must be a non-isolating curve. '''
		
		assert(not self.is_isolating())
		
		if k == 0: return self.triangulation.id_encoding()
		
		short, conjugation = self.conjugate_short()
		
		triangulation = short.triangulation
		# Grab the indices of the two edges we meet.
		e1, e2 = [edge_index for edge_index in short.triangulation.indices if short(edge_index) > 0]
		
		a, b, c, d = triangulation.square_about_edge(e1)
		# If the curve is going vertically through the square then ...
		if short(a) == 1 and short(c) == 1:
			# swap the labels round so it goes horizontally.
			e1, e2 = e2, e1
			a, b, c, d = triangulation.square_about_edge(e1)
		elif short(b) == 1 and short(d) == 1:
			pass
		
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
		# And e.index = e1 and b.index = d.index = e2.
		
		# Used to do:
		# twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [e1, e2]}, e1]))
		# return conjugation.inverse() * twist**k * conjugation
		
		twist_k = triangulation.encode([(e1, k)])
		return conjugation.inverse() * twist_k * conjugation
	
	def is_halftwistable(self):
		''' Return if this curve is a half twistable curve. '''
		
		# This is based off of self.encode_halftwist(). See the documentation there as to why this works.
		
		short, _ = self.conjugate_short()
		
		triangulation = short.triangulation
		
		e1, e2 = [edge_index for edge_index in short.triangulation.indices if short(edge_index) > 0]
		
		a, b, c, d = triangulation.square_about_edge(e1)
		if short(a) == 1 and short(c) == 1:
			e1, e2 = e2, e1
			a, b, c, d = triangulation.square_about_edge(e1)
		elif short(b) == 1 and short(d) == 1:
			pass
		
		_, _, z, w = triangulation.square_about_edge(a.label)
		_, _, x, y = triangulation.square_about_edge(c.label)
		
		return z == ~w or x == ~y
	
	
	def intersection(self, lamination):
		''' Return the geometric intersection between self and the given lamination.
		
		Currently assumes (and checks) that self is a non-isolating curve. '''
		
		assert(isinstance(lamination, Lamination))
		assert(lamination.triangulation == self.triangulation)
		
		if self.is_isolating():
			raise curver.AssumptionError('Can only compute geometric intersection number between a non-isolating curve and a lamination.')
		
		short_self, conjugator = self.conjugate_short()
		short_lamination = conjugator(lamination)
		
		triangulation = short_self.triangulation
		e1, e2 = [edge_index for edge_index in triangulation.indices if short_self(edge_index) > 0]
		# We might need to swap these edge indices so we have a good frame of reference.
		if triangulation.corner_lookup[e1].indices[2] != e2: e1, e2 = e2, e1
		
		a, b, c, d = triangulation.square_about_edge(e1)
		e = e1
		
		x = (short_lamination(a) + short_lamination(b) - short_lamination(e)) // 2
		y = (short_lamination(b) + short_lamination(e) - short_lamination(a)) // 2
		z = (short_lamination(e) + short_lamination(a) - short_lamination(b)) // 2
		x2 = (short_lamination(c) + short_lamination(d) - short_lamination(e)) // 2
		y2 = (short_lamination(d) + short_lamination(e) - short_lamination(c)) // 2
		z2 = (short_lamination(e) + short_lamination(c) - short_lamination(d)) // 2
		
		intersection_number = short_lamination(a) - 2 * min(x, y2, z)
		
		# Check that the other formula gives the same answer.
		assert(intersection_number == short_lamination(c) - 2 * min(x2, y, z2))
		
		return intersection_number
	
	def quasiconvex(self, other):
		''' Return a polynomial-sized quasiconvex subset of the curve complex that contains self and other. '''
		
		return NotImplemented
	
	def geodesic(self, other):
		''' Return a geodesic in the curve complex from self to other.
		
		The geodesic will always come from a tight geodesic. '''
		
		assert(isinstance(other, Curve))
		
		return NotImplemented
	
	def distance(self, other):
		''' Return the distance from self to other in the curve complex. '''
		
		return len(self.geodesic(other))
	
	def crush(self):
		''' Return the crush map. '''
		
		return NotImplemented


class MultiArc(Lamination):
	def is_multicurve(self):
		return False
	def is_multiarc(self):
		return True
	
	def boundary(self):
		''' Return the multicurve which is the boundary of a regular neighbourhood of this multiarc. '''
		
		short, conjugator = self.conjugate_short()
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
		
		boundary = short.triangulation.lamination(geometric).promote()
		
		return conjugator.inverse()(boundary)

class Arc(MultiArc):
	def is_arc(self):
		return True
	def components(self):
		return {self: 1}
	
	def encode_halftwist(self, k=1):
		''' Return an Encoding of a left half twist about a regular neighbourhood of this arc raised to the power k.
		
		Assumes (and checks) that this arc connects between distinct vertices.
		Currently, the boundary curve must also be non-isolating. '''
		
		boundary = self.boundary()
		if not isinstance(boundary, Curve):
			raise curver.AssumptionError('Arc connects a vertex to itself.')
		
		if k % 2 == 0:  # k is even so use a Dehn twist about the boundary.
			return boundary.encode_twist(k // 2)
		
		short_boundary, conjugation = boundary.conjugate_short()
		short = conjugation(self)
		
		if not short_boundary.weight() == 2:  # boundary.is_isolating().
			raise curver.AssumptionError('Boundary curve is isolating.')
		
		[arc_index] = [index for index in short.triangulation.indices if short(index) != 0]
		triangulation = short.triangulation
		
		# Let x be the edge of triangulation with index arc_index. We build out the neighbourhood
		# of x which, since boundary is short, looks like the following:
		#
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
		#  \         /
		#   \       /
		#    \x   y/
		#     \   /
		#      \ /
		#       #
		# Where d == ~b and x == ~y.
		
		c = ~triangulation.folded_boundary(arc_index)
		d, tilde_e, _, _ = triangulation.square_about_edge(c.label)
		e = ~tilde_e
		a, b, c, d = triangulation.square_about_edge(e.label)
		
		half_twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [b.index, e.index, c.index, x.index]}, b.index, e.index, c.index])
		
		# We accelerate large powers by replacing (T^1/2_self)^2 with T_self which includes acceleration.
		if abs(k) == 1:
			return conjugation.inverse() * half_twist**k * conjugation
		else:  # k is odd so we need to add in an additional half twist.
			# Note: k // 2 always rounds down, so even if k < 0 the additional half twist we need to do is positive.
			return conjugation.inverse() * short_boundary.encode_twist(k // 2) * half_twist * conjugation
	
	
	
	
	def intersection(self, lamination):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(lamination, Lamination))
		assert(lamination.triangulation == self.triangulation)
		
		# short_self = [0, 0, ..., 0, -1, 0, ..., 0]
		short_self, conjugator = self.conjugate_short()
		short_lamination = conjugator(lamination)
		
		return sum(b for a, b in zip(short_self, short_lamination) if a == -1 and b >= 0)
	
	

