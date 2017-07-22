
''' A module for representing laminations on Triangulations.

Provides one class: Curve. '''

import curver

INFTY = float('inf')

def sign(number):
	return 1 if number > 0 else 0 if number  == 0 else -1  # if number < 0.

class Lamination(object):
	''' This represents a lamination on an triangulation.
	
	It is given by a list of its geometric intersection numbers and a
	list of its algebraic intersection numbers with the (oriented) edges
	of underlying triangulation. Note that:
	     ^L
	     |
	-----|------> e
	     |
	has algebraic intersection +1.
	
	For an arc that runs parallel to an edge, we say that its algebraic intersection
	number is +1 if it runs in the same direction, that is:
	 --------------> e
	  \----------> L
	has algebraic intersection +1.
	
	Users should use Triangulation.lamination() to create laminations. '''
	def __init__(self, triangulation, geometric, algebraic):
		assert(isinstance(triangulation, curver.kernel.Triangulation))
		assert(isinstance(geometric, (list, tuple)))
		assert(isinstance(algebraic, (list, tuple)))
		# We should check that geometric / algebraic satisfies reasonable relations.
		
		self.triangulation = triangulation
		self.zeta = self.triangulation.zeta
		self.geometric = list(geometric)
		self.algebraic = list(algebraic)
		assert(len(self.geometric) == self.zeta)
		assert(len(self.algebraic) == self.zeta)
	
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
	
	def __getitem__(self, item):
		''' Return the algebraic measure assigned to item. '''
		if isinstance(item, curver.kernel.Edge):
			return self.algebraic[item.index] * item.sign()
		else:
			normalised = curver.kernel.norm(item)
			return self.algebraic[normalised] * (1 if normalised == item else -1)
	
	def __len__(self):
		return self.zeta
	
	def __eq__(self, other):
		return self.triangulation == other.triangulation and \
			all(v == w for v, w in zip(self.geometric, other.geometric)) and \
			all(v == w for v, w in zip(self.algebraic, other.algebraic))
	def __ne__(self, other):
		return not (self == other)
	
	def __hash__(self):
		# This should be done better.
		return hash(tuple(self.geometric) + tuple(self.algebraic))
	
	def __add__(self, other):
		if isinstance(other, Lamination):
			if other.triangulation != self.triangulation:
				raise ValueError('Laminations must be on the same triangulation to add them.')
			
			# Haken sum.
			geometric = [x + y for x, y in zip(self, other)]
			algebraic = [x + y for x, y in zip(self.algebraic, other.algebraic)]
			return Lamination(self.triangulation, geometric, algebraic)
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
			
			# Haken sum.
			geometric = [x - y for x, y in zip(self, other)]
			algebraic = [x - y for x, y in zip(self.algebraic, other.algebraic)]
			return Lamination(self.triangulation, geometric, algebraic)
		else:
			return NotImplemented
	def __mul__(self, other):
		geometric = [other * x for x in self]
		algebraic = [other * x for x in self]
		return Lamination(self.triangulation, geometric, algebraic)
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
		''' Returns a list of isometries taking this lamination to itself. '''
		
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
				algebraic = [0 if i != index else sign(lamination[i]) for i in range(self.zeta)]
				component,multiplicity = Arc(lamination.triangulation, geometric, algebraic), abs(lamination(index))
				comp[arc] = multiplicity
				lamination = lamination - multiplicity * arc
		
		# Now in each triangle lamination looks like one of the following types:
		# 0) Empty    # 1) One arc  # 2) Two arcs  # 3) Three arcs
		#     /\      #     /\      #     /\       #     /\
		#    /  \     #    /  \     #    /  \      #    /--\
		#   /    \    #   /\   \    #   /\  /\     #   /\  /\
		#  /      \   #  /  |   \   #  /  ||  \    #  /  ||  \
		#  --------   #  --------   #  --------    #  --------
		#
		# 0a), 1a) or 2a) which are the same as 0), 1) and 2) but with an extra arc. For example, 2a):
		#     /|\
		#    / | \
		#   /\ | /\
		#  /  |||  \
		#  ---------
		
		# We will subdivide the type 3) triangles. This type is determined by the fact that all of its
		# dual weights (a + b - c) / 2 are positive.
		
		def subdivide(triangle):
			weights = [self(index) for index in triangle.indices]
			return all(weights[(i+1)%3] + weights[(i+2)%3] - weights[i] > 0 for i in range(3))
		
		zeta = lamination.zeta  # The number of edges in the new triangulaton.
		geometric = [None] * zeta
		algebraic = [None] * zeta
		triangles = []
		for triangle in lamination.triangulation:
			if subdivide(triangle):
				p, q, r = [curver.kernel.Edge(label) for label in triangle.labels]
				s, t, u = [curver.kernel.Edge(zeta), curver.kernel.Edge(zeta+1), curver.kernel.Edge(zeta+2)]
				triangles.append(curver.kernel.Triangle([p, ~u, t]))
				triangles.append(curver.kernel.Triangle([q, ~s, u]))
				triangles.append(curver.kernel.Triangle([r, ~t, s]))
				
				# Record intersections with new edges.
				geometric.append(??)
				geometric.append(??)
				geometric.append(??)
				algebraic.append(??)
				algebraic.append(??)
				algebraic.append(??)
				
				zeta += 3
				pass
			else:
				triangles.append(curver.kernel.Triangle([curver.kernel.Edge(label) for label in triangle.labels]))
			
			for index in triangle.indices:
				geometric[index] = lamination(index)
				algebraic[index] = lamination[index]
		
		T = curver.kernel.Triangulation(triangles)
		lamination = Lamination(T, geometric, algebraic)
		
		
		def project(lamination):
			''' Project a good lamination on T back to one on self.triangulation. '''
			return Lamination(self.triangulation, lamination.geometric[:self.zeta], lamination.algebraic[:self.zeta])
		
		encoding = T.id_encoding()
		
		to_flip = None
		
		while not lamination.is_empty():
			# Remove all the obvious arcs.
			for index in lamination.triangulation.indices:
				if lamination(index) < 0:
					geometric = [0 if i != index else -1 for i in range(lamination.zeta)]
					algebraic = [0 if i != index else sign(lamination[i]) for i in range(lamination.zeta)]
					component,multiplicity = Arc(lamination.triangulation, geometric, algebraic), abs(lamination(index))
					comp[project(encoding.inverse()(arc))] = multiplicity
					lamination = lamination - multiplicity * arc
			
			if to_flip is None:
				???
			
			move = lamination.triangulation.encode_flip(to_flip)
			encoding = encoding * move
			lamination = move(lamination)
			
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
	
	def is_twistable(self):
		''' Return if this curve is a twistable curve. '''
		
		# This is based off of self.encode_twist(). See the documentation there as to why this works.
		if not self.is_curve(): return False
		
		short_curve, _ = self.conjugate_short()
		
		return short_curve.weight() == 2
	
	def is_halftwistable(self):
		''' Return if this curve is a half twistable curve. '''
		
		# This is based off of self.encode_halftwist(). See the documentation there as to why this works.
		
		short_curve, _ = self.conjugate_short()
		# We used to start with:
		#   if not self.is_twistable(): return False
		# But this wasted a lot of cycles repeating the calculation twice.
		if not short_curve.weight() == 2:
			return False
		
		triangulation = short_curve.triangulation
		
		e1, e2 = [edge_index for edge_index in short_curve.triangulation.indices if short_curve(edge_index) > 0]
		
		a, b, c, d = triangulation.square_about_edge(e1)
		if short_curve(a) == 1 and short_curve(c) == 1:
			e1, e2 = e2, e1
			a, b, c, d = triangulation.square_about_edge(e1)
		elif short_curve(b) == 1 and short_curve(d) == 1:
			pass
		
		_, _, z, w = triangulation.square_about_edge(a.label)
		_, _, x, y = triangulation.square_about_edge(c.label)
		
		return z == ~w or x == ~y
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve raised to the power k.
		
		This curve must be a twistable curve. '''
		
		assert(self.is_twistable())
		
		if k == 0: return self.triangulation.id_encoding()
		
		short_curve, conjugation = self.conjugate_short()
		
		triangulation = short_curve.triangulation
		# Grab the indices of the two edges we meet.
		e1, e2 = [edge_index for edge_index in short_curve.triangulation.indices if short_curve(edge_index) > 0]
		
		a, b, c, d = triangulation.square_about_edge(e1)
		# If the curve is going vertically through the square then ...
		if short_curve(a) == 1 and short_curve(c) == 1:
			# swap the labels round so it goes horizontally.
			e1, e2 = e2, e1
			a, b, c, d = triangulation.square_about_edge(e1)
		elif short_curve(b) == 1 and short_curve(d) == 1:
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
	
	def encode_halftwist(self, k=1):
		''' Return an Encoding of a left half twist about this curve raised to the power k.
		
		This curve must be a half-twistable curve. '''
		
		assert(self.is_halftwistable())
		
		if k % 2 == 0:  # k is even so use a Dehn twist
			return self.encode_twist(k // 2)
		
		# This first section is the same as in self.encode_flip.
		
		short_curve, conjugation = self.conjugate_short()
		
		triangulation = short_curve.triangulation
		e1, e2 = [edge_index for edge_index in short_curve.triangulation.indices if short_curve(edge_index) > 0]
		
		a, b, c, d = triangulation.square_about_edge(e1)
		# If the curve is going vertically through the square then ...
		if short_curve(a) == 1 and short_curve(c) == 1:
			# swap the labels round so it goes horizontally.
			e1, e2 = e2, e1
			a, b, c, d = triangulation.square_about_edge(e1)
		elif short_curve(b) == 1 and short_curve(d) == 1:
			pass
		
		# Get some more edges.
		_, _, z, w = triangulation.square_about_edge(a.label)
		_, _, x, y = triangulation.square_about_edge(c.label)
		
		# But now we have to go one further and worry about a, b, c, d Vs. c, d, a, b.
		# We want it so that x == ~y.
		if z.index == w.index:
			a, b, c, d = c, d, a, b
			w, x, y, z = y, z, w, x
		
		# So we now have:
		#       #
		#      / ^
		#     /   \
		#    /w   z\
		#   /       \
		#  V         \
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
		#  \         ^
		#   \       /
		#    \x   y/
		#     \   /
		#      V /
		#       #
		# Where e.index = e1 and b.index = d.index = e2,
		# and additionally x.index = y.index.
		
		half_twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [e1, e2, c.index, x.index]}, e2, e1, c.index])
		
		# We accelerate large powers by replacting (T^1/2_self)^2 with T_self which includes acceleration.
		if abs(k) == 1:
			return conjugation.inverse() * half_twist**k * conjugation
		else:  # k is odd so we need to add in an additional half twist.
			# Note: k // 2 always rounds down, so even if k < 0 the additional half twist we need to do is positive.
			return conjugation.inverse() * short_curve.encode_twist(k // 2) * half_twist * conjugation
	
	def intersection(self, lamination):
		''' Return the geometric intersection between self and the given lamination.
		
		Currently assumes (and checks) that self is a twistable curve. '''
		
		assert(isinstance(lamination, Lamination))
		assert(lamination.triangulation == self.triangulation)
		
		if not self.is_twistable():
			raise curver.AssumptionError('Can only compute geometric intersection number between a twistable curve and a curve.')
		
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

class Arc(MultiArc):
	def is_arc(self):
		return True
	
	def intersection(self, lamination):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(lamination, Lamination))
		assert(lamination.triangulation == self.triangulation)
		
		# short_self = [0, 0, ..., 0, -1, 0, ..., 0]
		short_self, conjugator = self.conjugate_short()
		short_lamination = conjugator(lamination)
		
		return sum(b for a, b in zip(short_self, short_lamination) if a == -1 and b >= 0)

