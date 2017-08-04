
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
			return Lamination(self.triangulation, geometric)
		elif other == 0:  # So we can use sum.
			return self
		else:
			return NotImplemented
	def __radd__(self, other):
		return self + other
	def __mul__(self, other):
		geometric = [other * x for x in self]
		return Lamination(self.triangulation, geometric)
	def __rmul__(self, other):
		return self * other
	def __sub__(self, other):
		# Haken sum.
		if isinstance(other, Lamination):
			if other.triangulation != self.triangulation:
				raise ValueError('Laminations must be on the same triangulation to add them.')
			
			geometric = [x - y for x, y in zip(self.geometric, other.geometric)]
			return Lamination(self.triangulation, geometric)
		else:
			return NotImplemented
	
	def weight(self):
		''' Return the geometric intersection of this leaf with its underlying triangulation. '''
		
		return sum(max(weight, 0) for weight in self)
	
	def dual_weight(self, edge):
		''' Return the number of component of this lamination dual to the given edge.
		
		Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		corner = self.triangulation.corner_lookup[edge.label]
		weights = [self(edge) for edge in corner]
		return dual_weight(*weights)
	
	def shorten(self):
		''' Return an encoding obtained by shortening each component in turn together with the image of self. '''
		
		conjugator = self.triangulation.id_encoding()
		for component in self.components():
			_, conj = conjugator(component).shorten()
			conjugator = conj * conjugator
		
		return conjugator(self), conjugator
	
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
	
	def skeleton(self):
		''' Return the lamination obtained by collapsing parallel components. '''
		
		return sum(component for component in self.components()).promote()
	
	def peek_component(self):
		''' Return one component of this Lamination. '''
		
		return self.components().pop()
	
	def intersection(self, lamination):
		''' Return the geometric intersection number between this lamination and the given one. '''
		
		assert(isinstance(lamination, Lamination))
		
		return sum(m1 * m2 * c1.intersection(c2) for c1, m1 in self.mcomponents() for c2, m2 in lamination.mcomponents())
	
	def is_disjoint(self, lamination):
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
		return [sum(components).promote() for i in range(len(components)) for sub in permutations(components, i)]

class MultiCurve(Lamination):
	''' A Lamination in which every component is a ClosedLeaf. '''
	def is_multicurve(self):
		return True
	def is_multiarc(self):
		return False
	
	def encode_twist(self, k=1):
		h = self.triangulation.id_encoding()
		for curve, multiplicity in self.mcomponents():
			h = curve.encode_twist(k * multiplicity) * h
		
		return h
	
	def boundary_neighbourhood_union(other):
		''' Return \partial N(self \cup other). '''
		assert(isinstance(other, MultiCurve))  # Could probably do any Lamination.
		
		crush = self.crush()
		lift = crush.inverse()
		other_prime = crush(other)
		m_prime = other_prime.boundary()
		return lift(m_prime)  # = m.
	
	def tight_paths(self, other, length):
		''' Return the set of all tight paths from self to other that are of the given length.
		
		From Algorithm 3 of Paper 3. '''
		assert(isinstance(other, MultiCurve))
		assert(length >= 0)
		
		if length == 0:
			return set([(self,)]) if self == other else set()
		elif length == 1:
			return set([(self, other)]) if self.intersection(other) == 0 and self.is_disjoint(other) else set()
		elif length == 2:
			m = self.boundary_neighbourhood_union(other)  # m = \partial N(self \cup other).
			return set([(self, m, other)]) if not m.is_empty() and self.is_disjoint(m) and other.is_disjoint(m) else set()
		else:  # length >= 3.
			crush = self.crush()
			lift = crush.inverse()
			b_prime = crush(other)
			A_1 = set()
			for triangulation in b_prime.eplore_ball(2*self.zeta*length + 2*self.zeta):
				for submultiarc in triangulation.sublaminations():
					m_prime = submultiarc.boundary()
					m = lift(m_prime)
					A_1.add(m)
			
			P = set()
			for a_1 in A_1:
				for multipath in a_1.tight_paths(other, length-1):  # Recurse.
					a_2 = multipath[0]
					if self.boundary_neighbourhood_union(a_2) == a_1:  # (self,) + multipath is tight:
						P.add((self,) + multipath)
			
			return P
	
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
	
	def shorten(self):
		''' Return an encoding which maps this curve to a curve with as little weight as possible together with its image. '''
		
		# Repeatedly flip to reduce the weight of this curve as much as possible.
		# TODO: 3) Make polynomial-time by taking advantage of spiralling.
		
		curve = self
		conjugator = curve.triangulation.id_encoding()
		
		weight_history = [INFTY, INFTY, curve.weight()]
		# If we ever fail to make progress more than once then the curve is as short as it's going to get.
		while weight_history[-1] < weight_history[-3]:
			# Find the flip which decreases our weight the most.
			flips = [curve.triangulation.encode_flip(index) for index in curve.triangulation.indices if curve.triangulation.is_flippable(index)]
			flip = min(flips, key=lambda flip: flip(curve).weight())
			
			conjugator = flip * conjugator
			curve = flip(curve)
			weight_history.append(curve.weight())
		
		return curve, conjugator
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve raised to the power k.
		
		Currently, this must be a non-isolating curve. '''
		
		if k == 0: return self.triangulation.id_encoding()
		
		short, conjugator = self.shorten()
		
		if short.weight() == 2:  # curve is non-isolating.
			triangulation = short.triangulation
			# Grab the indices of the two edges we meet.
			e1, e2 = [edge_index for edge_index in short.triangulation.indices if short(edge_index) > 0]
			
			a, b, c, d, e = triangulation.square(e1)
			# If the curve is going vertically through the square then ...
			if short(a) == 1 and short(c) == 1:
				# swap the labels round so it goes horizontally.
				e1, e2 = e2, e1
				a, b, c, d, e = triangulation.square(e1)
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
			# And e.index = e1 and b == ~d
			
			twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [e1, e2]}, e1])
			return conjugator.inverse() * twist**k * conjugator
			
			# Once Spiral is working we can do:
			twist_k = triangulation.encode([(e1, k)])
			return conjugator.inverse() * twist_k * conjugator
		else:  # curve is isolating.
			raise curver.AssumptionError('Curve is isolating.')  # TODO: 4) Handle isolating case.
	
	def intersection(self, other):
		''' Return the geometric intersection between self and the given other.
		
		Currently assumes (and checks) that self is a non-isolating curve. '''
		
		assert(isinstance(other, (Arc, Curve)))
		if isinstance(other, Arc):
			return other.intersection(self)
		
		assert(other.triangulation == self.triangulation)
		
		short, conjugator = self.shorten()
		short_other = conjugator(other)
		
		if short.weight() == 2:
			triangulation = short.triangulation
			e1, e2 = [index for index in triangulation.indices if short(index) > 0]
			# We might need to swap these edge indices so we have a good frame of reference.
			if triangulation.corner_lookup[e1].indices[2] != e2: e1, e2 = e2, e1
			
			a, b, c, d, e = triangulation.square(e1)
			
			x = (short_other(a) + short_other(b) - short_other(e)) // 2
			y = (short_other(b) + short_other(e) - short_other(a)) // 2
			z = (short_other(e) + short_other(a) - short_other(b)) // 2
			x2 = (short_other(c) + short_other(d) - short_other(e)) // 2
			y2 = (short_other(d) + short_other(e) - short_other(c)) // 2
			z2 = (short_other(e) + short_other(c) - short_other(d)) // 2
			
			return short_other(a) - 2 * min(x, y2, z)  # = short_other(c) - 2 * min(x2, y, z2))
		else:
			# TODO: 4) Implement LP to find intersection for general configuration.
			raise curver.AssumptionError('Currently can only compute geometric intersection number between a non-isolating ClosedLeaf and a Leaf.')
	
	def quasiconvex(self, other):
		''' Return a polynomial-sized K--quasiconvex subset of the curve complex that contains self and other. '''
		
		assert(isinstance(other, Curve))
		short, conjugator = self.shorten()
		
		train_track = conjugator(other).train_track()
		_, conjugator_tt = train_track.shorten()
		encodings = [conjugator_tt[i:] for i in range(len(conjugator_tt)+1)]
		return [conjugator.inverse()(encoding.inverse()(encoding(train_track).vertex_cycle())) for encoding in encodings]
	
	def all_tight_geodesic_multicurves(self, other):
		''' Return a set that contains all multicurves in any tight geodesic from self to other.
		
		From the first half of Algorithm 4 of Paper 3. '''
		
		assert(isinstance(other, Curve))
		
		guide = self.quasiconvex(other)  # U.
		L = 6*curver.kernel.constants.QUASICONVEXITY + 2  # See Richard's paper!?!
		return set(multicurve for length in range(L+1) for c1 in guide for c2 in guide for path in c1.tight_paths(c2, length) for multicurve in path)
	
	def tight_geodesic(self, other):
		''' Return a tight geodesic in the (multi)curve complex from self to other.
		
		From the second half of Algorithm 4 of Paper 3. '''
		
		assert(isinstance(other, Curve))
		
		vertices = list(self.all_tight_geodesic_multicurves(other))
		edges = [(i, j) for i in range(len(vertices)) for j in range(i) if vertices[i].intersection(vertices[j]) == 0 and vertices[i].is_disjoint(vertices[j])]
		
		G = networkx.Graph(edges)  # Build graph.
		indices = networkx.shortest_path(G, vertices.index(self), vertices.index(other))  # Find a geodesic from self to other.
		geodesic = [vertices[index] for index in indices]  # Get the geodesic, however this might not be tight.
		
		for i in range(1, len(geodesic)-1):
			geodesic[i] = geodesic[i-1].boundary_neighbourhood_union(geodesic[i+1])  # Tighten.
		
		return tuple(geodesic)
	
	def geodesic(self, other):
		''' Return a geodesic in the curve complex from self to other.
		
		The geodesic will always come from a tight geodesic.
		From Algorithm 5 of Paper 3. '''
		
		return tuple(multicurve.peek_component() for multicurve in self.tight_geodesic(other))
	
	def distance(self, other):
		''' Return the distance from self to other in the curve complex. '''
		
		# Could use self.tight_geodesic(other).
		return len(self.geodesic(other)) - 1
	
	def crush(self):
		''' Return the crush map associated to this Curve. '''
		
		# TODO: 2) Implement this and the associated Move.
		
		return NotImplemented

class MultiArc(Lamination):
	''' A Lamination in which every component is an OpenLeaf. '''
	def is_multicurve(self):
		return False
	def is_multiarc(self):
		return True
	
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
		
		X = set()
		for encoding in short.triangulation.all_encodings(radius):
			T = encoding.target_triangulation.as_lamination()
			X.add(conjugator.inverse()(encoding.inverse()(T)))
		
		return X

class Arc(MultiArc):
	''' A MultiArc with a single component. '''
	def mcomponents(self):
		return [(self, 1)]
	
	def shorten(self):
		''' Return an encoding which maps this arc to a arc with as little weight as possible together with its image.
		
		Uses Mosher's arguement. '''
		
		# TODO: 3) Make polynomial-time by taking advantage of spiralling.
		
		arc = self
		conjugator = arc.triangulation.id_encoding()
		
		while arc.weight() > 0:
			labels = [label for label in arc.triangulation.labels if arc.dual_weight(label) < 0]
			label = labels[0]
			flip = arc.triangulation.encode_flip(label)
			
			conjugator = flip * conjugator
			arc = flip(arc)
		
		return arc, conjugator
	
	def intersection(self, other):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(other, Leaf))
		assert(other.triangulation == self.triangulation)
		
		# short = [0, 0, ..., 0, -1, 0, ..., 0]
		short, conjugator = self.shorten()
		short_other = conjugator(other)
		
		return sum(b for a, b in zip(short, short_other) if a == -1 and b >= 0)
	
	def encode_halftwist(self, k=1):
		''' Return an Encoding of a left half twist about a regular neighbourhood of this arc raised to the power k.
		
		Assumes (and checks) that this arc connects between distinct vertices.
		Currently, the boundary curve must also be non-isolating. '''
		
		boundary = self.boundary()
		assert(not boundary.is_empty())  # Surface == S_{0, 3}
		
		if not isinstance(boundary, Curve):  # isinstance(boundary, MultiCurve):
			raise curver.AssumptionError('Arc connects a vertex to itself.')
		
		if k % 2 == 0:  # k is even so use a Dehn twist about the boundary.
			return boundary.encode_twist(k // 2)
		
		short_boundary, conjugator = boundary.shorten()
		short = conjugator(self)
		
		if short_boundary.weight() == 2:  # boundary is non-isolating.
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
			_, tilde_e, _, _, _ = triangulation.square(c.label)
			a, b, c, d, e = triangulation.square(~tilde_e.label)
			
			half_twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [b.index, e.index, c.index, x.index]}, b.index, e.index, c.index])
			
			# We accelerate large powers by replacing (T^1/2_self)^2 with T_self which includes acceleration.
			if abs(k) == 1:
				return conjugator.inverse() * half_twist**k * conjugator
			else:  # k is odd so we need to add in an additional half twist.
				# Note: k // 2 always rounds down, so even if k < 0 the additional half twist we need to do is positive.
				return conjugator.inverse() * short_boundary.encode_twist(k // 2) * half_twist * conjugator
		else:  # boundary is isolating.
			raise curver.AssumptionError('Boundary curve is isolating.')  # TODO: 4) Handle isolating case, use Will Worden's code.

