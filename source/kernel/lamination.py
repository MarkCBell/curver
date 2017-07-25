
''' A module for representing laminations on Triangulations.

Provides one class: Curve. '''

import curver

INFTY = float('inf')

class Lamination(object):
	''' This represents a lamination on an triangulation.
	
	Users can use Triangulation.lamination() as an alias for Lamination.from_weights(triangulation, ...). '''
	def __init__(self, triangulation, components):
		assert(isinstance(components, dict))
		
		self.triangulation = triangulation
		self.zeta = self.triangulation.zeta
		self.components = components
		# Could sanity check:
		# assert(all(c1.intersection(c2) == 0 for c1 in self for c2 in self))
		self.geometric = [sum(multiplicty * component(index) for component, multiplicity in self) for index in triangulation.indices]
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(self.geometric)
	def __iter__(self):
		return iter(self.components.items())
	def __len__(self):
		return sum(self.components.values())  # Total number of components.
	def __call__(self, item):
		''' Return the geometric measure assigned to item. '''
		# Could use self.geometric.
		return sum(multiplicity * component(item) for component, multiplicity in self)
	def __eq__(self, other):
		return self.triangulation == other.triangulation and self.components == other.components
	def __ne__(self, other):
		return not (self == other)
	def __hash__(self):
		return hash(self.geometric)
	
	def weight(self):
		''' Return the geometric intersection of this leaf with its underlying triangulation. '''
		
		return sum(multiplicity * component.weight() for component, multiplicity in self)
	
	def shorten(self):
		''' Return an encoding obtained by shortening each component in turn together with the image of self. '''
		
		conjugator = self.triangulation.id_encoding()
		for component, _ in self:
			conj, _ = conjugator(component).shorten()
			conjugator = conj * conjugator
		
		return conjugator, conjugator(self)
	
	@classmethod
	def from_weights(cls, triangulation, weights):
		''' Return the Lamination determined by the specified weights.
		
		Note: We currently install arbitary orientations on each component. '''
		
		return NotImplemented
		
		components = dict()
		
		# This is not really a Leaf but that has the right methods to make our lives easy.
		leaves = curver.kernel.Leaf(triangulation, weights, dict())
		
		# TODO: Modify geometric to remove any peripheral components before starting.
		# Until this is done MultiArc.boundary() wont work correctly.
		
		# max(min(leaves.dual_weight(leaves.triangulation.corner_lookup[edge.label].indices[1]) for edge in vertex_class), 0) for 
		# peripheral = dict((vertex, min(dual_weight(corner) for corner in self.corner_class_of_vertex(vertex))) for vertex in self.vertices)
		# geometric = [geometric[i] - sum(peripheral[v] for v in self.vertices_of_edge(i)) // 2 for i in range(self.zeta)]
		
		# Remove all the obvious arcs. This reduces the number of cases we have to consider later.
		for index in leaves.triangulation.indices:
			if leaves(index) < 0:
				geometric = [0 if i != index else -1 for i in range(self.zeta)]
				component, multiplicity = curver.kernel.OpenLeaf(lamination.triangulation, geometric, orientations).with_orientation(), abs(lamination(index))
				components[arc] = multiplicity
				lamination = lamination - multiplicity * arc
		
		# Remove all the obvious curves. This reduces the number of places we have to look for new curves later.
		for index in lamination.triangulation.indices:
			if lamination.triangulation.is_flippable(index):
				a, b, c, d, e = lamination.triangulation.square_about_edge(index)
				if b == ~d:
					multiplicity = ??
					if multiplicity > 0:
						curve = curver.kernel.ClosedLeaf(lamination.triangulation, geometric).with_orientation()
						components[curve] = multiplicity
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
		# WLOG: If there is a terminal normal arc in this triangle then it goes through side p.
		
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
				pass  # !?!
		
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
				component, multiplicity = curver.kernel.OpenLeaf(lamination.triangulation, geometric).with_orientation(), abs(lamination(to_flip))
				components[project(encoding.inverse()(component))] = multiplicity
				lamination = lamination - multiplicity * arc
			
			# The only places where a short curve could appear is across an edge adjacent to the one we just flipped.
			for edge in lamination.triangulation.square_about_edge(to_flip):
				if lamination.triangulation.is_flippable(edge.label):
					a, b, c, d, e = lamination.triangulation.square_about_edge(edge.label)
					if b == ~d:
						multiplicity = ??
						if multiplicity > 0:
							component = curver.kernel.ClosedLeaf(lamination.triangulation, geometric).with_orientation()
							components[project(encoding.inverse()(component))] = multiplicity
							lamination = lamination - multiplicity * curve
			
			# Accelerate!!
			
			to_flip = ??
		
		return Lamination(triangulation, components).promote()
	
	def is_multicurve(self):
		return all(isinstance(component, curver.kernel.ClosedLeaf) for component, _ in self)
	
	def is_curve(self):
		return self.is_multicurve() and len(self) == 1
	
	def is_multiarc(self):
		return all(isinstance(component, curver.kernel.OpenLeaf) for component, _ in self)
	
	def is_arc(self):
		return self.is_multiarc() and len(self) == 1
	
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
	
	def is_empty(self):
		''' Return if this lamination has no components. '''
		
		return len(self) == 0
	
	def skeleton(self):
		''' Return the lamination obtained by collapsing parallel components. '''
		
		return Lamination(self.triangulation, {component: 1 for component, _ in self})
	
	def intersection(self, lamination):
		''' Return the geometric intersection number between this lamination and the given one. '''
		
		assert(isinstance(lamination, Lamination))
		
		return sum(m1 * m2 * c1.intersection(c2) for c1, m1 in self for c2, m2 in lamination)

class MultiCurve(Lamination):
	''' A Lamination in which every component is a ClosedLeaf. '''
	def is_multicurve(self):
		return True
	def is_multiarc(self):
		return False

class Curve(Multicurve):
	''' A MultiCurve with a single component. '''
	
	def is_isolating(self):
		''' Return if this curve is an isolating curve.
		
		That is, is there a component of S - self that does not contain a puncture. '''
		
		# This is based off of self.encode_twist(). See the documentation there as to why this works.
		
		short, _ = self.shorten()
		
		return short.weight() == 2
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve raised to the power k.
		
		Currently, this must be a non-isolating curve. '''
		
		assert(not self.is_isolating())
		
		if k == 0: return self.triangulation.id_encoding()
		
		short, conjugation = self.shorten()
		
		triangulation = short.triangulation
		# Grab the indices of the two edges we meet.
		e1, e2 = [edge_index for edge_index in short.triangulation.indices if short(edge_index) > 0]
		
		a, b, c, d, e = triangulation.square_about_edge(e1)
		# If the curve is going vertically through the square then ...
		if short(a) == 1 and short(c) == 1:
			# swap the labels round so it goes horizontally.
			e1, e2 = e2, e1
			a, b, c, d, e = triangulation.square_about_edge(e1)
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
		
		twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [e1, e2]}, e1]))
		return conjugation.inverse() * twist**k * conjugation
		
		# Once Spiral is working we can do:
		twist_k = triangulation.encode([(e1, k)])
		return conjugation.inverse() * twist_k * conjugation
	
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
	''' A Lamination in which every component is an OpenLeaf. '''
	def is_multicurve(self):
		return False
	def is_multiarc(self):
		return True
	
	def boundary(self):
		''' Return the multicurve which is the boundary of a regular neighbourhood of this multiarc. '''
		
		# TODO: Should build correct orientations on boundary components.
		
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
		
		boundary = short.triangulation.lamination(geometric).promote()
		
		return conjugator.inverse()(boundary)

class Arc(MultiArc):
	''' A MultiArc with a single component. '''
	
	def encode_halftwist(self, k=1):
		''' Return an Encoding of a left half twist about a regular neighbourhood of this arc raised to the power k.
		
		Assumes (and checks) that this arc connects between distinct vertices.
		Currently, the boundary curve must also be non-isolating. '''
		
		boundary = self.boundary()
		if not isinstance(boundary, Curve):  # isinstance(boundary, MultiCurve):
			raise curver.AssumptionError('Arc connects a vertex to itself.')
		
		if boundary.is_empty():  # Surface == S_{0, 3}:
			return  # !?!
		
		if k % 2 == 0:  # k is even so use a Dehn twist about the boundary.
			return boundary.encode_twist(k // 2)
		
		short_boundary, conjugation = boundary.shorten()
		short = conjugation(self)
		
		if not short_boundary.weight() == 2:  # not boundary.is_isolating().
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
		_, tilde_e, _, _, _ = triangulation.square_about_edge(c.label)
		a, b, c, d, e = triangulation.square_about_edge(~tilde_e.label)
		
		half_twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [b.index, e.index, c.index, x.index]}, b.index, e.index, c.index])
		
		# We accelerate large powers by replacing (T^1/2_self)^2 with T_self which includes acceleration.
		if abs(k) == 1:
			return conjugation.inverse() * half_twist**k * conjugation
		else:  # k is odd so we need to add in an additional half twist.
			# Note: k // 2 always rounds down, so even if k < 0 the additional half twist we need to do is positive.
			return conjugation.inverse() * short_boundary.encode_twist(k // 2) * half_twist * conjugation
	

