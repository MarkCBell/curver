
''' A module for representing laminations on Triangulations.

Provides one class: Curve. '''

import networkx
from itertools import permutations

import curver

INFTY = float('inf')

class Lamination(object):
	''' This represents a lamination on an triangulation.
	
	Users can use Triangulation.lamination(). '''
	def __init__(self, triangulation, components):
		assert(isinstance(components, dict))
		
		self.triangulation = triangulation
		self.zeta = self.triangulation.zeta
		self.components = components
		# Could sanity check:
		# assert(all(c1.intersection(c2) == 0 for c1 in self for c2 in self))
		self.geometric = tuple(sum(multiplicty * component(index) for component, multiplicity in self) for index in triangulation.indices)
	
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
	
	def is_empty(self):
		''' Return if this lamination has no components. '''
		
		return len(self) == 0
	
	def is_multicurve(self):
		return not self.is_empty() and all(isinstance(component, curver.kernel.ClosedLeaf) for component, _ in self)
	
	def is_curve(self):
		return self.is_multicurve() and len(self) == 1
	
	def is_multiarc(self):
		return not self.is_empty() and all(isinstance(component, curver.kernel.OpenLeaf) for component, _ in self)
	
	def is_arc(self):
		return self.is_multiarc() and len(self) == 1
	
	def promote(self):
		elif self.is_multicurve():
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
		
		return Lamination(self.triangulation, {component: 1 for component, _ in self}).promote()
	
	def peek_component(self):
		''' Return one component of this Lamination. '''
		
		return Lamination(self.triangulation, {self.components.keys()[0]: 1}).promote()
	
	def intersection(self, lamination):
		''' Return the geometric intersection number between this lamination and the given one. '''
		
		assert(isinstance(lamination, Lamination))
		
		return sum(m1 * m2 * c1.intersection(c2) for c1, m1 in self for c2, m2 in lamination)
	
	def is_disjoint(self, lamination):
		''' Return that self does not share any components with the given Lamination. '''
		assert(isinstance(lamination, Lamination))
		
		return not any(component in self.components for component, _ in other)
	
	def arcs_and_curves(self):
		''' Return the set of Arcs and Curves that appear within self. '''
		
		return set(Lamination(self.triangulation, {component: 1}).promote() for component, _ in self)
	
	def sublaminations(self):
		''' Return all sublaminations that appear within self. '''
		return [Lamination(self.triangulation, dict(sub)).promote() for i in range(len(self.components)) for sub in permutations(self, i)]

class MultiCurve(Lamination):
	''' A Lamination in which every component is a ClosedLeaf. '''
	def is_multicurve(self):
		return True
	def is_multiarc(self):
		return False
	
	def encode_twist(self, k=1):
		h = self.triangulation.id_encoding()
		for curve, multiplicity in self:
			h = Curve(self.triangulation, {curve: 1}).encode_twist(k * multiplicity) * h
		
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
		
		for curve in self.arcs_and_curves():
			h = g(curve).crush()  # Map forward under crushes first.
			g = h * g
		
		return g

class Curve(Multicurve):
	''' A MultiCurve with a single component. '''
	
	def encode_twist(self, k=1):
		''' Return an Encoding of a left Dehn twist about this curve raised to the power k.
		
		Currently, this must be a non-isolating curve. '''
		
		if k == 0: return self.triangulation.id_encoding()
		
		short, conjugation = self.shorten()
		
		if short.weight() == 2:  # curve is non-isolating.
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
			
			twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [e1, e2]}, e1])
			return conjugation.inverse() * twist**k * conjugation
			
			# Once Spiral is working we can do:
			twist_k = triangulation.encode([(e1, k)])
			return conjugation.inverse() * twist_k * conjugation
		else:  # curve is isolating.
			raise curver.AssumptionError('Curve is isolating.')  # TODO: Handle isolating case.
	
	def quasiconvex(self, other):
		''' Return a polynomial-sized K--quasiconvex subset of the curve complex that contains self and other. '''
		
		assert(isinstance(other, Curve))
		
		# TODO: Implement train track splitting sequence.
		
		return NotImplemented
	
	def all_tight_geodesic_multicurves(self, other):
		''' Return a set that contains all multicurves in any tight geodesic from self to other.
		
		From the first half of Algorithm 4 of Paper 3. '''
		
		assert(isinstance(other, Curve))
		
		guide = self.quasiconvex(other)  # U.
		L = 6*curver.kernel.constants.QUASICONVEXITY + 2
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
		
		# TODO: Implement this and the associated Move.
		
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
		
		boundary = short.triangulation.lamination(geometric)
		
		return conjugator.inverse()(boundary)
	
	def explore_ball(self, radius):
		''' Extend this MultiArc to a triangulation and return all triangulations within the ball of the given radius of that one.
		
		Note that this is only well defined if self is filling. '''
		
		short, conjugator = self.shorten()
		
		X = set()
		for encoding in short.triangulation.all_encodings(radius):
			T = encoding.target_triangulation.as_lamination()
			X.add(conjugator.inverse()(encoding.inverse()(T)))
		
		return X

class Arc(MultiArc):
	''' A MultiArc with a single component. '''
	
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
		
		short_boundary, conjugation = boundary.shorten()
		short = conjugation(self)
		
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
			_, tilde_e, _, _, _ = triangulation.square_about_edge(c.label)
			a, b, c, d, e = triangulation.square_about_edge(~tilde_e.label)
			
			half_twist = triangulation.encode([{i: i for i in triangulation.indices if i not in [b.index, e.index, c.index, x.index]}, b.index, e.index, c.index])
			
			# We accelerate large powers by replacing (T^1/2_self)^2 with T_self which includes acceleration.
			if abs(k) == 1:
				return conjugation.inverse() * half_twist**k * conjugation
			else:  # k is odd so we need to add in an additional half twist.
				# Note: k // 2 always rounds down, so even if k < 0 the additional half twist we need to do is positive.
				return conjugation.inverse() * short_boundary.encode_twist(k // 2) * half_twist * conjugation
		else:  # boundary is isolating.
			raise curver.AssumptionError('Boundary curve is isolating.')  # TODO: Handle isolating case, use Will Worden's code.

