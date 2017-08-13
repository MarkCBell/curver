
''' A module for representing laminations on Triangulations.

Provides: Lamination, Shortenable. '''

from itertools import permutations

import curver
from curver.kernel.utilities import memoize  # Special import needed for decorating.

INFTY = float('inf')

def dual_weight(a, b, c):
	''' Return the dual weight of a triple.
	
	This is (b + c - a) / 2 when this is non-negative and b + c - a otherwise. '''
	
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
		# TODO: 2) Could save components if they have already been computed.
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
	
	def is_empty(self):
		''' Return if this lamination has no components. '''
		
		return not any(self)  # len(self) == 0
	
	def is_multicurve(self):
		''' Return if this lamination is actually a multicurve. '''
		
		return not self.is_empty() and all(isinstance(component, curver.kernel.Curve) for component in self.components())
	
	def is_curve(self):
		''' Return if this lamination is actually a curve. '''
		
		return self.is_multicurve() and len(self) == 1
	
	def is_multiarc(self):
		''' Return if this lamination is actually a multiarc. '''
		
		return not self.is_empty() and all(isinstance(component, curver.kernel.Arc) for component in self.components())
	
	def is_arc(self):
		''' Return if this lamination is actually a multiarc. '''
		
		return self.is_multiarc() and len(self) == 1
	
	def promote(self):
		''' Return this lamination in its finest form. '''
		
		# Could move cache across.
		if self.is_multicurve():
			if self.is_curve():
				return curver.kernel.Curve(self.triangulation, self.geometric)
			else:
				return curver.kernel.MultiCurve(self.triangulation, self.geometric)
		elif self.is_multiarc():
			if self.is_arc():
				return curver.kernel.Arc(self.triangulation, self.geometric)
			else:
				return curver.kernel.MultiArc(self.triangulation, self.geometric)
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
		
		return sum(multiplicity * component.intersection(lamination) for component, multiplicity in self.mcomponents())
	
	def no_common_component(self, lamination):
		''' Return that self does not share any components with the given Lamination. '''
		
		assert(isinstance(lamination, Lamination))
		
		self_components = self.components()
		return not any(component in self_components for component, _ in lamination.components())
	
	def train_track(self):
		''' Return the train track underlying this lamination. '''
		# In each triangle where this lamination looks like:
		# We introduce new edges to subdivide a triangle (p, q, r) as follows:
		#           / \                       /^\
		#          /   \                     / | \
		#         /     \                   /  |  \
		#        /-------\                 /   |s(i)
		#       /         \     ===>>     /    |    \
		#      /\         /\           r /    / \    \ q
		#     /  \       /  \           /   /     \   \
		#    /    |     |    \         /  /t(j) u(k)\  \
		#   /     |     |     \       /</             \>\
		#  /-------------------\     /-------------------\
		#                                      p
		# So that afterwards every complementary region can reach a vertex.
		
		geometric = list(self.geometric)
		triangles = []
		num_subdivided = 0  # Number of subdivided triangles.
		for triangle in self.triangulation:
			dual_weights = [self.dual_weight(label) for label in triangle.labels]
			if all(weight > 0 for weight in dual_weights):  # Type 3).
				p, q, r = [curver.kernel.Edge(label) for label in triangle.labels]
				s, t, u = [curver.kernel.Edge(i) for i in range(self.zeta + 3*num_subdivided, self.zeta + 3*num_subdivided + 3)]
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
			if isinstance(component, curver.kernel.Curve):
				components.append((curver.kernel.Curve(self.triangulation, component.geometric[:self.zeta]), multiplicity))
			elif isinstance(component, curver.kernel.Arc):
				components.append((curver.kernel.Arc(self.triangulation, component.geometric[:self.zeta]), multiplicity))
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
		return sum([multiplicity * component for component, multiplicity in self.mcomponents() if isinstance(component, curver.kernel.Arc)], empty)
	
	def multicurve(self):
		''' Return the maximal MultiCurve contained within this lamination. '''
		
		empty = self.triangulation.empty_lamination()
		return sum([multiplicity * component for component, multiplicity in self.mcomponents() if isinstance(component, curver.kernel.Curve)], empty)
	
	def boundary(self):
		''' Return the boundary of a regular neighbourhood of this lamination. '''
		
		multiarc = self.multiarc()  # Might be empty.
		multicurve = self.multicurve()  # Might be empty.
		return multicurve + (multiarc if multiarc.is_empty() else multiarc.boundary())

class Shortenable(Lamination):
	''' A special lamination that we can put into a canonical form. '''
	def is_short(self):
		''' Return if this lamination is short, that is, if it has as little weight as possible. '''
		
		return NotImplemented
	
	def shorten_score(self, edge):
		''' Return a score describing how good flipping this edge is for reducing the weight of this lamination.
		
		The higher the score, the better this flip is for reducing weight. '''
		
		return NotImplemented
	
	def shorten(self):
		''' Return an encoding which maps this lamination to a short one, together with its image. '''
		
		# TODO: 2) Make polynomial-time by taking advantage of spiralling.
		
		lamination = self
		conjugator = lamination.triangulation.id_encoding()
		
		extra = []
		while not lamination.is_short():
			edge = max(extra + lamination.triangulation.edges, key=lamination.shorten_score)
			# This edge is always flippable.
			
			move = lamination.triangulation.encode_flip(edge)
			conjugator = move * conjugator
			lamination = move(lamination)
			
			# TODO: 2) Accelerate!!
			a, b, c, d, e = lamination.triangulation.square(edge)
			extra = [a, d]
		
		return lamination, conjugator

