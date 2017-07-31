
def dual_weight(a, b, c):
	correction = min(a + b - c, b + c - a, c + a - b)
	return (a + b - c + min(correction, 0)) // 2

class Leaf(object):
	''' A 1-d object drawn on a Triangulation and determined by its number of intersections with each Edge. '''
	def __init__(self, triangulation, geometric):
		assert(isinstance(triangulation, curver.kernel.Triangulation))
		assert(isinstance(geometric, (list, tuple)))  # Maps indices to measures.
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
	def __hash__(self):
		return hash(tuple(self.geometric))
	def __eq__(self, other):
		return self.triangulation == other.triangulation and \
			self.geometric == other.geometric
	def __ne__(self, other):
		return not (self == other)
	def __call__(self, other):
		''' Return the geometric measure assigned to item. '''
		if isinstance(other, curver.kernel.Edge): other = other.label
		
		return self.geometric[curver.kernel.norm(other)]
	
	def weight(self):
		''' Return the geometric intersection of this leaf with its underlying triangulation. '''
		
		# Negative weights correspond to leaves that run parallel to an edge so don't contribute to the geometric intersection number.
		return sum(max(x, 0) for x in self)
	
	def dual_weight(self, label):
		''' Return the number of component of this leaf dual to the given edge.
		
		Note that when there is a terminal normal arc then we record this weight with a negative sign. '''
		
		corner = self.trinagulation.corner_lookup[label]
		weights = [self(edge) for edge in corner]
		return dual_weight(*weights)
	
	def arc_side(self, triangle):
		dual_weights = self.dual_weights(triangle)
		arc = min(dual_weights)
		if arc >= 0:  # No arc.
			return None
		else:
			return dual_weights.index(arc)


class ClosedLeaf(Leaf):
	def shorten(self):
		''' Return an encoding which maps this leaf to a leaf with as little weight as possible together with its image. '''
		
		# Repeatedly flip to reduce the weight of this leaf as much as possible.
		# TODO: 3) Make polynomial-time by taking advantage of spiralling.
		
		leaf = self
		conjugation = leaf.triangulation.id_encoding()
		
		weight_history = [INFTY, INFTY, leaf.weight()]
		# If we ever fail to make progress more than once then the leaf is as short as it's going to get.
		while weight_history[-1] < weight_history[-3]:
			# Find the flip which decreases our weight the most.
			flips = [leaf.triangulation.encode_flip(index) for index in leaf.triangulation.indices if leaf.triangulation.is_flippable(index)]
			flip = min(flips, key=lambda flip: flip(leaf).weight())
			
			conjugation = flip * conjugation
			leaf = flip(leaf)
			weight_history.append(leaf.weight())
		
		return leaf, conjugation
	
	def intersection(self, leaf):
		''' Return the geometric intersection between self and the given lamination.
		
		Currently assumes (and checks) that self is a non-isolating curve. '''
		
		assert(isinstance(leaf, Leaf))
		if isinstance(leaf, OpenLeaf):
			return leaf.intersection(self)
		
		assert(leaf.triangulation == self.triangulation)
		
		short, conjugator = self.shorten()
		short_leaf = conjugator(leaf)
		
		if short.weight() == 2:
			triangulation = short.triangulation
			e1, e2 = [index for index in triangulation.indices if short(index) > 0]
			# We might need to swap these edge indices so we have a good frame of reference.
			if triangulation.corner_lookup[e1].indices[2] != e2: e1, e2 = e2, e1
			
			a, b, c, d, e = triangulation.square_about_edge(e1)
			
			x = (short_leaf(a) + short_leaf(b) - short_leaf(e)) // 2
			y = (short_leaf(b) + short_leaf(e) - short_leaf(a)) // 2
			z = (short_leaf(e) + short_leaf(a) - short_leaf(b)) // 2
			x2 = (short_leaf(c) + short_leaf(d) - short_leaf(e)) // 2
			y2 = (short_leaf(d) + short_leaf(e) - short_leaf(c)) // 2
			z2 = (short_leaf(e) + short_leaf(c) - short_leaf(d)) // 2
			
			return short_leaf(a) - 2 * min(x, y2, z)  # = short_leaf(c) - 2 * min(x2, y, z2))
		else:
			# TODO: 4) Implement LP to find intersection for general configuration.
			raise curver.AssumptionError('Currently can only compute geometric intersection number between a non-isolating ClosedLeaf and a Leaf.')

class OpenLeaf(Leaf):
	def parallel(self):
		''' Return the Edge that this Leaf is parallel to or None if it is not parallel to any Edge. '''
		if all(x >= 0 for x in self):  # Leaf is not parallel to an Edge.
			return None
		
		edge_index = self.geometric.index(-1)
		if self[edge_index] == +1:
			edge_label = edge_index
		else:  # self[edge_index] == -1:
			edge_label = ~edge_index
		
		return self.triangulation.edge_lookup[edge_label]
		
	def shorten(self):
		''' Return an encoding which maps this leaf to a leaf with as little weight as possible together with its image.
		
		Uses Mosher's arguement. '''
		
		# TODO: 3) Make polynomial-time by taking advantage of spiralling.
		
		leaf = self
		conjugation = leaf.triangulation.id_encoding()
		
		while leaf.weight() > 0:
			labels = [label for label in leaf.triangulation.labels if leaf.dual_weight(label) < 0]
			label = labels[0]
			flip = leaf.triangulation.encode_flip(label)
			
			conjugation = flip * conjugation
			leaf = flip(leaf)
		
		return leaf, conjugation
	
	def intersection(self, leaf):
		''' Return the geometric intersection between self and the given lamination. '''
		
		assert(isinstance(leaf, Leaf))
		assert(leaf.triangulation == self.triangulation)
		
		# short = [0, 0, ..., 0, -1, 0, ..., 0]
		short, conjugator = self.shorten()
		short_leaf = conjugator(leaf)
		
		return sum(b for a, b in zip(short, short_leaf) if a == -1 and b >= 0)

