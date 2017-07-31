
''' A module for representing basic ways of changing triangulations.

Provides three classes: Isometry, EdgeFlip and LinearTransformation.

Perhaps in the future we will add a Spiral move so that curves can be
shortened in polynomial time. '''

import curver

class Move(object):
	''' Implements closedleaf and openleaf which apply this move to ClosedLeaf and OpenLeaf respectively. '''
	def __repr__(self):
		return str(self)
	def __invert__(self):
		return self.inverse()
	def __call__(self, other):
		if isinstance(other, curver.kernel.Lamination):
			return other.__class__(self.target_triangulation, {self(component): multiplicity for multiplicity, component in other})
		elif isinstance(other, curver.kernel.ClosedLeaf):
			return self.closedleaf(other)
		elif isinstance(other, curver.kernel.OpenLeaf):
			return self.openleaf(other)
		else:
			return NotImplemented
	
	def encode(self):
		''' Return the Encoding induced by this move. '''
		
		return curver.kernel.Encoding([self])

class Isometry(Move):
	''' This represents an isometry from one Triangulation to another.
	
	Triangulations can create the isometries between themselves and this
	is the standard way users are expected to create these. '''
	def __init__(self, source_triangulation, target_triangulation, label_map):
		''' This represents an isometry from source_triangulation to target_triangulation.
		
		It is given by a map taking each edge label of source_triangulation to a label of target_triangulation. '''
		
		assert(isinstance(source_triangulation, curver.kernel.Triangulation))
		assert(isinstance(target_triangulation, curver.kernel.Triangulation))
		assert(isinstance(label_map, dict))
		
		self.source_triangulation = source_triangulation
		self.target_triangulation = target_triangulation
		self.zeta = self.source_triangulation.zeta
		
		self.label_map = dict(label_map)
		
		# If we are missing any labels then use a depth first search to find the missing ones.
		# Hmmm, should always we do this just to check consistency?
		for i in self.source_triangulation.labels:
			if i not in self.label_map:
				raise curver.AssumptionError('This label_map not defined on edge %d.' % i)
		
		self.index_map = dict((i, curver.kernel.norm(self.label_map[i])) for i in self.source_triangulation.indices)
		# Store the inverses too while we're at it.
		self.inverse_label_map = dict((self.label_map[label], label) for label in self.source_triangulation.labels)
		self.inverse_index_map = dict((index, curver.kernel.norm(self.inverse_label_map[index])) for index in self.source_triangulation.indices)
	
	def __str__(self):
		return 'Isometry ' + str([self.target_triangulation.edge_lookup[self.label_map[i]] for i in self.source_triangulation.indices])
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.label_map))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		if not all(self.label_map[i] == i for i in self.source_triangulation.indices):  # If self is not the identity isometry.
			return {i: self.label_map[i] for i in self.source_triangulation.labels}
		else:
			return None
	
	def closedleaf(self, leaf):
		geometric = [leaf(self.inverse_index_map[index]) for index in self.source_triangulation.indices]
		return curver.kernel.ClosedLeaf(self.target_triangulation, geometric)
	
	def openleaf(self, leaf):
		geometric = [leaf(self.inverse_index_map[index]) for index in self.source_triangulation.indices]
		return curver.kernel.OpenLeaf(self.target_triangulation, geometric)
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		return Isometry(self.target_triangulation, self.source_triangulation, self.inverse_label_map)

class EdgeFlip(Move):
	''' Represents the change to a curve caused by flipping an edge. '''
	def __init__(self, source_triangulation, target_triangulation, edge_label):
		assert(isinstance(source_triangulation, curver.kernel.Triangulation))
		assert(isinstance(target_triangulation, curver.kernel.Triangulation))
		
		self.source_triangulation = source_triangulation
		self.target_triangulation = target_triangulation
		self.edge_label = edge_label
		self.edge_index = curver.kernel.norm(self.edge_label)
		self.zeta = self.source_triangulation.zeta
		assert(self.source_triangulation.is_flippable(self.edge_index))
		
		self.square = self.source_triangulation.square_about_edge(self.edge_label)
	
	def __str__(self):
		return 'Flip %s%d' % ('' if self.edge_index == self.edge_label else '~', self.edge_index)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.edge_label))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return self.edge_label
	
	def closedleaf(self, leaf):
		a, b, c, d, e = self.square
		
		# Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
		geometric = list(leaf.geometric)
		geometric[self.edge_index] = max(leaf(a) + leaf(c), leaf(b) + leaf(d)) - leaf(e)
		
		return curver.kernel.ClosedLeaf(self.target_triangulation, geometric)
	
	def openleaf(self, leaf):
		return NotImplemented
		a, b, c, d, e = self.square
		
		# Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
		geometric = list(leaf.geometric)
		
		if leaf.parallel() is not None:
			geometric[self.edge_index] = 1 if leaf.parallel() == self.edge_index else 0
		else:  # leaf is not parallel to an edge.
			if leaf(e) >= leaf(a) + leaf(b) and leaf(a) >= leaf(d) and leaf(b) >= leaf(c):  # CASE: A(ab)
				geometric[self.edge_index] = leaf(a) + leaf(b) - leaf(e)
			elif leaf(e) >= leaf(c) + leaf(d) and leaf(d) >= leaf(a) and leaf(c) >= leaf(b):  # CASE: A(cd)
				geometric[self.edge_index] = leaf(c) + leaf(d) - leaf(e)
			#elif leaf(e) <= 0 and leaf(a) >= leaf(b) and leaf(d) >= leaf(c):  # CASE: D(ad)
			#	geometric[self.edge_index] = leaf(a) + leaf(d) - leaf(e)
			#elif leaf(e) <= 0 and leaf(b) >= leaf(a) and leaf(c) >= leaf(d):  # CASE: D(bc)
			#	geometric[self.edge_index] = leaf(b) + leaf(c) - leaf(e) 
			elif leaf(a) >= leaf(b) + leaf(e) and leaf(d) >= leaf(c) + leaf(e):  # CASE: N(ad)
				geometric[self.edge_index] = leaf(a) + leaf(d) - 2*leaf(e)
			elif leaf(b) >= leaf(a) + leaf(e) and leaf(c) >= leaf(d) + leaf(e):  # CASE: N(bc)
				geometric[self.edge_index] = leaf(b) + leaf(c) - 2*leaf(e)
			elif leaf(a) + leaf(b) >= leaf(e) and leaf(b) + leaf(e) >= 2*leaf(c) + leaf(a) and leaf(a) + leaf(e) >= 2*leaf(d) + leaf(b):  # CASE: N(ab)
				geometric[self.edge_index] = (leaf(a) + leaf(b) - leaf(e)) // 2
			elif leaf(c) + leaf(d) >= leaf(e) and leaf(d) + leaf(e) >= 2*leaf(a) + leaf(c) and leaf(c) + leaf(e) >= 2*leaf(b) + leaf(d):  # CASE: N(cd)
				geometric[self.edge_index] = (leaf(c) + leaf(d) - leaf(e)) // 2
			else:
				geometric[self.edge_index] = max(leaf(a) + leaf(c), leaf(b) + leaf(d)) - leaf(e)
		
		return curver.kernel.OpenLeaf(self.target_triangulation, geometric)
	
	def inverse(self):
		''' Return the inverse of this map. '''
		
		return EdgeFlip(self.target_triangulation, self.source_triangulation, ~self.edge_label)

class Spiral(Move):
	''' This represents a spiral around a short curve. '''
	def __init__(self, source_triangulation, target_triangulation, edge_label, power):
		''' This represents spiralling around a short curve passing through edge_label.
		
		The number of spirals is determined by power and the compact form of this move
		means that the amount of work to compute the image of a curve under this move
		is logorithmic in the power.
		
		Because this is a mapping class, source_triangulation and target_triangulation should be equal. '''
		
		assert(isinstance(source_triangulation, curver.kernel.Triangulation))
		assert(isinstance(target_triangulation, curver.kernel.Triangulation))
		assert(source_triangulation == target_triangulation)
		
		self.source_triangulation = source_triangulation
		self.target_triangulation = target_triangulation
		self.zeta = self.source_triangulation.zeta
		
		self.edge_label = edge_label
		self.edge_index = curver.kernel.norm(self.edge_label)
		# Find a, b, c & d automatically.
		a, b, c, d, e = self.square
		assert(b == ~d)
		# Assert that b == d.
		self.power = power
	
	def __str__(self):
		return 'Spiral^%d %s' % (self.power, self.edge_label)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.edge_label, self.power))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return (self.edge_label, self.power)
	
	def mat(self, config, t):
		k = abs(self.power)
		if k == 0: return lambda (a, b, c, d): (a, b, c, d)
		
		def F(n):  # Note F(0) == Id.
			return lambda (a, b, c, d): (a, b, (n+1)*c - n*d, n*c + (1-n)*d)
		G = lambda (a, b, c, d): (a, b, a+b-d, c)
		
		# Take k steps in the graph. Leave unstable after t steps.
		if config == 1:  # Stable.
			M = lambda x: F(k)(x)
		elif config == 2:  # Transition.
			M = lambda x: F(k-1)(G(x))
		elif config == 3:  # Unstable.
			if k <= t:  # k steps in unstable
				M = lambda x: F(k)(x)
			elif k == t+1:  # t steps in unstable, then into transition.
				M = lambda x: G(F(t)(x))
			elif k == t+2:  # t steps in unstable, into transition, then out of transition.
				M = lambda x: G(G(F(t)(x)))
			else:  # k > t+2:  # t steps in unstable, into transition, out of transition, then the remaining steps in stable.
				M = lambda x: F(k - (t + 2))(G(G(F(t)(x))))
		
		return M
	
	def closedleaf(self, leaf):
		# We will begin with an easy case so we can later assume self.power != 0.
		
		a, b, c, d, e = self.square
		ai, bi, ci, di, ei = [self(edge) for edge in self.square]
		
		# Determine the number of strands passing through the annulus.
		xi = max(bi - ei, ai + ci - bi - ei, ei - bi)
		# Use that to determine the configuration we're in.
		# There are three possible states:
		#  1: Stable
		#  2: Transitioning
		#  3: Unstable
		# We use this slightly unusual ordering to ensure this gives the
		# same preference (a + c >= b + d) as EdgeFlip.
		state = 2 if xi == ai+ci-bi-ei else 1 if xi == bi-ei else 3
		
		k = abs(self.power)
		
		# Compute action on a, b, c, e.
		# Note that if self.power < 0 then (instead of using the ordering a, c, b, e) we use a, c, e, d.
		# This allows us to do two calculations with only one matrix.
		# Additionally d == b so we dont need to compute new_d.
		
		if self.power < 0:
			# Reverse the configuration if we are taking a negative power.
			state = 4 - state
			bi, ei = ei, bi
		
		# The maximum number of times you can perform the unstable state.
		t = min(max((2*bi - ai - ci) // (2*(ei - bi)) + 1, 0), k) if ei != bi else k
		M = self.mat(state, t)  # t only really matters if in the unstable, that is, if state == 3.
		_, _, new_bi, new_ei = M([ai, ci, bi, ei])
		
		if self.power < 0:
			# Reverse the configuration if we are taking a negative power.
			new_bi, new_ei = new_ei, new_bi
		
		geometric = [new_bi if index == b.index else new_e if index == e.index else leaf(index) for index in self.triangle.indices]
		return curver.kernel.ClosedLeaf(self.target_triangulation, geometric)
	
	def openleaf(self, leaf):
		return NotImplemented
		
		return curver.kernel.OpenLeaf(self.target_triangulation, geometric)
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		# inverse_corner_map = dict((self(corner), corner) for corner in self.corner_map)
		return Spiral(self.target_triangulation, self.source_triangulation, self.edge_label, -self.power)

