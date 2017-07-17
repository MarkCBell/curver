
''' A module for representing basic ways of changing triangulations.

Provides three classes: Isometry, EdgeFlip and LinearTransformation.

Perhaps in the future we will add a Spiral move so that curves can be
shortened in polynomial time. '''

import curver

class Move(object):
	def __repr__(self):
		return str(self)
	
	def encode(self):
		''' Return the Encoding induced by this isometry. '''
		
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
		self.inverse_label_map = dict((self.label_map[i], i) for i in self.source_triangulation.labels)
		self.inverse_index_map = dict((i, curver.kernel.norm(self.inverse_label_map[i])) for i in self.source_triangulation.indices)
		self.inverse_signs = dict((i, +1 if self.inverse_index_map[i] == self.inverse_label_map[i] else -1) for i in self.source_triangulation.indices)
	
	def __str__(self):
		return 'Isometry ' + str([self.target_triangulation.edge_lookup[self.label_map[i]] for i in self.source_triangulation.indices])
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.label_map))
	def __len__(self):
		return 1  # The number of pieces of this move.
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		if not all(self.label_map[i] == i for i in self.source_triangulation.indices):  # If self is not the identity isometry.
			return {i: self.label_map[i] for i in self.source_triangulation.labels}
		else:
			return None
	
	def apply_geometric(self, vector):
		return [vector[self.inverse_index_map[i]] for i in range(self.zeta)]
	def apply_algebraic(self, vector):
		return [vector[self.inverse_index_map[i]] * self.inverse_signs[i] for i in range(self.zeta)]
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		# inverse_corner_map = dict((self(corner), corner) for corner in self.corner_map)
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
	def __len__(self):
		return 2  # The number of pieces of this move.
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return self.edge_label
	
	def apply_geometric(self, vector):
		a, b, c, d = self.square
		m = max(vector[a.index] + vector[c.index], vector[b.index] + vector[d.index]) - vector[self.edge_index]
		return [vector[i] if i != self.edge_index else m for i in range(self.zeta)]
	def apply_algebraic(self, vector):
		a, b, c, d = self.square
		m = b.sign() * vector[b.index] + c.sign() * vector[c.index]
		return [vector[i] if i != self.edge_index else m for i in range(self.zeta)]
	
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
		self.square = self.source_triangulation.square_about_edge(self.edge_label)
		a, b, c, d = self.square
		assert(b == ~d)
		# Assert that b == d.
		self.power = power
	
	def __str__(self):
		return 'Spiral^%d %s' % (self.power, self.edge_label)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.edge_label, self.power))
	def __len__(self):
		return abs(self.power) + 3  # The number of pieces of this move.
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return (self.edge_label, self.power)
	
	def mat(self, config, t):
		k = abs(self.power)
		if k == 0: return curver.kernel.id_matrix(4)
		
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
			else:  # t+2 < k:  # t steps in unstable, into transition, out of transition, then the remaining steps in stable.
				M = lambda x: F(k - (t + 2))(G(G(F(t)(x))))
		
		return M
	
	def apply_geometric(self, vector):
		# We will begin with an easy case so we can later assume self.power != 0.
		
		ai, bi, ci, di = [edge.index for edge in self.square]
		ei = self.edge_index
		a, b, c, d = [vector[edge.index] for edge in self.square]
		e = vector[self.edge_index]
		
		# Determine the number of strands passing through the annulus.
		x = max(b - e, a + c - b - e, e - b)
		# Use that to determine the configuration we're in.
		# There are three possible stats:
		#  1: Stable
		#  2: Transitioning
		#  3: Unstable
		# We use this slightly unusual ordering to ensure this gives the
		# same preference (a + c >= b + d) as EdgeFlip.
		state = 2 if x == a+c-b-e else 1 if x == b-e else 3
		
		k = abs(self.power)
		
		# Compute action on a, b, c, e.
		# Note that if self.power > 0 then we use the ordering a, c, b, e
		# and otherwise we use a, c, e, d. This allows us to do two calculations
		# with only one matrix.
		# Additionally d == b so we dont need to compute new_d.
		
		# The maximum number of times you can perform the unstable state.
		# WLOG 0 <= t <= k.
		
		if self.power > 0:
			config = state
			t = min(max((2*b - a - c) // (2*(e - b)) + 1, 0), k) if e != b else k
			M = self.mat(config, t)  # t only matters if in the unstable, that is, if config == 3.
			_, _, new_b, new_e = M([a, c, b, e])
		else:
			# Reverse the configuration if we are taking a negative power.
			config = 4 - state
			t = min(max((2*e - a - c) // (2*(b - e)) + 1, 0), k) if e != b else k
			M = self.mat(config, t)  # t only matters if in the unstable, that is, if config == 3.
			_, _, new_e, new_b = M([a, c, e, b])
		
		return [new_b if i == bi else new_e if i == ei else vector[i] for i in range(self.zeta)]
	
	def apply_algebraic(self, vector):
		a, b, c, d = self.square
		e = self.source_triangulation.edge_lookup[self.edge_label]
		new_b = vector[b.index] + b.sign() * c.sign() * vector[c.index] * self.power
		new_e = vector[e.index] - e.sign() * c.sign() * vector[c.index] * self.power
		return [new_b if i == b.index else new_e if i == self.edge_index else vector[i] for i in range(self.zeta)]
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		# inverse_corner_map = dict((self(corner), corner) for corner in self.corner_map)
		return Spiral(self.target_triangulation, self.source_triangulation, self.edge_label, -self.power)

