
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
			return self.apply_lamination(other)
		if isinstance(other, curver.kernel.HomologyClass):
			return self.apply_homology(other)
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
	
	def apply_lamination(self, lamination):
		geometric = [lamination(self.inverse_index_map[index]) for index in self.source_triangulation.indices]
		return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
	
	def apply_homology(self, homology_class):
		algebraic = [homology_class(self.inverse_label_map[index]) for index in self.source_triangulation.indices]
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		return Isometry(self.target_triangulation, self.source_triangulation, self.inverse_label_map)

class EdgeFlip(Move):
	''' Represents the change to a curve caused by flipping an edge. '''
	def __init__(self, source_triangulation, target_triangulation, edge):
		assert(isinstance(source_triangulation, curver.kernel.Triangulation))
		assert(isinstance(target_triangulation, curver.kernel.Triangulation))
		
		if isinstance(edge, curver.IntegerType): edge = self.source_triangulation.edge_lookup[edge]  # If given an integer instead.
		
		self.source_triangulation = source_triangulation
		self.target_triangulation = target_triangulation
		self.edge = edge
		self.zeta = self.source_triangulation.zeta
		assert(self.source_triangulation.is_flippable(self.edge))
		
		self.square = self.source_triangulation.square(self.edge)
	
	def __str__(self):
		return 'Flip %s' % self.edge
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.edge))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return self.edge.label
	
	def apply_lamination(self, lamination):
		a, b, c, d, e = self.square
		L = lamination  # Shorter name.
		
		# Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
		geometric = list(L.geometric)
		
		if L(e) >= L(a) + L(b) and L(a) >= L(d) and L(b) >= L(c):  # CASE: A(ab)
			geometric[e.index] = L(a) + L(b) - L(e)
		elif L(e) >= L(c) + L(d) and L(d) >= L(a) and L(c) >= L(b):  # CASE: A(cd)
			geometric[e.index] = L(c) + L(d) - L(e)
		elif L(e) <= 0 and L(a) >= L(b) and L(d) >= L(c):  # CASE: D(ad)
			geometric[e.index] = L(a) + L(d) - L(e)
		elif L(e) <= 0 and L(b) >= L(a) and L(c) >= L(d):  # CASE: D(bc)
			geometric[e.index] = L(b) + L(c) - L(e) 
		elif L(a) >= L(b) + L(e) and L(d) >= L(c) + L(e):  # CASE: N(ad)
			geometric[e.index] = L(a) + L(d) - 2*L(e)
		elif L(b) >= L(a) + L(e) and L(c) >= L(d) + L(e):  # CASE: N(bc)
			geometric[e.index] = L(b) + L(c) - 2*L(e)
		elif L(a) + L(b) >= L(e) and L(b) + L(e) >= 2*L(c) + L(a) and L(a) + L(e) >= 2*L(d) + L(b):  # CASE: N(ab)
			geometric[e.index] = (L(a) + L(b) - L(e)) // 2
		elif L(c) + L(d) >= L(e) and L(d) + L(e) >= 2*L(a) + L(c) and L(c) + L(e) >= 2*L(b) + L(d):  # CASE: N(cd)
			geometric[e.index] = (L(c) + L(d) - L(e)) // 2
		else:
			geometric[e.index] = max(L(a) + L(c), L(b) + L(d)) - L(e)
		
		return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
	
	def apply_homology(self, homology_class):
		a, b, c, d, e = self.square
		
		algebraic = list(homology_class)
		# Move the homology on e onto a & b.
		algebraic[a.index] -= a.sign() * homology_class(e)
		algebraic[b.index] -= b.sign() * homology_class(e)
		algebraic[e.index] = 0
		
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		''' Return the inverse of this map. '''
		
		return EdgeFlip(self.target_triangulation, self.source_triangulation, ~self.edge)



class Spiral(Move):
	''' This represents a spiral around a short curve. '''
	# TODO: 4) Completely redo.
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
	
	def apply_lamination(self, leaf):
		# TODO: 4) Make work on all Laminations, not just MultiCurves.
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
		return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
	
	def apply_homology(self, homology_class):
		algebraic = list(homology_class)
		algebraic[self.edge_index] = 0  # TODO 4)
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		# inverse_corner_map = dict((self(corner), corner) for corner in self.corner_map)
		return Spiral(self.target_triangulation, self.source_triangulation, self.edge_label, -self.power)

class Crush(Move):
	# TODO: 2) EVERYTHING!
	pass


