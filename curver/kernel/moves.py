
''' A module for representing basic ways of changing triangulations.
These moves can also track how laminations and homology classes move through those changes.

Provides: Move, Isometry, EdgeFlip.

Also provides Spiral whose goal is to allow laminations to be shortened in
polynomial-time. However this is very incomplete at the minute. '''

import curver

class Move(object):
	''' A basic move from one triangulation to another. '''
	def __init__(self, source_triangulation, target_triangulation):
		assert(isinstance(source_triangulation, curver.kernel.Triangulation))
		assert(isinstance(target_triangulation, curver.kernel.Triangulation))
		
		self.source_triangulation = source_triangulation
		self.target_triangulation = target_triangulation
		self.zeta = self.source_triangulation.zeta
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
	
	def inverse(self):
		''' Return the inverse of this move. '''
		
		return NotImplemented
	def apply_lamination(self, lamination):
		''' Return the lamination obtained by mapping the given lamination through this move. '''
		
		return NotImplemented
	def apply_homology(self, homology_class):
		''' Return the homology class obtained by mapping the given homology class through this move. '''
		
		return NotImplemented

class Isometry(Move):
	''' This represents an isometry from one Triangulation to another.
	
	Triangulations can create the isometries between themselves and this
	is the standard way users are expected to create these. '''
	def __init__(self, source_triangulation, target_triangulation, label_map):
		''' This represents an isometry from source_triangulation to target_triangulation.
		
		It is given by a map taking each edge label of source_triangulation to a label of target_triangulation. '''
		
		super(Isometry, self).__init__(source_triangulation, target_triangulation)
		
		assert(isinstance(label_map, dict))
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
		
		return Isometry(self.target_triangulation, self.source_triangulation, self.inverse_label_map)

class EdgeFlip(Move):
	''' Represents the change to a curve caused by flipping an edge. '''
	def __init__(self, source_triangulation, target_triangulation, edge):
		super(EdgeFlip, self).__init__(source_triangulation, target_triangulation)
		
		if isinstance(edge, curver.IntegerType): edge = self.source_triangulation.edge_lookup[edge]  # If given an integer instead.
		
		self.edge = edge
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
		L = lamination  # Shorter name.
		a, b, c, d, e = self.square
		ai, bi, ci, di, ei = [max(L(edge), 0) for edge in self.square]
		
		# Most of the new information matches the old, so we'll take a copy and modify the places that have changed.
		geometric = list(L.geometric)
		
		if L(e) >= ai + bi and ai >= di and bi >= ci:  # CASE: A(ab)
			geometric[e.index] = ai + bi - L(e)
		elif L(e) >= ci + di and di >= ai and ci >= bi:  # CASE: A(cd)
			geometric[e.index] = ci + di - L(e)
		elif L(e) <= 0 and ai >= bi and di >= ci:  # CASE: D(ad)
			geometric[e.index] = ai + di - L(e)
		elif L(e) <= 0 and bi >= ai and ci >= di:  # CASE: D(bc)
			geometric[e.index] = bi + ci - L(e)
		elif L(e) >= 0 and ai >= bi + L(e) and di >= ci + L(e):  # CASE: N(ad)
			geometric[e.index] = ai + di - 2*L(e)
		elif L(e) >= 0 and bi >= ai + L(e) and ci >= di + L(e):  # CASE: N(bc)
			geometric[e.index] = bi + ci - 2*L(e)
		elif ai + bi >= L(e) and bi + L(e) >= 2*ci + ai and ai + L(e) >= 2*di + bi:  # CASE: N(ab)
			geometric[e.index] = (ai + bi - L(e)) // 2
		elif ci + di >= L(e) and di + L(e) >= 2*ai + ci and ci + L(e) >= 2*bi + di:  # CASE: N(cd)
			geometric[e.index] = (ci + di - L(e)) // 2
		else:
			geometric[e.index] = max(ai + ci, bi + di) - L(e)
		
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
		return EdgeFlip(self.target_triangulation, self.source_triangulation, ~self.edge)



class Spiral(Move):
	''' This represents a spiral around a short curve. '''
	# TODO: 3) Completely redo.
	def __init__(self, source_triangulation, target_triangulation, edge_label, power):
		''' This represents spiralling around a short curve passing through edge_label.
		
		The number of spirals is determined by power and the compact form of this move
		means that the amount of work to compute the image of a curve under this move
		is logorithmic in the power.
		
		Because this is a mapping class, source_triangulation and target_triangulation should be equal. '''
		
		super(Spiral, self).__init__(source_triangulation, target_triangulation)
		
		assert(self.source_triangulation == self.target_triangulation)
		
		self.edge_label = edge_label
		self.edge_index = curver.kernel.norm(self.edge_label)
		self.square = self.source_triangulation.square(self.edge_label)
		a, b, c, d, e = self.square
		assert(b == ~d)  # Assert that b == d.
		self.power = power
	
	def __str__(self):
		return 'Spiral^%d %s' % (self.power, self.edge_label)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.edge_label, self.power))
	def package(self):
		''' Return a small amount of data such that self.source_triangulation.encode([data]) == self.encode(). '''
		
		return (self.edge_label, self.power)
	
	def apply_lamination(self, lamination):
		# TODO: 3) Make work on all Laminations, not just MultiCurves.
		# We will begin with an easy case so we can later assume self.power != 0.
		if self.power == 0: return lamination
		
		k = abs(self.power)
		
		a, b, c, d, e = self.square
		ai, bi, ci, di, ei = [lamination(edge) for edge in self.square]
		
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
		
		def F(n, X):  # Note F(0) == Id.
			''' Common move. '''
			(a, b, c, d) = X
			return a, b, (n+1)*c - n*d, n*c + (1-n)*d
		def G(X):
			''' Edge case. '''
			(a, b, c, d) = X
			return a, b, a+b-d, c
		
		X = (ai, ci, bi, ei)
		
		# Take k steps in the graph. Leave unstable after t steps.
		if state == 1:  # Stable.
			Y = F(k, X)
		elif state == 2:  # Transition.
			Y = F(k-1, G(X))
		elif state == 3:  # Unstable.
			if k <= t:  # k steps in unstable
				Y = F(k, X)
			elif k == t+1:  # t steps in unstable, then into transition.
				Y = G(F(t, X))
			elif k == t+2:  # t steps in unstable, into transition, then out of transition.
				Y = G(G(F(t, X)))
			else:  # k > t+2:  # t steps in unstable, into transition, out of transition, then the remaining steps in stable.
				Y = F(k - (t + 2), G(G(F(t, X))))
		
		_, _, new_bi, new_ei = Y  # t only really matters if in the unstable, that is, if state == 3.
		
		if self.power < 0:
			# Reverse the configuration if we are taking a negative power.
			new_bi, new_ei = new_ei, new_bi
		
		geometric = [new_bi if index == b.index else new_ei if index == e.index else lamination(index) for index in self.source_triangulation.indices]
		return lamination.__class__(self.target_triangulation, geometric)  # Avoids promote.
	
	def apply_homology(self, homology_class):
		algebraic = list(homology_class)
		algebraic[self.edge_index] = 0  # TODO 3)
		return curver.kernel.HomologyClass(self.target_triangulation, algebraic)
	
	def inverse(self):
		''' Return the inverse of this isometry. '''
		
		# inverse_corner_map = dict((self(corner), corner) for corner in self.corner_map)
		return Spiral(self.target_triangulation, self.source_triangulation, self.edge_label, -self.power)

