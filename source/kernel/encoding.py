
''' A module for representing and manipulating maps between Triangulations.

Provides one class: Encoding. '''

import curver

NT_TYPE_PERIODIC = 'Periodic'
NT_TYPE_REDUCIBLE = 'Reducible'  # Strictly this  means "reducible and not periodic".
NT_TYPE_PSEUDO_ANOSOV = 'Pseudo-Anosov'

class Encoding(object):
	''' This represents a map between two Triagulations.
	
	If it maps to and from the same triangulation then it represents
	a mapping class. This can be checked using self.is_mapping_class().
	
	The map is given by a sequence of EdgeFlips, LinearTransformations
	and Isometries which act from right to left.
	
	>>> import curver
	>>> S = curver.load('S_1_1')
	>>> aB = S.mapping_class('aB')
	>>> bA = S.mapping_class('bA')
	>>> ab = S.mapping_class('ab')
	>>> i = S.mapping_class('')
	>>> a = S.mapping_class('a')
	>>> a
	a
	>>> x = S.triangulation.encode([1])
	>>> x
	[Flip 1]
	'''
	def __init__(self, sequence):
		assert(isinstance(sequence, (list, tuple)))
		assert(len(sequence) > 0)
		assert(all(isinstance(item, curver.kernel.Move) for item in sequence))
		# We used to also test:
		#  assert(all(x.source_triangulation == y.target_triangulation for x, y in zip(sequence, sequence[1:])))
		# However this makes composing Encodings a quadratic time algorithm!
		
		self.sequence = sequence
		
		self.source_triangulation = self.sequence[-1].source_triangulation
		self.target_triangulation = self.sequence[0].target_triangulation
		self.zeta = self.source_triangulation.zeta
	
	def is_mapping_class(self):
		''' Return if this encoding is a mapping class.
		
		That is, if it maps to the triangulation it came from.
		
		>>> aB.is_mapping_class(), bA.is_mapping_class()
		(True, True)
		>>> x.is_mapping_class()
		False
		'''
		
		return self.source_triangulation == self.target_triangulation
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(self.sequence)  # A backup name.
	def __iter__(self):
		return iter(self.sequence)
	def __len__(self):
		return len(self.sequence)
	def package(self):
		''' Return a small amount of info that self.source_triangulation can use to reconstruct this triangulation. '''
		return [item.package() for item in self]
	def __reduce__(self):
		return (create_encoding, (self.source_triangulation, self.package()))
	
	def identify(self):
		''' Return a tuple of integers which uniquely determines this map.
		
		The tuple we return is the intersection numbers of the images of
		the key_curves under this map. This uniquely determines the map
		(assuming we know source_triangulation and target_triangulation)
		by Alexanders trick. '''
		
		return tuple(self(arc) for arc in self.source_triangulation.edge_arcs())
	
	def __eq__(self, other):
		if isinstance(other, Encoding):
			if self.source_triangulation != other.source_triangulation or \
				self.target_triangulation != other.target_triangulation:
				raise ValueError('Cannot compare Encodings between different triangulations.')
			
			return self.identify() == other.identify()
		else:
			return NotImplemented
	def __ne__(self, other):
		return not (self == other)
	def __hash__(self):
		return hash(self.identify())
	
	def __call__(self, other):
		if isinstance(other, curver.kernel.Leaf):
			if self.source_triangulation != other.triangulation:
				raise ValueError('Cannot apply an Encoding to a Leaf on a triangulation other than source_triangulation.')
			
			for item in reversed(self.sequence):
				other = item(other)
			
			return other
		elif isinstance(other, curver.kernel.Lamination):
			return other.__class__(self.target_triangulation, {self(component): multiplicity for multiplicity, component in other})
		else:
			return NotImplemented
	def __mul__(self, other):
		if isinstance(other, Encoding):
			if self.source_triangulation != other.target_triangulation:
				raise ValueError('Cannot compose Encodings over different triangulations.')
			
			return Encoding(self.sequence + other.sequence)
		else:
			return NotImplemented
	def __pow__(self, k):
		assert(self.is_mapping_class())
		
		if k == 0:
			return self.source_triangulation.id_encoding()
		elif k > 0:
			return Encoding(self.sequence * k)
		else:
			return self.inverse()**abs(k)
	
	def inverse(self):
		''' Return the inverse of this encoding.
		
		>>> aB.inverse() == bA, ab == ab.inverse(), i == i.inverse()
		(True, False, True)
		'''
		
		return Encoding([item.inverse() for item in reversed(self.sequence)])
	def __invert__(self):
		return self.inverse()
	
	def closing_isometries(self):
		''' Return all the possible isometries from self.target_triangulation to self.source_triangulation.
		
		These are the maps that can be used to close this into a mapping class. '''
		
		return self.target_triangulation.isometries_to(self.source_triangulation)
	
	def order(self):
		''' Return the order of this mapping class.
		
		If this has infinite order then return 0.
		
		This encoding must be a mapping class.
		
		>>> aB.order(), a.order()
		(0, 0)
		>>> i.order(), ab.order()
		(1, 6)
		'''
		
		assert(self.is_mapping_class())
		
		# We could do:
		# for i in range(1, self.source_triangulation.max_order + 1):
		#	if self**i == self.source_triangulation.id_encoding():
		#		return i
		# But this is quadratic in the order so instead we do:
		arcs = self.source_triangulation.edge_arcs()
		possible_orders = set(range(1, self.source_triangulation.max_order+1))
		for arc in arcs:
			arc_image = arc
			for i in range(1, max(possible_orders)+1):
				arc_image = self(arc_image)
				if arc_image != arc:
					possible_orders.discard(i)
					if not possible_orders: return 0  # No finite orders remain so we are infinite order.
		
		return min(possible_orders)
	
	def is_identity(self):
		''' Return if this encoding is the identity map.
		
		>>> i.is_identity()
		True
		>>> aB.is_identity()
		False
		'''
		
		return self.is_mapping_class() and self.order() == 1
	
	def is_periodic(self):
		''' Return if this encoding has finite order.
		
		This encoding must be a mapping class.
		
		>>> aB.is_periodic(), a.is_periodic()
		(False, False)
		>>> i.is_periodic(), ab.is_periodic()
		(True, True)
		'''
		
		return self.order() > 0
	
	def is_reducible(self):
		''' Return if this encoding is reducible and NOT periodic.
		
		This encoding must be a mapping class.
		
		>>> aB.is_reducible(), a.is_reducible()
		(False, True)
		>>> i.is_reducible(), ab.is_reducible()
		(False, False)
		'''
		
		return self.nielsen_thurston_type() == NT_TYPE_REDUCIBLE
	
	def is_pseudo_anosov(self):
		''' Return if this encoding is pseudo-Anosov.
		
		This encoding must be a mapping class.
		
		>>> aB.is_pseudo_anosov(), a.is_pseudo_anosov()
		(True, False)
		>>> i.is_pseudo_anosov(), ab.is_pseudo_anosov()
		(False, False)
		'''
		
		return self.nielsen_thurston_type() == NT_TYPE_PSEUDO_ANOSOV
	
	def nielsen_thurston_type(self):
		''' Return the Nielsen--Thurston type of this encoding.
		
		This encoding must be a mapping class.
		
		>>> ab.nielsen_thurston_type(), a.nielsen_thurston_type(), aB.nielsen_thurston_type()
		('Periodic', 'Reducible', 'Pseudo-Anosov')
		'''
		
		if self.is_periodic():
			return NT_TYPE_PERIODIC
		
		if self.asymptotic_translation_length() == 0:
			return NT_TYPE_REDUCIBLE
		
		return NT_TYPE_PSEUDO_ANOSOV
	
	def asymptotic_translation_length(self):
		return NotImplemented
	
	def is_conjugate_to(self, other):
		''' Return if this mapping class is conjugate to other.
		
		It would also be straightforward to check if self^i ~~ other^j
		for some i, j.
		
		Both encodings must be mapping classes.
		
		Currently assumes that at least one mapping class is pseudo-Anosov. '''
		
		assert(isinstance(other, Encoding))
		
		# Nielsen-Thurston type is a conjugacy invariant.
		if self.nielsen_thurston_type() != other.nielsen_thurston_type():
			return False
		
		return NotImplemented
	

def create_encoding(source_triangulation, sequence):
	''' Return the encoding defined by sequence starting at source_triangulation.
	
	This is only really here to help with pickling. Users should use
	source_triangulation.encode(sequence) directly. '''
	
	assert(isinstance(source_triangulation, curver.kernel.Triangulation))
	
	return source_triangulation.encode(sequence)

def doctest_globs():
	''' Return the globals needed to run doctest on this module. '''
	
	S = curver.load('S_1_1')
	aB = S.mapping_class('aB')
	bA = S.mapping_class('bA')
	ab = S.mapping_class('ab')
	i = S.mapping_class('')
	a = S.mapping_class('a')
	x = S.triangulation.encode([1])
	
	return {'aB': aB, 'bA': bA, 'ab': ab, 'i': i, 'a': a, 'x': x}

