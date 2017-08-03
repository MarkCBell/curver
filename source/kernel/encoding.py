
''' A module for representing and manipulating maps between Triangulations.

Provides: Encoding. '''

from fractions import Fraction

import curver

NT_TYPE_PERIODIC = 'Periodic'
NT_TYPE_REDUCIBLE = 'Reducible'  # Strictly this  means "reducible and not periodic".
NT_TYPE_PSEUDO_ANOSOV = 'Pseudo-Anosov'

class Encoding(object):
	''' This represents a map between two Triagulations.
	
	If it maps to and from the same triangulation then it represents
	a mapping class. This can be checked using self.is_mapping_class().
	
	The map is given by a sequence of Moves which act from right to left. '''
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
		
		That is, if it maps to the triangulation it came from. '''
		
		return self.source_triangulation == self.target_triangulation
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(self.sequence)  # A backup name.
	def __iter__(self):
		return iter(self.sequence)
	def __len__(self):
		return len(self.sequence)
	def __getitem__(self, value):
		if isinstance(value, slice):
			# It turns out that handling all slices correctly is really hard.
			# We need to be very careful with "empty" slices. As Encodings require
			# non-empty sequences, we have to return just the id_encoding. This
			# ensures the Encoding that we return satisfies:
			#   self == self[:i] * self[i:j] * self[j:]
			# even when i == j.
			
			start = 0 if value.start is None else value.start if value.start >= 0 else len(self) + value.start
			stop = len(self) if value.stop is None else value.stop if value.stop >= 0 else len(self) + value.stop
			if start == stop:
				if 0 <= start < len(self):
					return self.sequence[start].target_triangulation.id_encoding()
				else:
					raise IndexError('list index out of range')
			return Encoding(self.sequence[value])
		elif isinstance(value, curver.IntegerType):
			return self.sequence[value]
		else:
			return NotImplemented
	def package(self):
		''' Return a small amount of info that self.source_triangulation can use to reconstruct this triangulation. '''
		return [item.package() for item in self]
	def __reduce__(self):
		return (create_encoding, (self.source_triangulation, self.package()))
	
	def __eq__(self, other):
		if isinstance(other, Encoding):
			if self.source_triangulation != other.source_triangulation or \
				self.target_triangulation != other.target_triangulation:
				raise ValueError('Cannot compare Encodings between different triangulations.')
			
			return all(self(arc) == other(arc) for arc in self.source_triangulation.edge_arcs()) and \
				all(self(hc) == other(hc) for hc in self.source_triangulation.edge_homologies())
		else:
			return NotImplemented
	def __ne__(self, other):
		return not (self == other)
	def __hash__(self):
		# In fact this hash is perfect unless the surface is S_{1,1}.
		return hash(tuple(weight for arc in self.source_triangulation.edge_arcs() for weight in self(arc)))
	
	def __call__(self, other):
		if self.source_triangulation != other.triangulation:
			raise ValueError('Cannot apply an Encoding to something on a triangulation other than source_triangulation.')
		
		for item in reversed(self.sequence):
			other = item(other)
		
		return other
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
		''' Return the inverse of this encoding. '''
		
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
		
		This encoding must be a mapping class. '''
		
		assert(self.is_mapping_class())
		
		for i in range(1, self.source_triangulation.max_order + 1):
			if self**i == self.source_triangulation.id_encoding():
				return i
		return 0
	
	def is_identity(self):
		''' Return if this encoding is the identity map. '''
		
		return self.is_mapping_class() and self.order() == 1
	
	def is_periodic(self):
		''' Return if this encoding has finite order.
		
		This encoding must be a mapping class. '''
		
		return self.order() > 0
	
	def is_reducible(self):
		''' Return if this encoding is reducible and NOT periodic.
		
		This encoding must be a mapping class. '''
		
		return self.nielsen_thurston_type() == NT_TYPE_REDUCIBLE
	
	def is_pseudo_anosov(self):
		''' Return if this encoding is pseudo-Anosov.
		
		This encoding must be a mapping class. '''
		
		return self.nielsen_thurston_type() == NT_TYPE_PSEUDO_ANOSOV
	
	def nielsen_thurston_type(self):
		''' Return the Nielsen--Thurston type of this encoding.
		
		This encoding must be a mapping class. '''
		
		assert(self.is_mapping_class())
		
		if self.is_periodic():
			return NT_TYPE_PERIODIC
		
		if self.asymptotic_translation_length() == 0:
			return NT_TYPE_REDUCIBLE
		
		return NT_TYPE_PSEUDO_ANOSOV
	
	def asymptotic_translation_length(self):
		''' Return the asymptotic translation length of this mapping class on the curve complex.
		
		From Algorithm 6 of Paper 3. '''
		
		# TODO: 2) Fix these constants.
		N = 1  # Some constant.
		D = 1  # Bowditch bound on denominator.
		c = self.triangulation.edge_arcs()[0].boundary()  # A "short" curve.
		geodesic = c.geodesic((self**N)(c))
		m = geodesic[len(geodesic)//2]  # midpoint
		
		n = m.distance((self**N)(m))  # Numerator.
		d = N  # Denominator.
		return Fraction(n, d).limit_denominator(D)

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

