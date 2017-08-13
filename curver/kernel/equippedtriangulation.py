
''' A module for representing triangulations along with laminations and mapping classes on them.

Provides: EquippedTriangulation. '''

from random import choice
import re

import curver

class EquippedTriangulation(object):
	''' This represents a triangulation along with a collection of named laminations and mapping classes on it.
	
	Most importantly this object can construct a mapping class from a string descriptor.
	See self.mapping_class for additional information. '''
	def __init__(self, triangulation, laminations, mapping_classes):
		assert(isinstance(triangulation, curver.kernel.Triangulation))
		assert(isinstance(laminations, (dict, list, tuple)))
		assert(isinstance(mapping_classes, (dict, list, tuple)))
		
		self.triangulation = triangulation
		if isinstance(laminations, dict):
			assert(all(isinstance(key, str) for key in laminations))
			assert(all(isinstance(laminations[key], curver.kernel.Lamination) for key in laminations))
			assert(all(laminations[key].triangulation == self.triangulation for key in laminations))
			self.laminations = laminations
		else:
			assert(all(isinstance(lamination, curver.kernel.Lamination) for lamination in laminations))
			assert(all(lamination.triangulation == self.triangulation for lamination in laminations))
			self.laminations = dict(list(curver.kernel.utilities.name_objects(laminations)))
		
		if isinstance(mapping_classes, dict):
			assert(all(isinstance(key, str) for key in mapping_classes))
			assert(all(isinstance(mapping_classes[key], curver.kernel.Encoding) for key in mapping_classes))
			assert(all(mapping_classes[key].source_triangulation == self.triangulation for key in mapping_classes))
			assert(all(mapping_classes[key].is_mapping_class() for key in mapping_classes))
			assert(all(key.swapcase() not in mapping_classes for key in mapping_classes))
			
			self.pos_mapping_classes = dict(mapping_classes)
			self.neg_mapping_classes = dict((name.swapcase(), self.pos_mapping_classes[name].inverse()) for name in self.pos_mapping_classes)
			self.mapping_classes = dict(list(self.pos_mapping_classes.items()) + list(self.neg_mapping_classes.items()))
		else:
			assert(all(isinstance(mapping_class, curver.kernel.Encoding) for mapping_class in mapping_classes))
			assert(all(mapping_class.source_triangulation == self.triangulation for mapping_class in mapping_classes))
			assert(all(mapping_class.is_mapping_class() for mapping_class in mapping_classes))
			
			self.pos_mapping_classes = dict(list(curver.kernel.utilities.name_objects(mapping_classes)))
			self.neg_mapping_classes = dict((name.swapcase(), self.pos_mapping_classes[name].inverse()) for name in self.pos_mapping_classes)
			self.mapping_classes = dict(list(self.pos_mapping_classes.items()) + list(self.neg_mapping_classes.items()))
		
		self.zeta = self.triangulation.zeta
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		lam_keys = sorted(self.laminations.keys())
		pos_keys = sorted(self.pos_mapping_classes.keys())
		return 'Triangulation with laminations: %s and mapping classes: %s.' % (lam_keys, pos_keys)
	
	def random_word(self, length, positive=True, negative=True, letters=None):
		''' Return a random word of the required length.
		
		The letters to choose from can be specified or, alternatively, the set
		of positive, negative or all (default) mapping classes can be used by using the 
		flags postive and negative. '''
		
		if letters is None:
			pos_keys = sorted(self.pos_mapping_classes.keys())
			neg_keys = sorted(self.neg_mapping_classes.keys())
			
			if positive and negative:
				letters = pos_keys + neg_keys
			elif positive and not negative:
				letters = pos_keys
			elif not positive and negative:
				letters = neg_keys
			else:
				raise TypeError('At least one of positive and negative must be allowed.')
		
		return '.'.join(choice(letters) for _ in range(length))
	
	def decompose_word(self, word):
		''' Return a list of mapping_classes keys whose concatenation is word and the keys are chosen greedly.
		
		Raises a TypeError if the greedy decomposition fails. '''
		
		assert(isinstance(word, str))
		
		# By sorting the available keys, longest first, we ensure that any time we
		# get a match it is as long as possible.
		available_letters = sorted(self.mapping_classes, key=len, reverse=True)
		decomposition = []
		for subword in word.split('.'):
			while subword:
				for letter in available_letters:
					if subword.startswith(letter):
						decomposition.append(letter)
						subword = subword[len(letter):]
						break
				else:
					raise TypeError('After extracting %s, the remaining %s cannot be greedly decomposed as a concatination of self.mapping_classes.' % ('.'.join(decomposition), word))
		
		return decomposition
	
	def mapping_class(self, word):
		''' Return the mapping class corresponding to the given word or a random one of given length if given an integer.
		
		The given word is decomposed using self.decompose_word and the composition
		of the mapping classes involved is returned.
		
		Raises a TypeError if the word does not correspond to a mapping class. '''
		
		if not isinstance(word, str):
			word = self.random_word(word)
		
		name = word  # Record the current word so we can use it later to name the mapping class.
		
		# Remove any whitespace.
		word = word.replace(' ', '')
		
		# Check for balanced parentheses.
		counter = 0
		for letter in word:
			if letter == '(': counter += 1
			if letter == ')': counter -= 1
			if counter < 0: raise TypeError('Unbalanced parentheses.')
		if counter != 0: raise TypeError('Unbalanced parentheses.')
		
		# Expand out parenthesis powers.
		# This can fail with a TypeError.
		old_word = None
		while word != old_word:  # While a change was made.
			old_word = word
			for subword, power in re.findall(r'(\([\w_\.]*\))\^(-?\d+)', word):
				decompose = self.decompose_word(subword[1:-1])
				int_power = int(power)
				if int_power > 0:
					replacement = '.'.join(decompose) * int_power
				else:
					replacement = '.'.join(letter.swapcase() for letter in decompose[::-1]) * abs(int_power)
				word = word.replace(subword + '^' + power, replacement)
		
		# Remove any remaining parenthesis, these do not have a power and are treated as ^1
		word = word.replace('(', '').replace(')', '')
		
		# Expand out powers without parenthesis. Here we use a greedy algorithm and take the
		# longest mapping class occuring before the power. Note that we only do one pass and so
		# only all pure powers to be expanded once, that is 'aBBB^2^3' is not recognised.
		available_letters = sorted(self.mapping_classes, key=len, reverse=True)
		for letter in available_letters:
			for subword, power in re.findall(r'(%s)\^(-?\d+)' % letter, word):
				int_power = int(power)
				word = word.replace(subword + '^' + power, (letter if int_power > 0 else letter.swapcase()) * abs(int_power))
		
		# This can fail with a TypeError.
		sequence = [item for letter in self.decompose_word(word) for item in self.mapping_classes[letter]]
		return curver.kernel.Encoding(sequence) if len(sequence) > 0 else self.triangulation.id_encoding()
	
	def __call__(self, word):
		''' A shortcut for self.mapping_class(...). '''
		return self.mapping_class(word)
	
	def lamination(self, geometric):
		''' Return a new lamination on this surface assigning the specified weight to each edge. '''
		
		return self.triangulation.lamination(geometric)

