
''' A module of useful, generic functions; including input and output formatting. '''

import curver

from itertools import product, groupby
from string import ascii_lowercase

def string_generator(n, skip=None):
	''' Return a list of n usable names, none of which are in skip. '''
	
	assert(isinstance(n, curver.IntegerType))
	assert(skip is None or isinstance(skip, (list, tuple, dict, set)))
	
	skip = set() if skip is None else set(skip)
	if n < 1: return []
	
	alphabet = ascii_lowercase
	results = []
	i = 0
	while True:
		i += 1
		for letters in product(alphabet, repeat=i):
			word = ''.join(letters)
			if word not in skip:
				results.append(word)
			if len(results) >= n:
				return results

def name_objects(objects, skip=None):
	''' Return a list of pairs (name, object). '''
	
	assert(isinstance(objects, (list, tuple)))
	
	return zip(string_generator(len(objects), skip), objects)

class UnionFind(object):
	''' A fast union--find data structure. Given items must be hashable. '''
	def __init__(self, items):
		self.parent = dict((item, item) for item in items)
		self.rank = dict((item, 0) for item in items)
	def __iter__(self):
		''' Iterate through the groups of self. '''
		return iter(list(g) for _, g in groupby(sorted(self.parent, key=self), key=self))
	def __repr__(self):
		return str(self)
	def __str__(self):
		return ', '.join('{' + ', '.join(str(item) for item in g) + '}' for g in self)
	def __call__(self, x):
		''' Find the root of x. Two items are in the same group iff they have the same root. '''
		if self.parent[x] == x:
			return x
		else:
			self.parent[x] = self(self.parent[x])
			return self.parent[x]
	def union2(self, x, y):
		rx, ry = self(x), self(y)
		if self.rank[x] > self.rank[y]:
			self.parent[ry] = rx
		elif self.rank[x] < self.rank[y]:
			self.parent[rx] = ry
		elif rx != ry:
			self.parent[ry] = rx
			self.rank[rx] += 1
	def union(self, *arg):
		for item in arg:
			self.union2(arg[0], item)

