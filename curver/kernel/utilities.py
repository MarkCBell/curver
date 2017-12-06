
''' A module of useful, generic functions; including input and output formatting. '''

from functools import wraps
from itertools import product
from string import ascii_lowercase, ascii_uppercase, digits
from collections import defaultdict
import re

import curver

ALPHABET = digits + ascii_lowercase + ascii_uppercase + '+-'

def b64encode(n):
    ''' Return n in base 64. '''
    
    out = []
    while n:
        n, r = divmod(n, 64)
        out.append(ALPHABET[r])
    return ''.join(out)

def b64decode(strn):
    ''' Return the integer with base 64 encoding strn. '''
    
    return sum(ALPHABET.index(c) * 64**i for i, c in enumerate(strn))

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
        groups = defaultdict(list)
        for item in self.parent:
            groups[self(item)].append(item)
        return iter(groups.values())
    def __len__(self):
        return sum(1 if self.parent[item] == item else 0 for item in self.parent)
    def __repr__(self):
        return str(self)
    def __str__(self):
        return ', '.join('{' + ', '.join(str(item) for item in g) + '}' for g in self)
    def __call__(self, x):
        ''' Find the root of x. Two items are in the same group iff they have the same root. '''
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            x, self.parent[x] = self.parent[x], root
        return root
    def union2(self, x, y):
        ''' Combine the class containing x and the class containing y. '''
        rx, ry = self(x), self(y)
        if self.rank[x] > self.rank[y]:
            self.parent[ry] = rx
        elif self.rank[x] < self.rank[y]:
            self.parent[rx] = ry
        elif rx != ry:
            self.parent[ry] = rx
            self.rank[rx] += 1
    def union(self, *args):
        ''' Combine all of the classes containing the given items. '''
        if len(args) == 1: args = args[0]
        for item in args:
            self.union2(args[0], item)

def memoize(fast=False):
    ''' A decorator that memoizes a method of a class. '''
    def memoizer(function):
        ''' A decorator that memoizes a method of a class and knows if the method is fast. '''
        @wraps(function)
        def caching(self, _is_fast=False):
            ''' The cached version of function.
            
            Note that this docstring will be overwritten with functions docstring by the wraps decorator. '''
            if not hasattr(self, '_cache'):
                self._cache = dict()
            key = function.__name__
            if _is_fast:
                return fast or key in self._cache
            else:
                if key not in self._cache:
                    self._cache[key] = function(self)
                return self._cache[key]
        
        return caching
    return memoizer

def cyclic_slice(L, x, y):
    ''' Return the sublist of L from x (inclusive) to y (exclusive).
    
    L may be cyclically permuted if needed. '''
    i = L.index(x)
    L = L[i:] + L[:i]  # x is now L[0].
    j = L.index(y)
    return L[:j]

def maximum(iterable, key=lambda x: x, upper_bound=None):
    ''' Return the maximum of iterable but terminate early when given an upper_bound. '''

    def helper():
        ''' A generator that yeilds elements from iterable up to and including one such that key(item) >= upper_bound. '''
        for item in iterable:
            yield item
            if upper_bound is not None and key(item) >= upper_bound: return

    return max(helper(), key=key)

def alphanum_key(strn):
    ''' Return a list of string and number chunks from a string. '''
    
    blocks = []
    for chunk in re.split('([0-9]+)', strn):
        try:
            blocks.append(int(strn))
        except ValueError:
            blocks.append(strn)
    
    return blocks

def matrix_vector(M, v):
    ''' Return M*v where v is interpreted as a column vector. '''
    
    return (v * M.T).tolist()[0]

