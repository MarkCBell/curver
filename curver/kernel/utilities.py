
''' A module of useful, generic functions; including input and output formatting. '''

from functools import partial
from itertools import product
from string import ascii_lowercase
from collections import defaultdict
import re

import curver

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
        return len([item for item in self.parent if self.parent[item] == item])
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

class memoize(object):
    ''' Cache the return value of a method.
    
    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.
    
    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method. '''
    def __init__(self, func):
        self.func = func
    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return partial(self, obj)
    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res

def cyclic_slice(L, x, y):
    ''' Return the sublist of L from x (inclusive) to y (exclusive).
    
    L may be cyclically permuted if needed. '''
    i = L.index(x)
    L = L[i:] + L[:i]  # x is now L[0].
    j = L.index(y)
    return L[:j]

def maximum(iterable, key=lambda x: x, upper_bound=None):
    ''' Return the maximum of iterable but terminate early when given an upper_bound. '''
    
    best, best_value = None, None
    for item in iterable:
        value = key(item)
        if best_value is None or value > best_value:
            best, best_value = item, value
        if upper_bound is not None and best_value >= upper_bound:
            return best
    return best

def alphanum_key(strn):
    ''' Return a list of string and number chunks from a string. '''
    
    def tryint(strn):
        try:
            return int(strn)
        except:
            return strn
    
    return [tryint(chunk) for chunk in re.split('([0-9]+)', strn)]
