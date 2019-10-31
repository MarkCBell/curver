
''' A module of useful, generic functions; including input and output formatting. '''

from itertools import product
from string import ascii_lowercase, ascii_uppercase, digits
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
    
    assert isinstance(n, curver.IntegerType)
    assert skip is None or isinstance(skip, (list, tuple, dict, set))
    
    skip = set() if skip is None else set(skip)
    
    alphabet = ascii_lowercase
    results = []
    i = 0
    while len(results) < n:
        i += 1
        for letters in product(alphabet, repeat=i):
            word = ''.join(letters)
            if word not in skip:
                results.append(word)
            if len(results) >= n: break
    
    return results

def name_objects(objects, skip=None):
    ''' Return a list of pairs (name, object). '''
    
    assert isinstance(objects, (list, tuple))
    
    return zip(string_generator(len(objects), skip), objects)

def cyclic_slice(L, x, y=None):
    ''' Return the sublist of L from x (inclusive) to y (exclusive).
    
    L may be cyclically permuted if needed.
    y may be omitted in which case the entire list (cyclically permuted so x is first) is returned. '''
    
    i = L.index(x)
    L = L[i:] + L[:i]  # x is now L[0].
    j = None if y is None else L.index(y)
    return L[:j]

def minimal(iterable, lower_bound):
    ''' Return the minimal item of iterable but terminate early when given a lower_bound. '''

    def helper():
        ''' A generator that yields elements from iterable up to and including one such that key(item) <= lower_bound. '''
        for item in iterable:
            if item <= lower_bound:
                yield lower_bound
                return
            yield item

    return min(helper())

def maximum(iterable, key=lambda x: x, upper_bound=None):
    ''' Return the maximum of iterable but terminate early when given an upper_bound. '''

    def helper():
        ''' A generator that yields elements from iterable up to and including one such that key(item) >= upper_bound. '''
        for item in iterable:
            yield item
            if upper_bound is not None and key(item) >= upper_bound: return

    return max(helper(), key=key)

def alphanum_key(strn):
    ''' Return a list of string and number chunks from a string. '''
    
    blocks = []
    for chunk in re.split('([0-9]+)', strn):
        try:
            blocks.append(int(chunk))
        except ValueError:
            blocks.append(chunk)
    
    return blocks

class Half(object):
    ''' A class for representing 1/2 in such a way that multiplication preserves types. '''
    def __mul__(self, other):
        if isinstance(other, curver.IntegerType):
            result = other // 2
        else:
            result = other / 2
        if 2*result != other:  # Sanity check.
            raise ValueError('{} is not halvable in its field'.format(other))
        return result
    def __str__(self):
        return '1/2'
    def __repr__(self):
        return str(self)
    def __rmul__(self, other):
        return self * other
    def __call__(self, other):
        return self * other

half = Half()

