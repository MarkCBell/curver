
''' A module of useful, generic functions; including input and output formatting. '''

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

def cyclic_slice(L, x, y=None):
    ''' Return the sublist of L from x (inclusive) to y (exclusive).
    
    L may be cyclically permuted if needed.
    If y is omitted or not present then the entire list (cyclically permuted so x is first) is returned. '''
    
    i = L.index(x)
    L = L[i:] + L[:i]  # x is now L[0].
    try:
        j = None if y is None else L.index(y)
    except ValueError:  # y is not present.
        j = None
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

def maxes(iterable, key=lambda x: x):
    ''' Return the list of items in iterable whose value is maximal. '''
    
    best = None
    ans = []
    for item in iterable:
        value = key(item)
        if best is None or value > best:
            ans = [item]
            best = value
        elif best is not None and value == best:
            ans.append(item)
    
    return ans

def alphanum_key(strn):
    ''' Return a list of string and number chunks from a string. '''
    
    blocks = []
    for chunk in re.split('([0-9]+)', strn):
        try:
            blocks.append(int(chunk))
        except ValueError:
            blocks.append(chunk)
    
    return blocks

def product(iterable, start=1, left=True):
    ''' Return the product of the items in iterable.
    
    If left then multiply items on the left, otherwise multiply on the right.
    If iterable is empty then return start. '''
    
    iterable = iter(iterable)
    
    try:
        result = next(iterable)
    except StopIteration:
        return start
    
    for item in iterable:
        result = item * result if left else result * item
    
    return result

class Half:
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

