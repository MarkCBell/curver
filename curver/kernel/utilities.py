
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

def maximin(*iterables):
    ''' Return the maximum of the minimum, terminating early.
    
    This is equivalent to: max(min(iterable) for iterable in iterables) '''
    
    iterables = iter(iterables)
    try:
        result = min(next(iterables))  # Get the first one through a full evaluation.
    except StopIteration:
        raise ValueError('max() arg is an empty sequence') from None
    
    for iterable in iterables:
        iterable = iter(iterable)
        try:
            best = next(iterable)
        except StopIteration:
            raise ValueError('min() arg is an empty sequence') from None
        
        if best <= result: continue
        
        for item in iterable:
            if item <= result:
                break
            if item < best:  # pylint: disable=consider-using-min-builtin
                best = item
        else:  # We never broke out, so best > result
            result = best
    
    return result

def maximum(iterable, key=lambda x: x, upper_bound=None):
    ''' Return the maximum of iterable but terminate early when given an upper_bound. '''

    def helper():
        ''' A generator that yields items from iterable, if it encounters an item s.t. key(item) >= upper_bound then yields it and exits. '''
        for item in iterable:
            yield item
            if key(item) >= upper_bound: return

    return max(iterable if upper_bound is None else helper(), key=key)

def maxes(iterable, key=lambda x: x):
    ''' Return the list of items in iterable whose value is maximal. '''
    
    iterable = iter(iterable)
    
    try:
        item = next(iterable)
    except StopIteration:
        raise ValueError('maxes() arg is an empty sequence') from None
    
    best = key(item)
    ans = [item]
    for item in iterable:
        value = key(item)
        if value > best:
            ans = [item]
            best = value
        elif value == best:
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
            raise ValueError(f'{other} is not halvable in its field')
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

