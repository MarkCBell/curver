
''' A module for representing permutations in Sym(N). '''

from bisect import bisect
from itertools import combinations
from math import factorial
import numpy as np

class Permutation(object):
    ''' This represents a permutation on 0, 1, ..., N-1. '''
    def __init__(self, perm):
        self.perm = perm
        assert set(self) == set(range(len(self)))
    
    def __str__(self):
        return str(self.perm)
    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, self.perm)
    def __getitem__(self, item):
        return self.perm[item]
    def __call__(self, item):
        return self[item]
    def __iter__(self):
        return iter(self.perm)
    def __len__(self):
        return len(self.perm)
    def __eq__(self, other):
        if isinstance(other, Permutation):
            if len(self) != len(other):
                raise ValueError('Cannot compare permutations defined over different number of elements.')
            
            return self.perm == other.perm
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(tuple(self.perm))
    
    def inverse(self):
        ''' Return the inverse of this permutation. '''
        
        return Permutation(sorted(range(len(self)), key=self))
    def __invert__(self):
        return self.inverse()
    
    @classmethod
    def from_dict(cls, dictionary, ordering=None):
        ''' Return a Permutation from a dictionary.
        
        The order of the elements within the dictionary may also be specified. '''
        
        if ordering is None: ordering = list(dictionary)
        
        index_lookup = dict((item, index) for index, item in enumerate(ordering))
        return cls([index_lookup[dictionary[item]] for item in ordering])
    
    def order(self):
        ''' Return the order of this permutation. '''
        
        identity = Permutation(list(range(len(self))))
        power = self
        i = 1
        while True:
            if power == identity: return i
            i += 1
            power = power * self
        
        raise RuntimeError('Should not be able to reach here.')
    
    def __mul__(self, other):
        if isinstance(other, Permutation):
            if len(self) != len(other):
                raise ValueError('Cannot compose permutations defined over different number of elements.')
            
            return Permutation([self(other(i)) for i in range(len(self))])
        else:
            return NotImplemented
    
    def __pow__(self, n):
        if n < 0: return (~self)**(-n)
        
        perm = self
        result = Permutation(list(range(len(self))))
        while n:
            if n % 2 == 1:
                result = result * perm
                n = n - 1
            perm = perm * perm
            n = n // 2
        return result
    
    @classmethod
    def from_index(cls, N, index):
        ''' Return the permutation in Sym(N) with the given index. '''
        
        P = []
        f = factorial(N)
        symbols = list(range(N))
        while symbols:
            f = f // len(symbols)
            i, index = divmod(index, f)
            P.append(symbols[i])
            symbols = symbols[:i] + symbols[i+1:]
        return cls(P)
    
    def index(self):
        ''' Return the index of this permutation in the (sorted) list of all permutations on this many symbols. '''
        
        symbols = sorted(self.perm)
        index = 0
        for p in self:
            i = bisect(symbols, p) - 1
            index = index * len(symbols) + i
            symbols = symbols[:i] + symbols[i+1:]
        return index
    
    def matrix(self):
        ''' Return the corresponding permutation matrix.
        
        That is, a matrix M such that M * e_i == e_{self[i]}. '''
        
        return np.array([[1 if i == j else 0 for j in self] for i in range(len(self))], dtype=object)
    
    def is_even(self):
        ''' Return whether this permutation is the composition of an even number of transposiions. '''
        
        return sum(1 if a > b else 0 for a, b in combinations(self, r=2) if a > b) % 2 == 0
    
    def cycle_lengths(self):
        ''' Return the sorted list of cycle lengths of this Permutation.
        
        This is a total conjugacy invariant. '''
        
        N = len(self)
        cycle_lengths = []
        seen = set()
        for i in range(N):
            if i not in seen:
                image = self(i)
                seen.add(image)
                for j in range(1, N+1):
                    if image == i:
                        cycle_lengths.append(j)
                        break
                    image = self(image)
                    seen.add(image)
        
        return sorted(cycle_lengths)
    
    def is_conjugate_to(self, other):
        ''' Return whether this permutation in conjugate to other.
        
        Two permutations are conjugate iff they have the same cycle lengths. '''
        
        assert isinstance(other, Permutation)
        if len(self) != len(other):
            return False
        
        return self.cycle_lengths() == other.cycle_lengths()

