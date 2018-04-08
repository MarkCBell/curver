
''' A module for decorators. '''

from functools import wraps

def memoize(function):
    ''' A decorator that memoizes a method of a class. '''
    @wraps(function)
    def caching(self):
        ''' The cached version of function.
        
        Note that this docstring will be overwritten with functions docstring by the wraps decorator. '''
        if not hasattr(self, '_cache'):
            self._cache = dict()
        key = function.__name__
        if key not in self._cache:
            self._cache[key] = function(self)
        return self._cache[key]
    
    return caching

def topological_invariant(function):
    ''' Mark this function as a topological invariant.
    
    This is allows it to be picked out by the TopologicalInvariant unittests. '''
    
    function.topological_invariant = True
    return function
