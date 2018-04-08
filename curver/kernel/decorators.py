
from functools import wraps

def memoize():
    ''' A decorator that memoizes a method of a class. '''
    def memoizer(function):
        ''' A decorator that memoizes a method of a class and knows if the method is fast. '''
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
    return memoizer

def topological_invariant(function):
    function.topological_invariant = True
    return function
