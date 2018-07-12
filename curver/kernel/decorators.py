
''' A module for decorators. '''

from decorator import decorator
from functools import wraps
import inspect

@decorator
def memoize(function, *args, **kwargs):
    ''' A decorator that memoizes a method of a class. '''
    
    inputs = inspect.getcallargs(function, *args, **kwargs)
    self = inputs['self']
    key = (function.__name__, frozenset(inputs.items()))
    
    if not hasattr(self, '_cache'):
        self._cache = dict()
    if key not in self._cache:
        self._cache[key] = function(self)
    return self._cache[key]

def topological_invariant(function):
    ''' Mark this function as a topological invariant.
    
    This is allows it to be picked out by the TopologicalInvariant unittests. '''
    
    function.topological_invariant = True
    return function

def ensure(*fs):
    ''' A decorator that specifies properties that the result of a functions should have. '''
    @decorator
    def wrapper(function, *args, **kwargs):
        ''' A decorator that checks that the result of a function has properties fs. '''
        
        result = function(*args, **kwargs)
        data = type('data', (), inspect.getcallargs(function, *args, **kwargs))  # pylint: disable=deprecated-method
        data.result = result
        
        for f in fs:
            assert f(data)
        return result
    
    return wrapper

