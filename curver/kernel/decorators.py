
''' A module for decorators. '''

import inspect
from decorator import decorator

@decorator
def memoize(function, *args, **kwargs):
    ''' A decorator that memoizes a function. '''
    
    inputs = inspect.getcallargs(function, *args, **kwargs)  # pylint: disable=deprecated-method
    self = inputs.pop('self', function)  # We test whether function is a method by looking for a `self` argument. If not we store the cache in the function itself.
    
    if not hasattr(self, '_cache'):
        self._cache = dict()
    key = (function.__name__, frozenset(inputs.items()))
    if key not in self._cache:
        try:
            self._cache[key] = function(*args, **kwargs)
        except Exception as error:  # pylint: disable=broad-except
            if isinstance(error, KeyboardInterrupt):
                raise
            self._cache[key] = error
    
    result = self._cache[key]
    if isinstance(result, Exception):
        raise result
    else:
        return result

def memoizable(cls):
    ''' A class decorator that add the 'set_cache' method to a class. '''
    
    def set_cache(self, function, answer, *args, **kwargs):
        ''' Set self._cache so that self.function(*args, **kwargs) returns `answer`. '''
        
        inputs = inspect.getcallargs(function, *args, **kwargs)  # pylint: disable=deprecated-method
        self = inputs.pop('self')
        
        if not hasattr(self, '_cache'):
            self._cache = dict()
        key = (function.__name__, frozenset(inputs.items()))
        self._cache[key] = answer
    
    setattr(cls, 'set_cache', set_cache)
    return cls

def topological_invariant(function):
    ''' Mark this function as a topological invariant.
    
    This is allows it to be picked out by the TopologicalInvariant unittests. '''
    
    function.topological_invariant = True
    return function

def conjugacy_invariant(function):
    ''' Mark this function as a conjugacy invariant.
    
    This is allows it to be picked out by the ConjugacyInvariant unittests. '''
    
    function.conjugacy_invariant = True
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

def decorate_all(function):
    ''' A decorator that applies a function, most likely another decorator, to all public methods of a class. '''
    def decorate(cls):
        ''' A class decorator that applies the given function to every public method. '''
        
        for attr, method in inspect.getmembers(cls):  # there's propably a better way to do this
            if not attr.startswith('_') and inspect.isfunction(method):
                setattr(cls, attr, function(method))
        return cls
    return decorate

