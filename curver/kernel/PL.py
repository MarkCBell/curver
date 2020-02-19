
''' A module for representing and manipulating partial linear functions. '''

import numpy as np
import realalg

class PartialLinearFunction(object):
    ''' A Linear function defined on a linear subset of Euclidean space. '''
    def __init__(self, action, condition):
        assert isinstance(action, np.ndarray)
        assert isinstance(condition, np.ndarray)
        self.action = np.array([row for index, row in enumerate(action) if index == 0 or row.any()])
        self.condition = np.array([row for index, row in enumerate(condition) if index == 0 or row.any()])
    
    def __str__(self):
        return 'Action: {}\nCondition: {}'.format(self.action, self.condition)
    def __repr__(self):
        return str(self)
    def __eq__(self, other):
        return np.array_equal(self.action, other.action) and np.array_equal(self.condition, other.condition)
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash((tuple(self.action.flatten()), tuple(self.condition.flatten())))
    
    def __call__(self, item):
        if (self.condition.dot(item) < 0).any():
            raise ValueError('Cannot apply a PartialLinearFunction outside of the domain specified by its condition matrix.')
        
        return list(self.action.dot(item))
    
    def __mul__(self, other):
        if other is None:
            return self
        
        assert isinstance(other, PartialLinearFunction)
        
        return PartialLinearFunction(self.action.dot(other.action), np.concatenate([other.condition, self.condition.dot(other.action)]))
    
    def eigenvector(self):
        ''' Return an `interesting` (eigenvalue, eigenvector) pair where self acts on the eigenvector.
        
        See realalg for the definition of `interesting`.
        
        Raises a ValueError if it cannot find an interesting vectors in the cone. '''
        
        for eigenvalue, eigenvector in realalg.eigenvectors(self.action):
            if (self.condition.dot(eigenvector) >= 0).all():
                return eigenvalue, eigenvector
        
        raise ValueError('No interesting eigenvalues in cell')

