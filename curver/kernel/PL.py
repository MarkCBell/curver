
import numpy as np

import curver

class PartialLinearFunction(object):
    def __init__(self, action, condition):
        assert isinstance(action, np.ndarray)
        assert isinstance(condition, np.ndarray)
        self.action = action
        self.condition = condition
    
    def __call__(self, item):
        if (self.condition.dot(item) < 0).any():
            raise ValueError('Cannot apply a PartialLinearFunction outside of the domain specified by its condition matrix.')
        
        return list(self.condition.dot(item))
    
    def __mul__(self, other):
        assert isinstance(other, PartialLinearFunction)
        
        return PartialLinearFunction(self.action.dot(other.action), np.concatenate([other.condition, self.condition.dot(other.action)]))
    
    def eigenvector(self):
        pass
