
''' A module for representing and manipulating partial linear functions. '''

import numpy as np
import cypari

import curver

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
        ''' Return the (eigenvalue, eigenvector) of an `interesting` eigenvector of self.action which lives inside of the cone defined by self.condition.
        
        An eigenvector is interesting if its corresponding eigenvalue is:
          - real,
          - greater than 1,
          - irrational, and
          - bigger than all of its Galois conjugates.
        
        Raises a ValueError if it cannot find an interesting vectors in the cone. '''
        
        x = cypari.pari('x')
        
        M = cypari.pari.matrix(*self.action.shape, entries=self.action.flatten())
        
        for polynomial in M.charpoly().factor()[0]:
            degree = int(polynomial.poldegree())
            # TODO: Make this also check for eigenvectors of eigenvalue 1.
            # This is necessary to make projective_invariant_lamination be able to find invariant multicurves.
            if degree > 1:  # It must be irrational to be interesting.
                try:
                    K = curver.kernel.RealNumberField([int(polynomial.polcoeff(i)) for i in range(degree+1)])  # It must be real to be interesting.
                except ValueError:  # No real roots.
                    continue
                
                if K.lmbda > 1:  # It must be > 1 to be interesting.
                    # Compute the kernel:
                    a = x.Mod(polynomial)
                    kernel_basis = (M - a).matker()
                    
                    if len(kernel_basis) == 1:  # Rank 1 kernel.
                        eigenvalue = K.lmbda
                        eigenvector = np.array([K([entry.lift().polcoeff(i) for i in range(degree)]) for entry in kernel_basis[0]], dtype=object)
                        assert np.array_equal(self.action.dot(eigenvector), eigenvalue * eigenvector)
                        if (self.condition.dot(eigenvector) >= 0).all():
                            return eigenvalue, eigenvector
                    else:
                        pass  # We can't handle higher rank kernels yet.
        
        raise ValueError('No interesting eigenvalues in cell')

