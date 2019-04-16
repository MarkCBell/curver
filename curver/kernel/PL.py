
import numpy as np
import cypari

import curver

class PartialLinearFunction(object):
    def __init__(self, action, condition):
        assert isinstance(action, np.ndarray)
        assert isinstance(condition, np.ndarray)
        self.action = action
        self.condition = condition
    
    def __str__(self):
        return 'Action: {}\nCondition: {}'.format(self.action, self.condition)
    def __repr__(self):
        return str(self)
    
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
        
        Raises a ValueError if  if it cannot find an interesting vectors in C. '''
        
        x = cypari.pari('x')
        
        M = cypari.pari.matrix(*self.action.shape, self.action.flatten())
        
        for polynomial in M.charpoly().factor()[0]:
            degree = int(polynomial.poldegree())
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
                        if (self.condition.dot(eigenvector) > 0).all():
                            return eigenvalue, eigenvector
                    else:
                        pass  # We can't handle higher rank kernels yet.
        
        raise ValueError('No interesting eigenvalues in cell')

