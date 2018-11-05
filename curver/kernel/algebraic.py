
from functools import total_ordering
from fractions import Fraction
import sympy as sp
import cypari as cp
import numpy as np
from math import log2 as log

import curver

sp_x = sp.Symbol('x')
cp_x = cp.pari('x')

def exp(x): return 2**x
def log_plus(x): return log(max(1, abs(x)))

def sp_polynomial(coefficients):
    return sp.Poly(coefficients[::-1], sp_x)

def cp_polynomial(coefficients):
    return cp.pari(' + '.join('{}*x^{}'.format(coefficient, index) for index, coefficient in enumerate(coefficients)))

class RealNumberField(object):
    def __init__(self, coefficients, index=-1):  # List of integers and / or Fractions, integer index
        self.coefficients = [Fraction(coefficient) for coefficient in coefficients]
        self.sp_polynomial = sp_polynomial(self.coefficients)
        self.cp_polynomial = cp_polynomial(self.coefficients)
        assert self.cp_polynomial.polisirreducible()
        self.degree = self.cp_polynomial.poldegree()
        self.length = sum(log_plus(coefficient.numerator) + log_plus(coefficient.denominator) for coefficient in self.coefficients)
        self.sp_place = self.sp_polynomial.real_roots()[index]  # Get the real root, will raise an IndexError if no real roots.
    
    def __str__(self):
        return 'QQ(x) / <<{}>> embedding x |--> {}'.format(self.cp_polynomial.lift(), self.N())
    def __repr__(self):
        return str(self)
    def __call__(self, coefficients):
        return RealAlgebraic.from_coefficients(self, coefficients)
    def N(self, prec=64):
        return sp.N(self.sp_place, prec)

@total_ordering
class RealAlgebraic(object):
    def __init__(self, field, cp_mod):
        self.field = field
        self.cp_mod = cp_mod
        self.cp_polynomial = self.cp_mod.lift()
        self.degree = self.cp_polynomial.poldegree()
        self.coefficients = [Fraction(int(self.cp_polynomial.polcoeff(i).numerator()), int(self.cp_polynomial.polcoeff(i).denominator())) for i in range(self.degree+1)]
        self.sp_polynomial = sp_polynomial(self.coefficients)
        self.length = sum(log_plus(coefficient.numerator) + log_plus(coefficient.denominator) + index * self.field.length for index, coefficient in enumerate(self.coefficients))
        self.sp_place = self.sp_polynomial(self.field.sp_place)
    @classmethod
    def from_coefficients(cls, field, coefficients):
        return cls(field, cp_polynomial(coefficients).Mod(field.cp_polynomial))
    @classmethod
    def from_rational(cls, field, rational):
        return cls(field, cp_polynomial([rational]).Mod(field.cp_polynomial))
    def __str__(self):
        return str(self.N())
    def __repr__(self):
        return str(self)
    def __add__(self, other):
        if isinstance(other, RealAlgebraic):
            return RealAlgebraic(self.field, self.cp_mod + other.cp_mod)
        elif isinstance(other, Fraction) or isinstance(other, curver.IntegerType):
            return self + RealAlgebraic.from_rational(self.field, other)
        else:
            return NotImplemented
    def __radd__(self, other):
        return self + other
    def __sub__(self, other):
        return self + (-other)
    def __neg__(self):
        return RealAlgebraic(self.field, -self.cp_mod)
    def __mul__(self, other):
        if isinstance(other, RealAlgebraic):
            return RealAlgebraic(self.field, self.cp_mod * other.cp_mod)
        elif isinstance(other, Fraction) or isinstance(other, curver.IntegerType):
            return self * RealAlgebraic.from_rational(self.field, other)
        else:
            return NotImplemented
    def __rmul__(self, other):
        return self * other
    def __div__(self, other):
        return self._truediv__(other)
    def __truediv__(self, other):
        if other == 0:
            raise ZeroDivisionError('division by zero')
        if isinstance(other, RealAlgebraic):
            return RealAlgebraic(self.field, self.cp_mod / other.cp_mod)
        elif isinstance(other, Fraction) or isinstance(other, curver.IntegerType):
            return self / RealAlgebraic.from_rational(self.field, other)
        else:
            return NotImplemented
    def __rdiv__(self, other):
        return self.__rtruediv__(other)
    def __rtruediv__(self, other):
        if isinstance(other, Fraction) or isinstance(other, curver.IntegerType):
            return RealAlgebraic.from_rational(self.field, other) / self
        else:
            return NotImplemented
    def N(self, prec=64):
        return sp.N(self.sp_place, prec)
    def __float__(self):
        return float(self.N())
    def __int__(self):
        return int(float(self))
    def sign(self):
        scaled_approx = self.N(2*int(self.length+1)) * exp(int(self.length+1))
        if scaled_approx < -1:
            return -1
        elif scaled_approx > 1:
            return +1
        else:
            return 0
    def __eq__(self, other):
        return (self - other).sign() == 0
    def __gt__(self, other):
        return (self - other).sign() == +1

