
import math

EPSILON = 10**-8

class Vector2(object):
    # Warning: This needs to be updated if the interals of this class ever change.
    __slots__ = ['x', 'y']
    def __init__(self, x, y):
        self.x, self.y = x, y
    def __str__(self):
        return "(%s, %s)" % (self.x, self.y)
    def __repr__(self):
        return str(self)
    def __reduce__(self):
        # Having __slots__ means we need to pickle manually.
        return (self.__class__, self.to_tuple())
    def approx(self, other, epsilon=EPSILON):
        return (self - other).norm2() < epsilon
    def __eq__(self, other):
        raise TypeError('Susceptible to FPE.')
        return self.x == other.x and self.y == other.y
    def __ne__(self, other):
        return not self == other
    def __iter__(self):
        return iter([self.x, self.y])
    def __getitem__(self, item):
        if item == 0: return self.x
        if item == 1: return self.y
    def __setitem__(self, item, value):
        if item == 0: self.x = value
        if item == 1: self.y = value
    def __neg__(self):
        return Vector2(-self.x, -self.y)
    def __add__(self, other):
        return Vector2(self.x + other.x, self.y + other.y)
    def __radd__(self, other):
        return self + other
    def __sub__(self, other):
        return self + -other
    def __mul__(self, other):
        return Vector2(self.x * other, self.y * other)
    def __rmul__(self, other):
        return self * other
    def __div__(self, other):
        return Vector2(self.x / other, self.y / other)
    def __truediv__(self, other):
        return self.__div__(other)
    def dot(self, other):
        return self.x * other.x + self.y * other.y
    def cross(self, other):
        return self.x * other.y - self.y * other.x
    def norm2(self):
        return self.dot(self)
    def norm(self):
        return math.sqrt(self.norm2())
    def to_tuple(self):
        return (self.x, self.y)
    def rotate(self, angle):
        # Angle is in degrees.
        s, c = math.sin(math.pi * angle / 180.0), math.cos(math.pi * angle / 180.0)
        return Vector2(c * self.x - s * self.y, s * self.x + c * self.y)

