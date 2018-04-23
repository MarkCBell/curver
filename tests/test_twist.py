
from hypothesis import given, settings, assume
import hypothesis.strategies as st
import pickle
import unittest

import curver
import strategies

class TestTwist(unittest.TestCase):
    @given(strategies.curves(), st.integers())
    @settings(max_examples=3)
    def test_pickle(self, curve, power):
        twist_i = curve.encode_twist(power)
        
        self.assertEqual(twist_i, pickle.loads(pickle.dumps(twist_i)))
    
    @given(strategies.curves(), st.integers(), st.integers())
    @settings(max_examples=2)
    def test_powers(self, curve, power1, power2):
        assume(power1 != 0)
        assume(power2 != 0)
        assume(power1 + power2 != 0)
        twist_i = curve.encode_twist(power1)
        twist_j = curve.encode_twist(power2)
        twist_ij = curve.encode_twist(power1 + power2)
        twist_neg_i = curve.encode_twist(-power1)
        
        self.assertEqual(twist_i * twist_j, twist_j * twist_i)  # Commute.
        self.assertEqual(twist_i * twist_j, twist_ij)  # Additive.
        self.assertEqual(twist_neg_i, ~twist_i)  # Inverse.
    
    @given(strategies.curves(), st.integers(min_value=-10, max_value=10))
    def test_encoding(self, curve, power):
        twist = curve.encode_twist(power)
        encoding = curver.kernel.Encoding([subitem for item in twist for subitem in (item.encoding**item.power if isinstance(item, curver.kernel.Twist) else [item])])
        self.assertEqual(twist, encoding)

class TestHalfTwist(unittest.TestCase):
    @given(strategies.arcs().filter(lambda a: a.connects_distinct_vertices()), st.integers(min_value=-10, max_value=10))
    def test_encoding(self, arc, power):
        assume(power != 0)
        halftwist = arc.encode_halftwist(power)
        encoding = curver.kernel.Encoding([subitem for item in halftwist for subitem in (item.encoding**item.power if isinstance(item, curver.kernel.HalfTwist) else [item])])
        self.assertEqual(halftwist, encoding)

