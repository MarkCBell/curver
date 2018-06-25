
from hypothesis import given, settings, assume
import hypothesis.strategies as st
import pickle
import unittest

import strategies

class TestTwist(unittest.TestCase):
    @given(strategies.curves(), st.integers())
    @settings(max_examples=3)
    def test_pickle(self, curve, power):
        assume(power != 0)
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

class TestHalfTwist(unittest.TestCase):
    @given(strategies.arcs().filter(lambda a: a.connects_distinct_vertices()), st.integers(), st.integers())
    @settings(max_examples=2)
    def test_powers(self, arc, power1, power2):
        assume(power1 != 0)
        assume(power2 != 0)
        assume(power1 + power2 != 0)
        htwist_i = arc.encode_halftwist(power1)
        htwist_j = arc.encode_halftwist(power2)
        htwist_ij = arc.encode_halftwist(power1 + power2)
        htwist_neg_i = arc.encode_halftwist(-power1)
        
        self.assertEqual(htwist_i * htwist_j, htwist_j * htwist_i)  # Commute.
        self.assertEqual(htwist_i * htwist_j, htwist_ij)  # Additive.
        self.assertEqual(htwist_neg_i, ~htwist_i)  # Inverse.

