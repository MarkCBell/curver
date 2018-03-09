
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import unittest

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
        twist_i = curve.encode_twist(power1)
        twist_j = curve.encode_twist(power2)
        twist_ij = curve.encode_twist(power1 + power2)
        twist_neg_i = curve.encode_twist(-power1)
        
        self.assertEqual(twist_i * twist_j, twist_j * twist_i)  # Commute.
        self.assertEqual(twist_i * twist_j, twist_ij)  # Additive.
        self.assertEqual(twist_neg_i, ~twist_i)  # Inverse.

