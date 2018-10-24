
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
    
    @given(st.data())
    def test_intersection(self, data):
        # From Proposition 3.2 of FarbMarg12.
        a = data.draw(strategies.curves())
        b = data.draw(strategies.curves(a.triangulation))
        k = data.draw(st.integers())
        
        self.assertEqual(a.encode_twist(power=k)(b).intersection(b), abs(k) * a.intersection(b)**2)
    
    @given(st.data())
    def test_conjugation(self, data):
        # From Fact 3.7 of FarbMarg12.
        a = data.draw(strategies.curves())
        f = data.draw(strategies.mapping_classes(a.triangulation))
        
        self.assertEqual(f(a).encode_twist(), f * a.encode_twist() * f.inverse())

class TestHalfTwist(unittest.TestCase):
    @given(strategies.arcs().filter(lambda a: a.connects_distinct_vertices()), st.integers(), st.integers())
    @settings(max_examples=2)
    def test_powers(self, arc, power1, power2):
        htwist_i = arc.encode_halftwist(power1)
        htwist_j = arc.encode_halftwist(power2)
        htwist_ij = arc.encode_halftwist(power1 + power2)
        htwist_neg_i = arc.encode_halftwist(-power1)
        
        self.assertEqual(htwist_i * htwist_j, htwist_j * htwist_i)  # Commute.
        self.assertEqual(htwist_i * htwist_j, htwist_ij)  # Additive.
        self.assertEqual(htwist_neg_i, ~htwist_i)  # Inverse.

