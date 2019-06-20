from hypothesis import given, settings, assume
import hypothesis.strategies as st
import unittest

import strategies

import curver

class TestLoad(unittest.TestCase):
    @given(st.integers(min_value=0, max_value=3), st.integers(min_value=1, max_value=5))
    def test_pair(self, g, n):
        assume(2 - 2*g - n < 0)
        S = curver.load(g, n)
        T = S.triangulation
        self.assertEqual(S.genus, g)
        self.assertEqual(S.num_vertices, n)
    
