from hypothesis import given, settings, assume
import hypothesis.strategies as st
import unittest

import strategies

import curver

class TestLoad(unittest.TestCase):
    @given(st.integers(min_value=0), st.integers(min_value=1))
    def test_pair(self, g, n):
        assume(2 - 2*g - n < 0)
        S = curver.load(g, n)
        T = S.triangulation
        self.assertEqual(S.genus, g)
        self.assertEqual(S.num_vertices, n)
    
