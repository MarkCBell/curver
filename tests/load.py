
import unittest

from hypothesis import given, assume
import hypothesis.strategies as st

import curver

class TestLoad(unittest.TestCase):
    @given(st.integers(min_value=0, max_value=3), st.integers(min_value=1, max_value=5))
    def test_pair(self, g, p):
        assume(2 - 2*g - p < 0)
        T = curver.load(g, p).triangulation
        self.assertTrue(T.is_connected())
        
        surface = T.surface()
        self.assertEqual(len(surface), 1)
        
        S = list(surface.values())[0]
        self.assertEqual(S.g, g)
        self.assertEqual(S.p, p)

