
from hypothesis import given, assume
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver

@st.composite
def braid_groups(draw):
    p = draw(st.integers(min_value=3, max_value=5))
    return curver.load(0, p)

@st.composite
def mcgs_genus(draw):
    g = draw(st.integers(min_value=1, max_value=2))
    p = draw(st.integers(min_value=1, max_value=5))
    return curver.load(g, p)

@st.composite
def mcgs(draw):
    return draw(st.one_of(braid_groups(), mcgs_genus()))


@pytest.mark.slow
class TestMCG(unittest.TestCase):
    @given(mcgs())
    def test_pickle(self, mcg):
        strn = pickle.dumps(mcg)
        self.assertEqual(mcg, pickle.loads(strn))
    
    @given(mcgs())
    def test_braid_relation(self, mcg):
        p = mcg.triangulation.num_vertices
        g = (mcg.triangulation.euler_characteristic + p - 2) // 2
        for i in range(p - 2):
            self.assertEqual(mcg('s_%ds_%ds_%d' % (i, (i+1) % p, i)), mcg('s_%ds_%ds_%d' % ((i+1) % p, i, (i+1) % p)))
        for i in range(g):
            self.assertEqual(mcg('a_%db_%da_%d' % (i, i, i)), mcg('b_%da_%db_%d' % (i, i, i)))
        for i in range(g-1):
            self.assertEqual(mcg('b_%dc_%db_%d' % (i, i, i)), mcg('c_%db_%dc_%d' % (i, i, i)))
        if g > 0:
            for i in range(1, p):
                self.assertEqual(mcg('p_%db_%dp_%d' % (i, g-1, i)), mcg('b_%dp_%db_%d' % (g-1, i, g-1)))

