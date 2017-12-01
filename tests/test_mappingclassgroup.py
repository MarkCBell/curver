
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
    
    @given(braid_groups())
    def test_braid_relation(self, braid_group):
        p = braid_group.triangulation.num_vertices
        for i in range(p):
            self.assertEqual(braid_group('s_%ds_%ds_%d' % (i, (i+1) % p, i)), braid_group('s_%ds_%ds_%d' % ((i+1) % p, i, (i+1) % p)))
    
    @given(mcgs_genus())
    def test_braid_relation(self, mcg):
        pass

