
from hypothesis import given, assume, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver

MCGS = [curver.load(g, p) for g in range(0, 5) for p in range(1, 5) if 6*g + 3*p - 6 >= 3]

@st.composite
def mcgs(draw):
    return draw(st.sampled_from(MCGS))


@pytest.mark.slow
class TestMCG(unittest.TestCase):
    @given(mcgs())
    @settings(max_examples=2, deadline=None)
    def test_pickle(self, mcg):
        self.assertEqual(mcg, pickle.loads(pickle.dumps(mcg)))
    
    @given(st.data())
    @settings(max_examples=50, deadline=None)
    def test_curve_relation(self, data):
        mcg = data.draw(mcgs())
        name1 = data.draw(st.sampled_from(sorted(mcg.curves)))
        name2 = data.draw(st.sampled_from(sorted(mcg.curves)))
        curve1 = mcg.curves[name1]
        curve2 = mcg.curves[name2]
        self.assertTrue(
            (curve1.intersection(curve2) != 0 or mcg(name1 + name2) == mcg(name2 + name1)) or \
            (curve1.intersection(curve2) == 1 and mcg(name1 + name2 + name1) == mcg(name2 + name1 + name2))
            )

