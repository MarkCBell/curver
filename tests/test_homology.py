
from hypothesis import given
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_triangulation import triangulations

import curver

@st.composite
def homology_classes(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    algebraic = [draw(st.integers()) for _ in range(triangulation.zeta)]
    return curver.kernel.HomologyClass(triangulation, algebraic)


class TestHomologyClass(unittest.TestCase):
    @given(homology_classes())
    def test_pickle(self, hc):
        self.assertEqual(hc, pickle.loads(pickle.dumps(hc)))
    
    @given(st.data())
    def test_hash(self, data):
        hc1 = data.draw(homology_classes())
        hc2 = data.draw(homology_classes(hc1.triangulation))
        self.assertTrue(hash(hc1) != hash(hc2) or hc1 == hc2)
    
    @given(st.data())
    def test_canonical(self, data):
        hc1 = data.draw(homology_classes())
        hc2 = data.draw(homology_classes(hc1.triangulation))
        self.assertEqual(hc1.canonical() + hc2.canonical(), (hc1 + hc2).canonical())
        self.assertEqual(hc1.canonical() - hc2.canonical(), (hc1 - hc2).canonical())
        self.assertEqual(-(hc1.canonical()), (-hc1).canonical())
    
    @given(st.data())
    def test_abelian(self, data):
        T = data.draw(triangulations())
        hcs = data.draw(st.lists(elements=homology_classes(T), min_size=1))
        self.assertEqual(sum(hcs), sum(data.draw(st.permutations(hcs))))
    
    @given(st.data())
    def test_orientation(self, data):
        hc = data.draw(homology_classes())
        edge = data.draw(st.sampled_from(hc.triangulation.edges))
        self.assertEqual(hc(edge), -hc(~edge))

