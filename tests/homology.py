
from hypothesis import given
import hypothesis.strategies as st
import pickle
import unittest

import strategies

class TestHomologyClass(unittest.TestCase):
    @given(strategies.homology_classes())
    def test_pickle(self, hc):
        self.assertEqual(hc, pickle.loads(pickle.dumps(hc)))
    
    @given(st.data())
    def test_hash(self, data):
        hc1 = data.draw(strategies.homology_classes())
        hc2 = data.draw(strategies.homology_classes(hc1.triangulation))
        self.assertTrue(hc1 != hc2 or hash(hc1) == hash(hc2))
    
    @given(st.data())
    def test_canonical(self, data):
        hc1 = data.draw(strategies.homology_classes())
        hc2 = data.draw(strategies.homology_classes(hc1.triangulation))
        self.assertEqual(hc1.canonical() + hc2.canonical(), (hc1 + hc2).canonical())
        self.assertEqual(hc1.canonical() - hc2.canonical(), (hc1 - hc2).canonical())
        self.assertEqual(-(hc1.canonical()), (-hc1).canonical())
    
    @given(st.data())
    def test_image(self, data):
        hc1 = data.draw(strategies.homology_classes())
        hc2 = data.draw(strategies.homology_classes(hc1.triangulation))
        h = data.draw(strategies.mappings(hc1.triangulation))
        self.assertEqual(h(hc1 + hc2), h(hc1) + h(hc2))
    
    @given(st.data())
    def test_abelian(self, data):
        T = data.draw(strategies.triangulations())
        hcs = data.draw(st.lists(elements=strategies.homology_classes(T), min_size=1))
        self.assertEqual(sum(hcs), sum(data.draw(st.permutations(hcs))))
    
    @given(st.data())
    def test_orientation(self, data):
        hc = data.draw(strategies.homology_classes())
        edge = data.draw(st.sampled_from(hc.triangulation.edges))
        self.assertEqual(hc(edge), -hc(~edge))

