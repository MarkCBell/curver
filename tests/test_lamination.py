
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import unittest

from base_classes import TopologicalInvariant
import strategies

class TestLamination(TopologicalInvariant, unittest.TestCase):
    _strategy = staticmethod(strategies.laminations)
    
    @given(strategies.laminations())
    def test_pickle(self, lamination):
        self.assertEqual(lamination, pickle.loads(pickle.dumps(lamination)))
    
    @given(st.data())
    def test_hash(self, data):
        lamination1 = data.draw(strategies.laminations())
        lamination2 = data.draw(strategies.laminations(lamination1.triangulation))
        self.assertTrue(lamination1 != lamination2 or hash(lamination1) == hash(lamination2))
    
    @given(st.data())
    def test_orientation(self, data):
        lamination = data.draw(strategies.laminations())
        edge = data.draw(st.sampled_from(lamination.triangulation.edges))
        self.assertEqual(lamination(edge), lamination(~edge))

    @given(strategies.laminations())
    @settings(max_examples=20)
    def test_components(self, lamination):
        self.assertEqual(lamination.triangulation.sum([multiplicity * component for component, multiplicity in lamination.components().items()]), lamination)
        self.assertEqual(lamination.triangulation.disjoint_sum([multiplicity * component for component, multiplicity in lamination.components().items()]), lamination)
        for component in lamination.components():
            self.assertEqual(component.intersection(component), 0)
            self.assertEqual(component.components(), {component: 1})
    
    @given(st.data())
    @settings(max_examples=20)
    def test_components_image(self, data):
        lamination = data.draw(strategies.laminations())
        h = data.draw(strategies.mappings(lamination.triangulation))
        self.assertEqual(set(h(lamination).components()), {h(component) for component in lamination.components()})
    
    @given(st.data())
    @settings(max_examples=10)
    def test_intersection(self, data):
        lamination1 = data.draw(strategies.laminations())
        lamination2 = data.draw(strategies.laminations(lamination1.triangulation))
        h = data.draw(strategies.mappings(lamination1.triangulation))
        self.assertGreaterEqual(lamination1.intersection(lamination2), 0)
        self.assertEqual(lamination1.intersection(lamination2), lamination2.intersection(lamination1))
        self.assertEqual(lamination1.intersection(lamination2), h(lamination1).intersection(h(lamination2)))

