
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import unittest

from base_classes import TopologicalInvariant
import strategies

class TestLamination(TopologicalInvariant, unittest.TestCase):
    _strategy_name = 'laminations'
    
    @given(strategies.laminations())
    def test_pickle(self, lamination):
        self.assertEqual(lamination, pickle.loads(pickle.dumps(lamination)))
    
    @given(st.data())
    def test_hash(self, data):
        lamination1 = data.draw(strategies.laminations())
        lamination2 = data.draw(strategies.laminations(lamination1.triangulation))
        self.assertTrue(hash(lamination1) != hash(lamination2) or lamination1 == lamination2)
    
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
        encoding = data.draw(strategies.encodings(lamination.triangulation))
        self.assertEqual(set(encoding(lamination).components()), {encoding(component) for component in lamination.components()})
    
    @given(st.data())
    @settings(max_examples=10)
    def test_intersection(self, data):
        lamination1 = data.draw(strategies.laminations())
        lamination2 = data.draw(strategies.laminations(lamination1.triangulation))
        encoding = data.draw(strategies.encodings(lamination1.triangulation))
        self.assertGreaterEqual(lamination1.intersection(lamination2), 0)
        self.assertEqual(lamination1.intersection(lamination2), lamination2.intersection(lamination1))
        self.assertEqual(lamination1.intersection(lamination2), encoding(lamination1).intersection(encoding(lamination2)))

