
import pickle
import unittest

from hypothesis import given, settings
import hypothesis.strategies as st

from base_classes import TopologicalInvariant
import strategies

class TestLamination(TopologicalInvariant, unittest.TestCase):
    _strategy = staticmethod(strategies.laminations)
    
    def assertImplies(self, A, B):
        self.assertTrue(not A or B, msg='AssertionError: %s =/=> %s' % (A, B))
    
    @given(st.data())
    def test_pickle(self, data):
        lamination = data.draw(self._strategy())
        self.assertEqual(lamination, pickle.loads(pickle.dumps(lamination)))
    
    @given(st.data())
    def test_hash(self, data):
        lamination1 = data.draw(self._strategy())
        lamination2 = data.draw(self._strategy(lamination1.triangulation))
        self.assertImplies(lamination1 == lamination2, hash(lamination1) == hash(lamination2))
    
    @given(st.data())
    def test_orientation(self, data):
        lamination = data.draw(self._strategy())
        edge = data.draw(st.sampled_from(lamination.triangulation.edges))
        self.assertEqual(lamination(edge), lamination(~edge))

    @given(st.data())
    @settings(max_examples=20)
    def test_components(self, data):
        lamination = data.draw(self._strategy())
        self.assertEqual(lamination.triangulation.sum([multiplicity * component for component, multiplicity in lamination.components().items()]), lamination)
        self.assertEqual(lamination.triangulation.disjoint_sum([multiplicity * component for component, multiplicity in lamination.components().items()]), lamination)
        self.assertEqual(lamination.triangulation.disjoint_sum(lamination.components()), lamination)
        
        for component in lamination.components():
            self.assertEqual(component.intersection(component), 0)
            self.assertEqual(component.components(), {component: 1})
    
    @given(st.data())
    @settings(max_examples=20)
    def test_components_image(self, data):
        lamination = data.draw(self._strategy())
        h = data.draw(strategies.mappings(lamination.triangulation))
        self.assertEqual(set(h(lamination).components()), {h(component) for component in lamination.components()})
    
    @given(st.data())
    @settings(max_examples=10)
    def test_intersection(self, data):
        lamination1 = data.draw(self._strategy())
        lamination2 = data.draw(self._strategy(lamination1.triangulation))
        h = data.draw(strategies.mappings(lamination1.triangulation))
        self.assertGreaterEqual(lamination1.intersection(lamination2), 0)
        self.assertEqual(lamination1.intersection(lamination2), lamination2.intersection(lamination1))
        self.assertEqual(lamination1.intersection(lamination2), h(lamination1).intersection(h(lamination2)))
    
    @given(st.data())
    @settings(max_examples=20)
    def test_boundary_intersection(self, data):
        lamination = data.draw(self._strategy())
        self.assertEqual(lamination.intersection(lamination.boundary()), 0)

