
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import unittest

import strategies
import numpy as np

class TestEncoding(unittest.TestCase):
    _strategy = staticmethod(strategies.encodings)
    
    def assertEqualArray(self, M, N):
        self.assertTrue(np.array_equal(M, N), msg='AssertionError: %s != %s' % (M, N))
    def assertImplies(self, A, B):
        self.assertTrue(not A or B, msg='AssertionError: %s =/=> %s' % (A, B))
    
    @given(st.data())
    def test_pickle(self, data):
        h = data.draw(self._strategy())
        self.assertEqual(h, pickle.loads(pickle.dumps(h)))
        self.assertEqual(h.sequence, pickle.loads(pickle.dumps(h)).sequence)
    
    @given(st.data())
    def test_hash(self, data):
        g = data.draw(self._strategy())
        h = data.draw(self._strategy(g.source_triangulation))
        self.assertImplies(g == h, hash(g) == hash(h))
    
    @given(st.data())
    def test_slice(self, data):
        h = data.draw(self._strategy())
        i = data.draw(st.integers(min_value=0, max_value=len(h)))
        j = data.draw(st.integers(min_value=0, max_value=len(h)))
        i, j = sorted([i, j])
        self.assertEqual(h[:i] * h[i:j] * h[j:], h)
    
    @given(st.data())
    def test_package(self, data):
        h = data.draw(self._strategy())
        self.assertEqual(h, h.source_triangulation.encode(h.package()))

class TestMapping(TestEncoding):
    _strategy = staticmethod(strategies.mappings)
    
    @given(st.data())
    def test_homology_matrix(self, data):
        g = data.draw(self._strategy())
        h = data.draw(self._strategy(g.target_triangulation))
        self.assertEqualArray(h.homology_matrix().dot(g.homology_matrix()), (h * g).homology_matrix())
    
    @given(st.data())
    def test_intersection_matrix(self, data):
        h = data.draw(self._strategy())
        self.assertEqualArray(h.intersection_matrix().transpose(), (~h).intersection_matrix())
    
    @given(st.data())
    @settings(max_examples=10)
    def test_inverse(self, data):
        g = data.draw(self._strategy())
        h = data.draw(self._strategy(g.target_triangulation))
        self.assertEqual(~(~g), g)
        self.assertEqual(~g * ~h, ~(h * g))
    
    @given(st.data())
    def test_simplify(self, data):
        h = data.draw(self._strategy())
        self.assertEqual(h.simplify(), h)
    
    @given(st.data())
    def test_vertex_map(self, data):
        h = data.draw(self._strategy())
        vertex_map = h.vertex_map()
        self.assertEqual(sorted(vertex_map.keys()), sorted(h.source_triangulation.vertices))
        self.assertEqual(sorted(vertex_map.values()), sorted(h.target_triangulation.vertices))
    
    @given(st.data())
    def test_flip_mapping(self, data):
        h = data.draw(self._strategy())
        self.assertEqual(h, h.flip_mapping())

class TestMappingClass(TestMapping):
    _strategy = staticmethod(strategies.mapping_classes)
    
    @given(st.data())
    # @settings(max_examples=2)
    def test_identity(self, data):
        h = data.draw(self._strategy())
        self.assertEqual(h.is_identity(), h == h.source_triangulation.id_encoding())
    
    @given(st.data())
    # @settings(max_examples=2)
    def test_power(self, data):
        h = data.draw(self._strategy())
        c = data.draw(strategies.curves(h.source_triangulation))
        i = data.draw(st.integers(min_value=-10, max_value=10))
        self.assertEqual(h(c, power=i), (h**i)(c))
    
    @given(st.data())
    @settings(max_examples=2)
    def test_order(self, data):
        h = data.draw(self._strategy())
        self.assertLessEqual(h.order(), h.source_triangulation.max_order())
        self.assertEqual(h**(h.order()), h.source_triangulation.id_encoding())
    
    @given(st.data())
    @settings(max_examples=3)
    def test_conjugacy(self, data):
        h = data.draw(strategies.periodic_mapping_classes())
        self.assertTrue(h.is_periodic())
        
        f = data.draw(strategies.mapping_classes(h.source_triangulation, power_range=1))  # Don't make the word length too large.
        g = ~f * h * f
        self.assertTrue(g.is_conjugate_to(h))

