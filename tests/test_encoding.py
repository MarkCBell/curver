
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import unittest

import strategies
import numpy as np

class TestEncoding(unittest.TestCase):
    def assertEqualArray(self, M, N):
        self.assertTrue(np.array_equal(M, N), msg='AssertionError: %s != %s' % (M, N))
    
    @given(strategies.encodings())
    def test_pickle(self, h):
        self.assertEqual(h, pickle.loads(pickle.dumps(h)))
    
    @given(st.data())
    def test_slice(self, data):
        h = data.draw(strategies.encodings())
        i = data.draw(st.integers(min_value=0, max_value=len(h)))
        j = data.draw(st.integers(min_value=0, max_value=len(h)))
        i, j = sorted([i, j])
        self.assertEqual(h[:i] * h[i:j] * h[j:], h)
    
    @given(st.data())
    @settings(max_examples=20)
    def test_inverse(self, data):
        g = data.draw(strategies.encodings())
        h = data.draw(strategies.encodings(g.target_triangulation))
        self.assertEqual(~(~g), g)
        self.assertEqual(~g * ~h, ~(h * g))
    
    @given(st.data())
    def test_homology_matrix(self, data):
        g = data.draw(strategies.encodings())
        h = data.draw(strategies.encodings(g.target_triangulation))
        self.assertEqualArray(h.homology_matrix() * g.homology_matrix(), (h * g).homology_matrix())
    
    @given(strategies.encodings())
    def test_intersection_matrix(self, h):
        self.assertEqualArray(h.intersection_matrix().transpose(), (~h).intersection_matrix())
    
    @given(strategies.encodings())
    def test_simplify(self, h):
        self.assertEqual(h.simplify(), h)
    
    @given(strategies.encodings())
    def test_vertex_map(self, h):
        vertex_map = h.vertex_map()
        self.assertEqual(sorted(vertex_map.keys()), sorted(h.source_triangulation.vertices))
        self.assertEqual(sorted(vertex_map.values()), sorted(h.target_triangulation.vertices))

class TestMappingClass(unittest.TestCase):
    @given(st.data())
    @settings(max_examples=10)
    def test_hash(self, data):
        mcg = data.draw(strategies.mcgs())
        g = data.draw(strategies.mapping_classes(mcg))
        h = data.draw(strategies.mapping_classes(mcg))
        self.assertTrue(hash(g) != hash(h) or g == h)
    
    @given(strategies.mapping_classes())
    # @settings(max_examples=2)
    def test_identity(self, h):
        self.assertEqual(h.is_identity(), h == h.source_triangulation.id_encoding())
    
    @given(strategies.mapping_classes())
    @settings(max_examples=2)
    def test_order(self, h):
        self.assertLessEqual(h.order(), h.source_triangulation.max_order())
        self.assertEqual(h**(h.order()), h.source_triangulation.id_encoding())

