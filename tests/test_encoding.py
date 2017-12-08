
from hypothesis import given, settings, HealthCheck
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver
import strategies
import numpy as np

class TestEncoding(unittest.TestCase):
    @given(strategies.encodings())
    def test_pickle(self, h):
        self.assertEqual(h, pickle.loads(pickle.dumps(h)))
    
    @given(st.data())
    def test_slice(self, data):
        h = data.draw(strategies.encodings())
        i = data.draw(st.integers(min_value=0, max_value=len(h)))
        j = data.draw(st.integers(min_value=0, max_value=len(h)))
        i, j = sorted([i, j])  # Sorted.
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
        self.assertTrue(np.array_equal(h.homology_matrix() * g.homology_matrix(), (h * g).homology_matrix()))
    
    @given(strategies.encodings())
    def test_intersection_matrix(self, h):
        self.assertTrue(np.array_equal(h.intersection_matrix().transpose(), (~h).intersection_matrix()))

class TestMappingClass(unittest.TestCase):
    @given(st.data())
    @settings(max_examples=20)
    def test_hash(self, data):
        mcg = data.draw(strategies.mcgs())
        g = data.draw(strategies.mapping_classes(mcg))
        h = data.draw(strategies.mapping_classes(mcg))
        self.assertTrue(hash(g) != hash(h) or g == h)
    
    @given(strategies.mapping_classes())
    @settings(max_examples=2)
    def test_order(self, h):
        self.assertLessEqual(h.order(), h.source_triangulation.max_order())
        self.assertEqual(h**(h.order()), h.source_triangulation.id_encoding())

