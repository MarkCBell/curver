
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver
import strategies

class TestTriangulation(unittest.TestCase):
    @given(strategies.triangulations())
    def test_pickle(self, triangulation):
        self.assertEqual(triangulation, pickle.loads(pickle.dumps(triangulation)))
    
    @given(st.data())
    def test_hash(self, data):
        triangulation1 = data.draw(strategies.triangulations())
        triangulation2 = data.draw(strategies.triangulations())
        self.assertTrue(hash(triangulation1) != hash(triangulation2) or triangulation1 == triangulation2)
    
    @given(st.data())
    def test_flip(self, data):
        triangulation = data.draw(strategies.triangulations())
        edge = data.draw(st.sampled_from(triangulation.edges))
        self.assertEqual(triangulation.surface(), triangulation.flip_edge(edge).surface())
        
    @given(strategies.triangulations())
    def test_isometry(self, triangulation):
        self.assertTrue(triangulation.is_isometric_to(triangulation))
        identity = triangulation.id_isometry()
        isometries = triangulation.self_isometries()
        self.assertIn(identity, isometries)

    @given(strategies.triangulations())
    @settings(deadline=None)
    def test_sig(self, triangulation):
        self.assertEqual(triangulation, curver.triangulation_from_sig(triangulation.sig()))

