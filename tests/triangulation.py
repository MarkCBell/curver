
import pickle
import unittest

from hypothesis import given
import hypothesis.strategies as st

import curver
import strategies

class TestTriangulation(unittest.TestCase):
    def assertImplies(self, A, B):
        self.assertTrue(not A or B, msg='AssertionError: %s =/=> %s' % (A, B))
    
    @given(strategies.triangulations())
    def test_pickle(self, triangulation):
        self.assertEqual(triangulation, pickle.loads(pickle.dumps(triangulation)))
    
    @given(st.data())
    def test_hash(self, data):
        triangulation1 = data.draw(strategies.triangulations())
        triangulation2 = data.draw(strategies.triangulations())
        self.assertImplies(triangulation1 == triangulation2, hash(triangulation1) == hash(triangulation2))
    
    @given(strategies.triangulations())
    def test_isometry(self, triangulation):
        self.assertTrue(triangulation.is_isometric_to(triangulation))
        identity = triangulation.id_isometry()
        isometries = triangulation.self_isometries()
        self.assertIn(identity, isometries)
    
    @given(strategies.triangulations())
    def test_sig(self, triangulation):
        self.assertEqual(triangulation, curver.triangulation_from_sig(triangulation.sig()))
    
    @given(strategies.triangulations())
    def test_homology(self, triangulation):
        self.assertEqual(len(triangulation.homology_basis()), 1 - triangulation.euler_characteristic)  # Assumes connected.
    
    @given(strategies.triangulations())
    def test_connected(self, triangulation):
        for encoding in triangulation.all_encodings(1):
            self.assertEqual(triangulation.is_connected(), encoding.target_triangulation.is_connected())

