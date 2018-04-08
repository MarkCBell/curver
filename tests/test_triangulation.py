
from hypothesis import given
import hypothesis.strategies as st
import pickle
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
    
    @given(st.data())
    def test_relabel(self, data):
        triangulation = data.draw(strategies.triangulations())
        label_map = [i if data.draw(st.booleans()) else ~i for i in data.draw(st.permutations(range(triangulation.zeta)))]
        
        T = triangulation.relabel_edges(label_map)
        self.assertTrue(triangulation.is_isometric_to(T))
        self.assertTrue(any(all(isom.label_map[i] == label_map[i] for i in range(triangulation.zeta)) for isom in triangulation.isometries_to(T)))
    
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

