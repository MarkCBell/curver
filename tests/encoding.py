
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import unittest

import strategies
import numpy as np
import curver

class TestEncoding(unittest.TestCase):
    _strategy = staticmethod(strategies.encodings)
    
    def assertEqualArray(self, M, N):
        self.assertTrue(np.array_equal(M, N), msg='AssertionError: %s != %s' % (M, N))
    
    @given(st.data())
    def test_pickle(self, data):
        h = data.draw(self._strategy())
        self.assertEqual(h, pickle.loads(pickle.dumps(h)))
    
    @given(st.data())
    def test_hash(self, data):
        g = data.draw(self._strategy())
        h = data.draw(self._strategy(g.source_triangulation))
        self.assertTrue(g != h or hash(g) == hash(h))
    
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
    def test_identity_quotient(self, data):
        T = data.draw(strategies.triangulations())
        
        h = T.id_encoding()
        self.assertEqual(h.order(), 1)
        
        T_signature = [(2 - 2*g - n, 1, [(True, 0, 1) for _ in range(n)]) for (g, n) in T.surface().values()]
        self.assertEqual(h.quotient_orbifold_signature(), T_signature)
    
    @given(st.data())
    @settings(max_examples=2)
    def test_order(self, data):
        h = data.draw(self._strategy())
        self.assertLessEqual(h.order(), h.source_triangulation.max_order())
        self.assertEqual(h**(h.order()), h.source_triangulation.id_encoding())
    
    @given(st.data())
    def test_orbifold(self, data):
        # Periodic mapping classes.
        h = data.draw(st.sampled_from([
            curver.load(0, 6)('s_0.s_1.s_2.s_3.s_4'),
            curver.load(0, 6)('(s_0.s_1.s_2.s_3.s_4)^2'),
            curver.load(0, 6)('(s_0.s_1.s_2.s_3.s_4)^3'),
            curver.load(0, 6)('s_0.s_1.S_3.S_4'),
            curver.load(1, 1)('a_0.b_0'),
            curver.load(1, 1)('a_0.b_0.a_0'),
            curver.load(2, 1)('a_0.b_0.c_0.b_1'),
            curver.load(2, 1)('a_0.b_0.c_0.b_1.a_1'),
            curver.load(2, 2)('a_0.b_0.c_0.b_1.p_1'),
            ]))
        
        f = data.draw(strategies.mapping_classes(h.source_triangulation, power_range=1))  # Don't make the word length too large.
        g = ~f * h * f
        self.assertEqual(g.quotient_orbifold_signature(), h.quotient_orbifold_signature())

