
from hypothesis import given
import hypothesis.strategies as st
import pickle
import unittest

import strategies
import curver

class TestFiniteSubgroup(unittest.TestCase):
    @given(st.data())
    def test_trivial_subgroup(self, data):
        T = data.draw(strategies.triangulations())
        
        h = T.id_encoding()
        self.assertEqual(h.subgroup(), curver.kernel.FiniteSubgroup.from_generators({0: h}))
        
        T_signature = [(S.chi, 1, [(True, 1, [0], 1) for _ in range(S.p)]) for S in T.surface().values()]
        self.assertEqual(h.subgroup().quotient_orbifold_signature(), T_signature)
    
    @given(st.data())
    def test_conjugacy(self, data):
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
        G = g.subgroup()
        H = h.subgroup()
        self.assertEqual(G.quotient_orbifold_signature(), H.quotient_orbifold_signature())

