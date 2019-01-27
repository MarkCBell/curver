
from hypothesis import given, settings
import hypothesis.strategies as st
import pytest
import unittest

from fractions import Fraction
import strategies
import curver

class TestFiniteSubgroup(unittest.TestCase):
    @given(st.data())
    @settings(max_examples=5)
    def test_trivial_subgroup(self, data):
        T = data.draw(strategies.triangulations())
        
        h = T.id_encoding()
        self.assertEqual(h.subgroup(), curver.kernel.FiniteSubgroup.from_generators({0: h}))
        
        T_signature = [(S.chi, 1, [(True, 1, [0], 1) for _ in range(S.p)]) for S in T.surface().values()]
        self.assertEqual(h.subgroup().quotient_orbifold_signature(), T_signature)
    
    @given(st.data())
    @settings(max_examples=3)
    @pytest.mark.slow
    def test_conjugacy(self, data):
        h = data.draw(strategies.periodic_mapping_classes())
        
        f = data.draw(strategies.mapping_classes(h.source_triangulation, power_range=1))  # Don't make the word length too large.
        g = ~f * h * f
        G = g.subgroup()
        H = h.subgroup()
        self.assertEqual(G.quotient_orbifold_signature(), H.quotient_orbifold_signature())

    @given(st.integers(min_value=1, max_value=3))
    @settings(max_examples=3)
    @pytest.mark.slow
    def test_klein(self, genus):
        S = curver.load(genus, 2)
        
        g = S('(a_0.b_0.' + '.'.join('c_{}.b_{}'.format(i, i+1) for i in range(genus-1)) + '.p_1)^{}'.format(genus+1)).simplify()
        h = S('(a_0.b_0.' + '.'.join('c_{}.b_{}'.format(i, i+1) for i in range(genus-1)) + ')^{}.S_1'.format(2*genus+1)).simplify()
        
        K = curver.kernel.FiniteSubgroup.from_generators({'g': g, 'h': h})
        self.assertEqual(len(K), 4)
        
        signature = [(Fraction(-genus, 2), 1, sorted([(False, 2, ['hg' if genus % 2 == 0 else 'g'], 2), (True, 2, ['g'], 2)] + [(False, 2, ['h'], 2)] * (genus+1)))]
        self.assertEqual(K.quotient_orbifold_signature(), signature)
    
    @given(st.integers(min_value=1, max_value=2))
    @settings(max_examples=3)
    @pytest.mark.slow
    def test_dihedral(self, genus):
        S = curver.load(genus, 2)
        
        g = S('a_0.b_0.' + '.'.join('c_{}.b_{}'.format(i, i+1) for i in range(genus-1)) + '.p_1').simplify()
        h = S('(a_0.b_0.' + '.'.join('c_{}.b_{}'.format(i, i+1) for i in range(genus-1)) + ')^{}.S_1'.format(2*genus+1)).simplify()
        
        K = curver.kernel.FiniteSubgroup.from_generators({'g': g, 'h': h})
        self.assertEqual(len(K), 4*(genus+1))
        
        signature = [(Fraction(-genus, 2*(genus+1)), 1, [(False, 2, ['h'], 2*(genus+1)), (False, 2*(genus+1), ['hg'], 2), (True, 2*(genus+1), ['g' * (2*genus + 1)], 2)])]
        self.assertEqual(K.quotient_orbifold_signature(), signature)

