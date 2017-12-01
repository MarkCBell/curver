
from hypothesis import given
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver

@st.composite
def permutations(draw, N=None):
    if N is None: N = draw(st.integers(min_value=1, max_value=10))
    return curver.kernel.Permutation(draw(st.permutations(range(N))))


class TestPermutation(unittest.TestCase):
    @given(permutations())
    def test_pickle(self, perm):
        strn = pickle.dumps(perm)
        self.assertEqual(perm, pickle.loads(strn))
    
    @given(st.data())
    def test_hash(self, data):
        perm1 = data.draw(permutations())
        perm2 = data.draw(permutations(len(perm1)))
        self.assertTrue(hash(perm1) != hash(perm2) or perm1 == perm2)
    
    @given(permutations())
    def test_from_index(self, perm):
        self.assertEqual(perm, curver.kernel.Permutation.from_index(len(perm), perm.index()))
    
    @given(permutations())
    def test_equal(self, perm):
        self.assertEqual(perm, perm)
    
    @given(st.data())
    def test_inverse(self, data):
        perm1 = data.draw(permutations())
        perm2 = data.draw(permutations(len(perm1)))
        identity = curver.kernel.Permutation.from_index(len(perm1), 0)
        self.assertEqual(perm1 * ~perm1, identity)
        self.assertEqual(~perm1 * perm1, identity)
        self.assertEqual(~(perm1 * perm2), ~perm2 * ~perm1)
    
    @given(permutations())
    def test_involution(self, perm):
        identity = curver.kernel.Permutation.from_index(len(perm), 0)
        self.assertTrue(perm.order() <= 2 or perm != ~perm)
    
    @given(st.data())
    def test_even(self, data):
        perm1 = data.draw(permutations())
        perm2 = data.draw(permutations(len(perm1)))
        self.assertEqual(perm1.is_even() == perm2.is_even(), (perm1 * perm2).is_even())
        self.assertEqual(perm1.is_even() == perm2.is_even(), (perm2 * perm1).is_even())
    
    @given(permutations())
    def test_order(self, perm):
        identity = curver.kernel.Permutation.from_index(len(perm), 0)
        for i in range(1, perm.order()):
            self.assertNotEqual(perm**i, identity)
        self.assertEqual(perm**perm.order(), identity)

