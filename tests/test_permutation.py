
from hypothesis import given
import hypothesis.strategies as st
from math import factorial
import pickle
import unittest

import curver


@st.composite
def indices(draw, N=None):
    if N is None: N = draw(st.integers(min_value=1, max_value=10))
    index = draw(st.integers(min_value=0, max_value=factorial(N)-1))
    return N, index

@st.composite
def permutations(draw, N=None):
    N, index = draw(indices(N))
    return curver.kernel.Permutation.from_index(N, index)

class TestPermutations(unittest.TestCase):
    @given(indices())
    def test_from_index(self, data):
        N, index = data
        perm = curver.kernel.Permutation.from_index(N, index)
        self.assertEqual(len(perm), N)
        self.assertEqual(perm.index(), index)
    
    @given(permutations())
    def test_equal(self, perm):
        self.assertEqual(perm, perm)
    
    @given(permutations())
    def test_involution(self, perm):
        identity = curver.kernel.Permutation.from_index(len(perm), 0)
        self.assertTrue(perm.order() <= 2 or perm != ~perm)
    
    @given(permutations())
    def test_pickle(self, perm):
        strn = pickle.dumps(perm)
        self.assertEqual(perm, pickle.loads(strn))
    
    @given(permutations())
    def test_inverse(self, perm):
        identity = curver.kernel.Permutation.from_index(len(perm), 0)
        self.assertEqual(perm * ~perm, identity)
        self.assertEqual(~perm * perm, identity)
    
    @given(st.data())
    def test_even(self, data):
        perm1 = data.draw(permutations())
        perm2 = data.draw(permutations(len(perm1)))
        self.assertEqual(perm1.is_even() == perm2.is_even(), (perm1 * perm2).is_even())
        self.assertEqual(perm1.is_even() == perm2.is_even(), (perm2 * perm1).is_even())

