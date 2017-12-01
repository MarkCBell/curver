
from hypothesis import given
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_triangulation import triangulations

import curver

@st.composite
def laminations(draw, triangulation=None, min_weight=None, max_weight=None):
    if triangulation is None: triangulation = draw(triangulations())
    geometric = [None] * triangulation.zeta
    for index in triangulation.indices:
        _, a, b, _, c, d = triangulation.corner_lookup[index].indices + triangulation.corner_lookup[~index].indices
        a, b, c, d = [max(geometric[i], 0) if geometric[i] is not None else None for i in [a, b, c, d]]
        if all(w is not None for w in [a, b, c, d]):
            lower_ab, upper_ab = abs(a - b), a + b
            parity_ab = (a + b) % 2
            lower_cd, upper_cd = abs(c - d), c + d
            parity_cd = (c + d) % 2
            exclude_lower, exclude_upper = (max(lower_ab, lower_cd), min(upper_ab, upper_cd)) if parity_ab != parity_cd else (0, 0)
            e = draw(st.one_of(st.integers(max_value=exclude_lower), st.integers(min_value=exclude_upper)).map(
                lambda x: x + 1 if (x % 2 != parity_ab and lower_ab < x < upper_ab) or (x % 2 != parity_cd and lower_cd < x < upper_cd) else x
                ))
        elif all(w is not None for w in [a, b]):
            lower, upper = abs(a - b), a + b
            parity = (a + b) % 2
            e = draw(st.integers().map(lambda x: x + 1 if x % 2 != parity and lower < x < upper else x))
        elif all(w is not None for w in [c, d]):
            lower, upper = abs(c - d), c + d
            parity = (c + d) % 2
            e = draw(st.integers().map(lambda x: x + 1 if x % 2 != parity and lower < x < upper else x))
        else:
             e = draw(st.integers())
        geometric[index] = e
    return triangulation.lamination(geometric)


@pytest.mark.slow
class TestLamination(unittest.TestCase):
    @given(laminations())
    def test_pickle(self, lamination):
        strn = pickle.dumps(lamination)
        self.assertEqual(lamination, pickle.loads(strn))
    
    @given(st.data())
    def test_hash(self, data):
        lamination1 = data.draw(laminations())
        lamination2 = data.draw(laminations(lamination1.triangulation))
        self.assertTrue(hash(lamination1) != hash(lamination2) or lamination1 == lamination2)
    
    @given(st.data())
    def test_orientation(self, data):
        lamination = data.draw(laminations())
        edge = data.draw(st.sampled_from(lamination.triangulation.edges))
        self.assertEqual(lamination(edge), lamination(~edge))

