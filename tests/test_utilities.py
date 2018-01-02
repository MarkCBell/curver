
from hypothesis import given
import hypothesis.strategies as st
from string import ascii_lowercase
import unittest

import curver

class TestMisc(unittest.TestCase):
    @given(st.integers(min_value=0))
    def test_b64(self, integer):
        strn = curver.kernel.utilities.b64encode(integer)
        self.assertEqual(curver.kernel.utilities.b64decode(strn), integer)

    @given(st.integers(min_value=0, max_value=1000), st.sets(elements=st.text(ascii_lowercase)))
    def test_string_generator(self, n, skip):
        strns = curver.kernel.utilities.string_generator(n, skip)
        self.assertEqual(len(strns), n)
        self.assertFalse(skip.intersection(strns))
    
    @given(st.data())
    def test_maximum(self, data):
        bound = data.draw(st.integers())
        integers = data.draw(st.lists(elements=st.integers(max_value=bound), min_size=1))
        value = curver.kernel.utilities.maximum(integers, upper_bound=bound)
        self.assertEqual(value, min(max(integers), bound))

