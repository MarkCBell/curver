
from math import prod
import unittest

from hypothesis import given, assume, settings
import hypothesis.strategies as st

import curver

class TestMisc(unittest.TestCase):
    @given(st.integers(min_value=0))
    def test_b64(self, integer):
        strn = curver.kernel.utilities.b64encode(integer)
        self.assertEqual(curver.kernel.utilities.b64decode(strn), integer)
    
    @given(st.data())
    def test_cyclic_slice(self, data):
        items = data.draw(st.lists(elements=st.integers(), min_size=2))
        start = data.draw(st.sampled_from(items))
        end = data.draw(st.sampled_from(items))
        assume(start != end)
        
        L = curver.kernel.utilities.cyclic_slice(items, start, end)
        self.assertEqual(L[0], start)
        
        L = curver.kernel.utilities.cyclic_slice(items, start)
        self.assertEqual(L[0], start)
        self.assertEqual(len(L), len(items))
    
    @given(st.lists(elements=st.lists(elements=st.integers(), min_size=1), min_size=1))
    def test_maximin(self, iterables):
        result = curver.kernel.utilities.maximin(*iterables)
        self.assertEqual(result, max(min(iterable) for iterable in iterables))
    
    @given(st.data())
    def test_maximum(self, data):
        bound = data.draw(st.integers())
        integers = data.draw(st.lists(elements=st.integers(max_value=bound), min_size=1))
        value = curver.kernel.utilities.maximum(integers, upper_bound=bound)
        self.assertEqual(value, min(max(integers), bound))
    
    @given(st.lists(elements=st.integers(), min_size=1))
    def test_maxes(self, iterable):
        ''' Return the list of items in iterable whose value is maximal. '''
        
        result = curver.kernel.utilities.maxes(iterable)
        self.assertEqual(set(result), {max(iterable)})
        
        result = curver.kernel.utilities.maxes(iterable, key=lambda x: 1 if x >= 0 else 0)
        self.assertEqual(result, [item for item in iterable if item >= 0])
    
    @given(st.text())
    def test_alphanum_key(self, strn):
        result = curver.kernel.utilities.alphanum_key(strn)
        
        self.assertEqual(''.join(str(x) for x in result), strn)
    
    @given(st.lists(elements=st.integers(), min_size=1))
    def test_product(self, iterable):
        self.assertEqual(curver.kernel.utilities.product(iterable), prod(iterable))

class TestHalf(unittest.TestCase):
    ''' A class for representing 1/2 in such a way that multiplication preserves types. '''
    half = curver.kernel.utilities.half
    
    @given(st.integers())
    @settings(max_examples=100)
    def test_integer(self, integer):
        assume(integer % 2 == 0)
        
        self.assertEqual(self.half * integer, integer // 2)
        self.assertEqual(self.half * integer, self.half(integer))
    
    @given(st.fractions())
    @settings(max_examples=100)
    def test_fraction(self, fraction):
        self.assertEqual((self.half * fraction) * 2, fraction)
        self.assertEqual(self.half * fraction, self.half(fraction))

