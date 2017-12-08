
from hypothesis import given
import hypothesis.strategies as st
from hypothesis.stateful import RuleBasedStateMachine, rule
import pytest
from string import ascii_lowercase
import unittest

import curver
import strategies

class UnionFindRules(RuleBasedStateMachine):
    def __init__(self):
        super(UnionFindRules, self).__init__()
        self.initialize([])
    
    @rule(items=st.sets(elements=st.integers()))
    def initialize(self, items):
        self.__items = list(items)
        self.__union_find = curver.kernel.utilities.UnionFind(items)
    
    @rule(data=st.data())
    def find(self, data):
        a = data.draw(st.sampled_from(self.__items))
        self.__union_find(a)
    
    @rule(data=st.data())
    def iterate(self, data):
        assert set(self.__items) == set(sum(self.__union_find, []))
    
    @rule(data=st.data())
    def union2(self, data):
        a = data.draw(st.sampled_from(self.__items))
        b = data.draw(st.sampled_from(self.__items))
        self.__union_find.union2(a, b)
        assert self.__union_find(a) == self.__union_find(b)
    
    @rule(data=st.data())
    def union(self, data):
        items = data.draw(st.lists(elements=st.sampled_from(self.__items)))
        self.__union_find.union(items)
        for item in items:
            assert self.__union_find(item) == self.__union_find(items[0])

TestUnionFind = UnionFindRules.TestCase

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

