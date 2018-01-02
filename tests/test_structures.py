
import hypothesis.strategies as st
from hypothesis.stateful import RuleBasedStateMachine, rule

import curver

class UnionFindRules(RuleBasedStateMachine):
    def __init__(self):
        super(UnionFindRules, self).__init__()
        self.initialize([])
    
    @rule(items=st.sets(elements=st.integers()))
    def initialize(self, items):
        self.__items = list(items)
        self.__union_find = curver.kernel.UnionFind(items)
    
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

