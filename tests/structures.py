
import hypothesis.strategies as st
from hypothesis.stateful import Bundle, RuleBasedStateMachine, rule

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

class SLPRules(RuleBasedStateMachine):
    SLPs = Bundle('slps')
    
    @rule(target=SLPs, items=st.lists(elements=st.integers()))
    def newslp(self, items):
        return curver.kernel.StraightLineProgram(items)
    
    @rule(target=SLPs, slp=SLPs)
    def copy(self, slp):
        copy = curver.kernel.StraightLineProgram(slp)
        assert list(copy) == list(slp)
        return copy
    
    @rule(target=SLPs, slp1=SLPs, slp2=SLPs)
    def add(self, slp1, slp2):
        added = slp1 + slp2
        assert list(added) == list(slp1) + list(slp2)
        return added
    
    @rule(slp=SLPs, factor=st.integers(min_value=0, max_value=1000))
    def multiply(self, slp, factor):
        assert list(slp * factor) == list(slp) * factor
        assert list(factor * slp) == list(slp * factor)
    
    @rule(data=st.data())
    def getitem(self, data):
        slp = data.draw(self.SLPs.filter(lambda s: len(s) > 0))  # Non-empty.
        index = data.draw(st.integers(min_value=0, max_value=len(slp)-1))
        assert list(slp)[index] == slp[index]
    
    @rule(target=SLPs, slp=SLPs)
    def reverse(self, slp):
        rev = slp.reverse()
        assert list(rev) == list(slp)[::-1]
        assert list(reversed(slp)) == list(slp)[::-1]
        return rev

TestSLP = SLPRules.TestCase

