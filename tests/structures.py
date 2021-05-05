
from hypothesis import settings
import hypothesis.strategies as st
from hypothesis.stateful import Bundle, RuleBasedStateMachine, rule

import curver

@settings(max_examples=100)
class UnionFindRules(RuleBasedStateMachine):
    Unions = Bundle("unions")

    @rule(target=Unions, items=st.lists(elements=st.integers(), min_size=1, unique=True))
    def initialize(self, items):
        return curver.kernel.UnionFind(items)

    @rule(data=st.data(), union=Unions)
    def find(self, data, union):
        a = data.draw(st.sampled_from(union.items))
        union(a)

    @rule(union=Unions)
    def iterate(self, union):
        assert set(sum(union, [])) == set(union.items)

    @rule(union=Unions, data=st.data())
    def union2(self, data, union):
        a = data.draw(st.sampled_from(union.items))
        b = data.draw(st.sampled_from(union.items))
        orig_len = len(union)
        union.union2(a, b)
        assert union(a) == union(b)
        assert len(union) in (orig_len, orig_len - 1)

    @rule(data=st.data(), union=Unions)
    def union(self, data, union):
        items = data.draw(st.lists(elements=st.sampled_from(union.items)))
        union.union(*items)
        for item in items:
            assert union(item) == union(items[0])

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

