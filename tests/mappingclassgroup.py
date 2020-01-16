
import pickle
import unittest

from hypothesis import given, settings, assume
import hypothesis.strategies as st
import pytest

import strategies

class TestMCG(unittest.TestCase):
    @pytest.mark.skip('Slow.')
    @given(strategies.mcgs())
    @settings(max_examples=1)
    def test_pickle(self, mcg):
        self.assertEqual(mcg, pickle.loads(pickle.dumps(mcg)))
    
    @given(st.data())
    def test_curve_intersection(self, data):
        mcg = data.draw(strategies.mcgs())
        assume(mcg.curves)
        name1 = data.draw(st.sampled_from(sorted(mcg.curves)))
        name2 = data.draw(st.sampled_from(sorted(mcg.curves)))
        curve1 = mcg.curves[name1]
        curve2 = mcg.curves[name2]
        intersection = curve1.intersection(curve2)
        self.assertGreaterEqual(intersection, 0)
        self.assertLessEqual(intersection, 1)
    
    @given(st.data())
    def test_arc_intersection(self, data):
        mcg = data.draw(strategies.mcgs())
        assume(mcg.arcs)
        name1 = data.draw(st.sampled_from(sorted(mcg.arcs)))
        name2 = data.draw(st.sampled_from(sorted(mcg.arcs)))
        arc1 = mcg.arcs[name1]
        arc2 = mcg.arcs[name2]
        self.assertEqual(arc1.intersection(arc2), 0)
    
    @given(st.data())
    @settings(max_examples=30)
    def test_curve_relation(self, data):
        mcg = data.draw(strategies.mcgs())
        assume(mcg.curves)
        name1 = data.draw(st.sampled_from(sorted(mcg.curves)))
        name2 = data.draw(st.sampled_from(sorted(mcg.curves)))
        curve1 = mcg.curves[name1]
        curve2 = mcg.curves[name2]
        intersection = curve1.intersection(curve2)
        if intersection == 0:  # Commute.
            self.assertEqual(mcg(name1 + name2), mcg(name2 + name1))
        elif intersection == 1:  # Braid.
            self.assertEqual(mcg(name1 + name2 + name1), mcg(name2 + name1 + name2))
        else:
            pass

    @given(st.data())
    @settings(max_examples=50)
    def test_arc_relation(self, data):
        mcg = data.draw(strategies.mcgs())
        distinct_end_arcs = sorted(name for name, arc in mcg.arcs.items() if arc.connects_distinct_vertices)
        assume(distinct_end_arcs)
        name1 = data.draw(st.sampled_from(distinct_end_arcs))
        name2 = data.draw(st.sampled_from(distinct_end_arcs))  # Hmm, should we check arc1.intersection(arc2) == 0?
        arc1 = mcg.arcs[name1]
        arc2 = mcg.arcs[name2]
        num_distinct_vertices = len(arc1.vertices().union(arc2.vertices()))
        
        if num_distinct_vertices == 4:  # Commute.
            self.assertEqual(mcg(name1 + name2), mcg(name2 + name1))
        elif num_distinct_vertices == 3:  # Braid.
            self.assertEqual(mcg(name1 + name2 + name1), mcg(name2 + name1 + name2))
        else:
            pass
    
    @given(st.data())
    @settings(max_examples=2)
    def test_expand_word(self, data):
        mcg = data.draw(strategies.mcgs())
        word1 = mcg.random_word(data.draw(st.integers(min_value=0, max_value=10)))
        word2 = mcg.random_word(data.draw(st.integers(min_value=0, max_value=10)))
        power = data.draw(st.integers(min_value=-10, max_value=10))
        
        self.assertEqual(mcg(word1 + word2), mcg(word1) * mcg(word2))
        self.assertEqual(mcg('(%s)^%d' % (word1, power)), mcg(word1)**power)

