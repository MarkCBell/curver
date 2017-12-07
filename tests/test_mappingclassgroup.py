
from hypothesis import given, assume, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver
import strategies

class TestMCG(unittest.TestCase):
    @given(strategies.mcgs())
    @settings(max_examples=1)
    def test_pickle(self, mcg):
        self.assertEqual(mcg, pickle.loads(pickle.dumps(mcg)))
    
    @given(st.data())
    def test_curve_intersection(self, data):
        mcg = data.draw(strategies.mcgs())
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
        name1 = data.draw(st.sampled_from(sorted(mcg.arcs)))
        name2 = data.draw(st.sampled_from(sorted(mcg.arcs)))
        arc1 = mcg.arcs[name1]
        arc2 = mcg.arcs[name2]
        self.assertEqual(arc1.intersection(arc2), 0)
    
    @given(st.data())
    @settings(max_examples=50)
    def test_curve_relation(self, data):
        mcg = data.draw(strategies.mcgs())
        name1 = data.draw(st.sampled_from(sorted(mcg.curves)))
        name2 = data.draw(st.sampled_from(sorted(mcg.curves)))
        curve1 = mcg.curves[name1]
        curve2 = mcg.curves[name2]
        intersection = curve1.intersection(curve2)
        self.assertTrue(
            (intersection != 0 or mcg(name1 + name2) == mcg(name2 + name1)) or \
            (intersection == 1 and mcg(name1 + name2 + name1) == mcg(name2 + name1 + name2))
            )

    @given(st.data())
    @settings(max_examples=50)
    def test_arc_relation(self, data):
        mcg = data.draw(strategies.mcgs())
        name1 = data.draw(st.sampled_from(sorted(mcg.arcs)))
        name2 = data.draw(st.sampled_from(sorted(mcg.arcs)))
        arc1 = mcg.arcs[name1]
        arc2 = mcg.arcs[name2]
        assume(arc1.connects_distinct_vertices() and arc2.connects_distinct_vertices())
        num_distinct_vertices = len(set(arc1.vertices() + arc2.vertices()))
        
        # We have already tested that arc1.intersection(arc2) == 0.
        self.assertTrue(
            (num_distinct_vertices == 4 and mcg(name1 + name2) == mcg(name2 + name1)) or \
            (num_distinct_vertices == 3 and mcg(name1 + name2 + name1) == mcg(name2 + name1 + name2)) or \
            (num_distinct_vertices == 2)
            )

