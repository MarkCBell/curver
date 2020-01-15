
from hypothesis import given, settings
import hypothesis.strategies as st
import unittest

import strategies
from base_classes import TopologicalInvariant

class TestMultiArc(TopologicalInvariant, unittest.TestCase):
    _strategy = staticmethod(strategies.multiarcs)

class TestArc(TestMultiArc):
    _strategy = staticmethod(strategies.arcs)
    
    @given(st.data())
    @settings(max_examples=25)
    def test_boundary_intersection(self, data):
        arc = data.draw(self._strategy())
        boundary = arc.boundary()
        self.assertEqual(arc.intersection(boundary), 0)
    
    @given(st.data())
    @settings(max_examples=10)
    def test_halftwist(self, data):
        arc = data.draw(self._strategy().filter(lambda a: a.connects_distinct_vertices()))
        self.assertEqual(arc.boundary().encode_twist(), arc.encode_halftwist()**2)

