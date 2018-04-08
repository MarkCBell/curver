
from hypothesis import given, settings
import unittest

import strategies
from base_classes import TopologicalInvariant

class TestMultiArc(TopologicalInvariant, unittest.TestCase):
    _strategy_name = 'multiarcs'

class TestArc(TopologicalInvariant, unittest.TestCase):
    _strategy_name = 'arcs'
    
    @given(strategies.arcs())
    @settings(max_examples=25)
    def test_boundary_intersection(self, arc):
        boundary = arc.boundary()
        self.assertEqual(arc.intersection(boundary), 0)
    
    @given(strategies.arcs().filter(lambda a: a.connects_distinct_vertices()))
    @settings(max_examples=10)
    def test_halftwist(self, arc):
        self.assertEqual(arc.boundary().encode_twist(), arc.encode_halftwist()**2)

