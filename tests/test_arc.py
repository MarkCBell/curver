
from hypothesis import given, settings
import unittest

import strategies

class TestArc(unittest.TestCase):
    @given(strategies.arcs())
    @settings(max_examples=25)
    def test_boundary_intersection(self, arc):
        boundary = arc.boundary()
        self.assertEqual(arc.intersection(boundary), 0)
    
    @given(strategies.arcs().filter(lambda a: a.connects_distinct_vertices()))
    @settings(max_examples=50)
    def test_halftwist(self, arc):
        self.assertEqual(arc.boundary().encode_twist(), arc.encode_halftwist()**2)

