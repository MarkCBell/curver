
from hypothesis import given, settings
import hypothesis.strategies as st
import unittest

import strategies

class TestArc(unittest.TestCase):
    @given(strategies.arcs())
    @settings(max_examples=25)
    def test_boundary_intersection(self, arc):
        boundary = arc.boundary()
        self.assertEqual(arc.intersection(boundary), 0)
    
    @given(strategies.arcs().filter(lambda a: a.connects_distinct_vertices()))
    @settings(max_examples=10)
    def test_halftwist(self, arc):
        self.assertEqual(arc.boundary().encode_twist(), arc.encode_halftwist()**2)

class TestMultiArc(unittest.TestCase):
    @given(st.data())
    @settings(max_examples=10)
    def test_types(self, data):
        multiarc = data.draw(strategies.multiarcs())
        encoding = data.draw(strategies.encodings(multiarc.triangulation))
        self.assertEqual(multiarc.is_polygonalisation(), encoding(multiarc).is_polygonalisation())
        self.assertEqual(multiarc.is_triangulation(), encoding(multiarc).is_triangulation())
