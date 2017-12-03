
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_triangulation import triangulations

import curver

@st.composite
def arcs(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    edge = draw(st.sampled_from(triangulation.edges))
    return triangulation.lamination_from_cut_sequence([edge])


class TestArc(unittest.TestCase):
    @given(arcs())
    @settings(max_examples=50, deadline=None)
    def test_boundary_intersection(self, arc):
        boundary = arc.boundary()
        self.assertEqual(arc.intersection(boundary), 0)
    
    @given(arcs())
    @settings(max_examples=10, deadline=None)
    def test_halftwist(self, arc):
        if arc.connects_distinct_vertices():
            self.assertEqual(arc.boundary().encode_twist(), arc.encode_halftwist()**2)

