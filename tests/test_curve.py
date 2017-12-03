
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_triangulation import triangulations
from test_mappingclass import mapping_classes

import curver

@st.composite
def curves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    edge = draw(st.sampled_from(triangulation.edges))
    path = []
    seen = set()
    while edge not in seen:
        seen.add(edge)
        path.append(edge)
        edge = ~draw(st.sampled_from(triangulation.corner_lookup[edge.label].edges[1:]))
    start = path.index(edge)
    multicurve = triangulation.lamination_from_cut_sequence(path[start:])
    return multicurve.peek_component()

@st.composite
def multicurves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    curves = draw(st.lists(elements=curves(triangulation), min_size=1))
    return triangulation.sum(curves)


class TestCurve(unittest.TestCase):
    #@given(curves())
    #@settings(max_examples=10, deadline=None)
    #def test_boundary_intersection(self, curve):
    #    boundary = curve.boundary()
    #    self.assertEqual(curve.intersection(boundary), 0)
    
    @given(st.data())
    @settings(max_examples=10, deadline=None)
    def test_topological_type(self, data):
        h = data.draw(mapping_classes())
        curve = data.draw(curves(h.source_triangulation))
        self.assertEqual(curve.topological_type(), h(curve).topological_type())

