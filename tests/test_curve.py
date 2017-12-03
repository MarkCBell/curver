
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_triangulation import triangulations
from test_encoding import mapping_classes

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
    
    @pytest.skip('Not written')
    @given(st.data())
    def test_slope(self, data):
        curve = data.draw(curves())
        lamination = data.draw(laminations(curve.triangulation))
        assume(curve.intersection(lamination) > 0)
        slope = curve.slope(lamination)
        twist = curve.encode_twist()
        self.assertTrue(-1 <= slope <= 1 or curve.slope(twist(lamination)) == slope - 1)
    
    @pytest.skip('Not written')
    @given(st.data())
    def test_relative_twist(self, data):
        curve = data.draw(curves())
        lamination1 = data.draw(laminations(curve.triangulation))
        lamination2 = data.draw(laminations(curve.triangulation))
        assume(curve.intersection(lamination1) > 0)
        assume(curve.intersection(lamination2) > 0)
        h = data.draw(encodings(curve.triangulation))
        self.assertEqual(curve.relative_twist(lamination, lamination2), h(curve).relative_twist(h(lamination), h(lamination2)))
    
    @given(st.data())
    @settings(max_examples=10, deadline=None)
    def test_topological_type(self, data):
        curve = data.draw(curves())
        h = data.draw(encodings(curve.triangulation))
        self.assertEqual(curve.topological_type(), h(curve).topological_type())

