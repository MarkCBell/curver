
from hypothesis import given, assume, settings, HealthCheck
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver
import strategies

class TestCurve(unittest.TestCase):
    @given(strategies.curves())
    @settings(max_examples=50)
    def test_boundary_intersection(self, curve):
        boundary = curve.boundary()
        self.assertEqual(curve.intersection(boundary), 0)
    
    @given(st.data())
    @settings(max_examples=2)
    def test_sum(self, data):
        triangulation = data.draw(strategies.triangulations())
        multicurves = data.draw(st.lists(elements=strategies.multicurves(triangulation), min_size=2, max_size=3))
        self.assertIsInstance(multicurves[0] + multicurves[1], curver.kernel.MultiCurve)
        self.assertIsInstance(triangulation.sum(multicurves), curver.kernel.MultiCurve)
    
    @given(st.data())
    @settings(max_examples=10)
    def test_slope(self, data):
        curve = data.draw(strategies.curves())
        assume(not curve.is_peripheral())
        lamination = data.draw(strategies.laminations(curve.triangulation))
        assume(curve.intersection(lamination) > 0)
        slope = curve.slope(lamination)
        twist = curve.encode_twist()
        self.assertTrue(-1 <= slope <= 1 or curve.slope(twist(lamination)) == slope - 1)
    
    @pytest.mark.skip('Too rare')
    @given(st.data())
    def test_relative_twist(self, data):
        curve = data.draw(strategies.curves())
        assume(not curve.is_peripheral())
        lamination1 = data.draw(strategies.laminations(curve.triangulation))
        lamination2 = data.draw(strategies.laminations(curve.triangulation))
        assume(curve.intersection(lamination1) > 0)
        assume(curve.intersection(lamination2) > 0)
        h = data.draw(strategies.encodings(curve.triangulation))
        self.assertEqual(curve.relative_twisting(lamination1, lamination2), h(curve).relative_twisting(h(lamination1), h(lamination2)))
    
    @given(st.data())
    @settings(max_examples=20)
    def test_topological_type(self, data):
        curve = data.draw(strategies.curves())
        h = data.draw(strategies.encodings(curve.triangulation))
        self.assertEqual(curve.topological_type(), h(curve).topological_type())

