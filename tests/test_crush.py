
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_curve import curves
from test_lamination import laminations
from test_triangulation import triangulations

import curver

from collections import Counter


class TestCrush(unittest.TestCase):
    # Crush is not pickleable since it doesn't implement .package().
    #@given(curves())
    #@settings(max_examples=1, deadline=None)
    #def test_pickle(self, curve):
    #    crush = curve.crush()
    #    self.assertEqual(crush, pickle.loads(pickle.dumps(crush)))
    
    @given(curves())
    @settings(max_examples=10, deadline=None)
    def test_inverse(self, curve):
        crush = curve.crush()
        self.assertEqual(crush, ~(~crush))
    
    @given(curves())
    @settings(deadline=None)
    def test_lift(self, curve):
        crush = curve.crush()
        lift = crush.inverse()
        
        peripheral_curves = [lift.target_triangulation.curve_from_cut_sequence(vertex) for vertex in lift.target_triangulation.vertices]
        lifted_peripherals = [lift(lift.source_triangulation.curve_from_cut_sequence(vertex)) for vertex in lift.source_triangulation.vertices]
        self.assertEqual(Counter(lifted_peripherals), Counter(peripheral_curves + ([] if curve.is_peripheral() else [curve, curve])))
    
    @given(st.data())
    @settings(deadline=None)
    def test_twist(self, data):
        triangulation = data.draw(triangulations())
        curve = data.draw(curves(triangulation))
        lamination = data.draw(curves(triangulation))
        
        crush = curve.crush()
        twist = curve.encode_twist()
        
        self.assertEqual(crush(lamination), crush(twist(lamination)))

