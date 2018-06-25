
from hypothesis import given, settings
import hypothesis.strategies as st
import unittest

from base_classes import TopologicalInvariant
import strategies
import curver

class TestMultiCurve(TopologicalInvariant, unittest.TestCase):
    _strategy = staticmethod(strategies.multicurves)
    
    @given(strategies.multicurves())
    @settings(max_examples=20)
    def test_boundary_intersection(self, multicurve):
        boundary = multicurve.boundary()
        self.assertEqual(multicurve.intersection(boundary), 0)
    
    @given(st.data())
    @settings(max_examples=5)
    def test_boundary_union(self, data):
        multicurve = data.draw(strategies.multicurves().filter(lambda c: not c.is_peripheral()))  # Assume not peripheral.
        lamination = data.draw(strategies.laminations(multicurve.triangulation))
        boundary = multicurve.boundary_union(lamination)
        self.assertEqual(multicurve.intersection(boundary), 0)
        self.assertEqual(lamination.intersection(boundary), 0)
    
    @given(strategies.multicurves())
    @settings(max_examples=20)
    def test_crush(self, multicurve):
        crush = multicurve.crush()
        self.assertEqual(crush.source_triangulation.euler_characteristic, crush.target_triangulation.euler_characteristic)
    
    @given(strategies.multicurves())
    @settings(max_examples=20)
    def test_fills(self, multicurve):
        self.assertEqual(multicurve.is_filling(), multicurve.fills_with(multicurve))
    
    @given(strategies.multicurves())
    @settings(max_examples=20)
    def test_vertex_cycle(self, multicurve):
        vertex_cycle = multicurve.vertex_cycle()
        self.assertTrue(all(0 <= vertex_cycle.dual_weight(edge) <= max(multicurve.dual_weight(edge), 2) for edge in multicurve.triangulation.edges))
    
    @given(st.data())
    @settings(max_examples=2)
    def test_sum(self, data):
        triangulation = data.draw(strategies.triangulations())
        multicurves = data.draw(st.lists(elements=strategies.multicurves(triangulation), min_size=2, max_size=3))
        self.assertIsInstance(multicurves[0] + multicurves[1], curver.kernel.MultiCurve)
        self.assertIsInstance(triangulation.sum(multicurves), curver.kernel.MultiCurve)

class TestCurve(TestMultiCurve):
    _strategy = staticmethod(strategies.curves)
    
    def assertWithinOne(self, x, y):
        self.assertTrue(abs(x - y) <= 1, msg='AssertionError: |%s - %s| > 1' % (x, y))
    
    @given(strategies.curves())
    @settings(max_examples=50)
    def test_boundary_intersection(self, curve):
        boundary = curve.boundary()
        self.assertEqual(curve.intersection(boundary), 0)
    
    @given(st.data())
    @settings(max_examples=10)
    def test_slope(self, data):
        curve = data.draw(strategies.curves().filter(lambda c: not c.is_peripheral()))  # Assume not peripheral.
        lamination = data.draw(strategies.laminations(curve.triangulation).filter(lambda l: curve.intersection(l) > 0))  # Assume intersect.
        slope = curve.slope(lamination)
        twist = curve.encode_twist()
        self.assertTrue(-1 <= slope <= 1 or curve.slope(twist(lamination)) == slope - 1)
    
    @given(st.data())
    @settings(max_examples=1)
    def test_relative_twisting(self, data):
        curve = data.draw(strategies.curves().filter(lambda c: not c.is_peripheral()))  # Assume not peripheral.
        lamination1 = data.draw(strategies.laminations(curve.triangulation).filter(lambda l: curve.intersection(l) > 0))  # Assume intersect.
        power = data.draw(st.integers().filter(lambda p: p))
        lamination2 = curve.encode_twist(power)(lamination1)
        rel_twisting = curve.relative_twisting(lamination1, lamination2)
        self.assertWithinOne(rel_twisting, power)
        
        h = data.draw(strategies.mappings(curve.triangulation))
        self.assertWithinOne(h(curve).relative_twisting(h(lamination1), h(lamination2)), power)

