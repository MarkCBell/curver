
from hypothesis import given, settings
import hypothesis.strategies as st
import unittest

import curver
import strategies
from base_classes import TopologicalInvariant

class TestMultiArc(TopologicalInvariant, unittest.TestCase):
    _strategy = staticmethod(strategies.multiarcs)

    def test_reg_topological_type(self):
        # Regression test that the multiarcs are different.
        # #---#===#---#===#
        # |  /|   |\  |   |
        # | # |   | # |   |
        # |   |   |   |   |
        # #---#===#---#===#
        #
        # #---#===#---#===#
        # |  /|   |   |   |
        # | # |   | # |   |
        # |   |   |  \|   |
        # #---#===#---#===#
        # This was initially solved by adding in nodes of order two.
        T = curver.create_triangulation([(0, 1, 2), (~1, 3, 4), (~2, 5, ~3), (~4, 6, ~5), (~6, 7, 8), (~7, ~8, 9), (~9, 10, 11), (~10, 12, 13), (~11, 14, ~12), (~13, 15, ~14), (~15, 16, 17), (~16, ~17, ~0)])
        b = T([-1, 0, 0, -1, 0, -1, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0])
        x = b + T.edge_arc(11)
        y = b + T.edge_arc(13)
        self.assertNotEqual(x.topological_type(), y.topological_type())

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

