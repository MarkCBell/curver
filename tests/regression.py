
import unittest

import curver

class TestRegression(unittest.TestCase):
    def test_correct_component_of_topological_type_marked(self):
        S = curver.load(2, 8)
        c = S.curves['p_1'] + S.curves['p_3'] + S.curves['p_5'] + S.curves['p_7'] + S.arcs['s_2'].boundary() + S.arcs['s_4'].boundary() + S.arcs['s_6'].boundary() + S.arcs['s_2']
        d = S.curves['p_1'] + S.curves['p_3'] + S.curves['p_5'] + S.curves['p_7'] + S.arcs['s_2'].boundary() + S.arcs['s_4'].boundary() + S.arcs['s_6'].boundary() + S.arcs['s_4']
        self.assertNotEqual(c.topological_type(), d.topological_type())

    def test_topological_type_embedded_correctly(self):
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
