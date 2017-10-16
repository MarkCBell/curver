
from hypothesis import given
from hypothesis.strategies import lists, integers, composite, assume
import unittest
import curver

# Things to test:
#  Relations in the mapping class group.
#  Squares of half twists are twists about their boundary.

@composite
def laminations(draw, S):
    geometric = draw(lists(integers(min_value=0, max_value=100).map(lambda x: x * 2), min_size=S.zeta, max_size=S.zeta))
    return S.lamination(geometric)

@composite
def multicurves(draw, S):
    multicurve = draw(laminations2(S)).multicurve()
    assume(not multicurve.is_empty())
    return multicurve

@composite
def multiarcs(draw, S):
    multiarc = draw(laminations2(S)).multiarc()
    assume(not multiarc.is_empty())
    return multiarc

@composite
def curves(draw, S):
    multicurve = draw(multicurves(S))
    return multicurve.peek_component()

@composite
def arcs(draw, S):
    multiarc = draw(multiarcs(S))
    return multiarc.peek_component()


class TestS_1_1alt(unittest.TestCase):
    S = curver.load('S_1_1')
    def setUp(self):
        self.identity = self.S('')
        self.empty = self.S.triangulation.empty_lamination()
    @given(laminations(S))
    def test_mcomponents(self, lamination):
        self.assertEqual(lamination.triangulation.sum([mult*comp for comp, mult in lamination.mcomponents().items()]), lamination)
    def test_slope_twist(self, curve, lamination):
        assume(curve.intersection(lamination) > 0)
        slope = curve.slope(lamination)
        self.assert(-1 <= slope <= 1 or curve.slope(self.encode_twist()(lamination)) == slope - 1)
    def test_package(self, encoding):
        self.assertEqual(encoding, encoding.source_triangulation.encode(encoding.package()))
    def test_components(self, lamination, encoding):
        self.assertEqual(encoding(lamination).components(), {encoding(component) for component in lamination.components})  # !?!
    def test_crush(self, curve, lamination):
        crush = curve.crush()
        self.assertEqual(crush(lamination), crush(curve.encode_twist()(lamination)))
    def test_intersection(self, lamination, lamination2, encoding):
        self.assertGreaterEqual(lamination.intersection(lamination2), 0)
        self.assertEqual(lamination.intersection(lamination2), encoding(lamination).intersection(encoding(lamination)))
    def test_relative_twist(self, curve, lamination, lamination2, encoding):
        self.assertEqual(curve.relative_twist(lamination, lamination2), encoding(curve).relative_twist(encoding(lamination), encoding(lamination2)))
    def test_halftwist(self, arc):
        try:
            halftwist = arc.encode_halftwist()
            self.assertEqual(halftwist**2, arc.boundary.encode_twist())
        except curver.AssumptionError:
            pass

if __name__ == '__main__':
    unittest.main()

