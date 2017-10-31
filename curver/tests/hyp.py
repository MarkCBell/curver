
from hypothesis import given
from hypothesis.strategies import lists, integers, composite, assume, one_of
import hypothesis.strategies as st
import unittest
import curver

# Things to test:
#  Relations in the mapping class group.
#  Squares of half twists are twists about their boundary.

@composite
def laminations(draw, S):
    geometric = [None] * S.zeta
    for index in S.triangulation.indices:
        _, a, b, _, c, d = S.triangulation.corner_lookup[index].indices + S.triangulation.corner_lookup[~index].indices
        a, b, c, d = [max(i, 0) if i is not None else None for i in [a, b, c, d]]
        if all(geometric[i] is not None for i in [a, b, c, d]):
            lower_ab, upper_ab = abs(a - b), a + b
            parity_ab = (a + b) % 2
            lower_cd, upper_cd = abs(c - d), c + d
            parity_cd = (c + d) % 2
            exclude_lower, exclude_upper = (max(lower_ab, lower_cd), min(upper_ab, upper_cd))
            e = draw(one_of(integers(max_value=exclude_lower),integers(min_value=exclude_upper)).map(lambda x: x + parity_ab if lower_ab < x < upper_ab else x + parity_cd if lower_cd < x < upper_cd else x))
            pass
        elif all(geometric[i] is not None for i in [a, b]):
            lower, upper = abs(a - b), a + b
            parity = (a + b) % 2
            e = draw(integers().map(lambda x: x + parity if lower < x < upper else x))
        elif all(geometric[i] is not None for i in [c, d]):
            lower, upper = abs(c - d), c + d
            parity = (c + d) % 2
            e = draw(integers().map(lambda x: x + parity if lower < x < upper else x))
        else:
             e = draw(integers())
        geometric[index] = e
    print(geometric)
    return S.lamination(geometric)

@composite
def multicurves(draw, S):
    multicurve = draw(laminations(S)).multicurve()
    assume(not multicurve.is_empty())
    return multicurve

@composite
def multiarcs(draw, S):
    multiarc = draw(laminations(S)).multiarc()
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

@composite
def encodings(draw, S):
    return S.triangulation.id_encoding()
    return S.triangulation.encoding(flips)


class TestS_1_1alt(unittest.TestCase):
    S = curver.load('S_1_1')
    def setUp(self):
        self.identity = self.S('')
        self.empty = self.S.triangulation.empty_lamination()
    @given(laminations(S))
    def test_components(self, lamination):
        self.assertEqual(lamination.triangulation.sum([multiplicity * component for component, multiplicity in lamination.components().items()]), lamination)
        self.assertEqual(lamination.triangulation.disjoint_sum([multiplicity * component for component, multiplicity in lamination.components().items()]), lamination)
    @given(curves(S), laminations(S))
    def a_test_slope_twist(self, curve, lamination):
        assume(curve.intersection(lamination) > 0)
        slope = curve.slope(lamination)
        self.assertTrue(-1 <= slope <= 1 or curve.slope(self.encode_twist()(lamination)) == slope - 1)
    @given(encodings(S))
    def a_test_package(self, encoding):
        self.assertEqual(encoding, encoding.source_triangulation.encode(encoding.package()))
    @given(laminations(S), encodings(S))
    def a_test_components(self, lamination, encoding):
        self.assertEqual(set(encoding(lamination).components()), {encoding(component) for component in lamination.components})
    @given(curves(S), laminations(S))
    def a_test_crush(self, curve, lamination):
        crush = curve.crush()
        self.assertEqual(crush(lamination), crush(curve.encode_twist()(lamination)))
    @given(laminations(S), laminations(S), encodings(S))
    def a_test_intersection(self, lamination, lamination2, encoding):
        self.assertGreaterEqual(lamination.intersection(lamination2), 0)
        self.assertEqual(lamination.intersection(lamination2), encoding(lamination).intersection(encoding(lamination)))
    @given(curves(S), laminations(S), laminations(S), encodings(S))
    def a_test_relative_twist(self, curve, lamination, lamination2, encoding):
        self.assertEqual(curve.relative_twist(lamination, lamination2), encoding(curve).relative_twist(encoding(lamination), encoding(lamination2)))
    @given(arcs(S))
    def test_halftwist(self, arc):
        try:
            halftwist = arc.encode_halftwist()
            self.assertEqual(halftwist**2, arc.boundary.encode_twist())
        except curver.AssumptionError:
            pass

if __name__ == '__main__':
    unittest.main()

