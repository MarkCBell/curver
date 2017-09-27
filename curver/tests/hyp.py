
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
		self.assertEqual(sum([mult*comp for comp, mult in lamination.mcomponents()], self.empty), lamination)

if __name__ == '__main__':
	unittest.main()

