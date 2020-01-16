
from hypothesis import given, settings
import hypothesis.strategies as st

from lamination import TestLamination
import strategies

class TestMultiArc(TestLamination):
    _strategy = staticmethod(strategies.multiarcs)

class TestArc(TestMultiArc):
    _strategy = staticmethod(strategies.arcs)
    
    @given(st.data())
    @settings(max_examples=10)
    def test_halftwist(self, data):
        arc = data.draw(self._strategy().filter(lambda a: a.connects_distinct_vertices()))
        self.assertEqual(arc.boundary().encode_twist(), arc.encode_halftwist()**2)

# Remove the TestLamination class from this namespace to prevent py.test from
# running all of its tests again here.
del TestLamination

