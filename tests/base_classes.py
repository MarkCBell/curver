
from hypothesis import given
import hypothesis.strategies as st

import strategies

class TopologicalInvariant(object):
    @given(st.data())
    def test_topological_invariants(self, data):
        item = data.draw(self._strategy())
        h = data.draw(strategies.mappings(item.triangulation))
        image = h(item)
        for method_name in dir(item):
            method = getattr(item, method_name)
            if hasattr(method, 'topological_invariant'):
                image_method = getattr(image, method_name)
                item_property = method()
                image_property = image_method()
                self.assertEqual(item_property, image_property, msg='In %s: %s != %s' % (method_name, item_property, image_property))

