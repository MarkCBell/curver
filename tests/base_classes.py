
from hypothesis import given
import hypothesis.strategies as st

from . import strategies

class TopologicalInvariant(object):
    @given(st.data())
    def test_topological_invariants(self, data):
        item = data.draw(self._strategy())  # pylint: disable=no-member
        h = data.draw(strategies.mappings(item.triangulation))
        image = h(item)
        for method_name in dir(item):
            method = getattr(item, method_name)
            if hasattr(method, 'topological_invariant'):
                image_method = getattr(image, method_name)
                item_property = method()
                image_property = image_method()
                self.assertEqual(item_property, image_property, msg='In %s: %s != %s' % (method_name, item_property, image_property))  # pylint: disable=no-member

class ConjugacyInvariant(object):
    @given(st.data())
    def test_conjugacy_invariants(self, data):
        item = data.draw(self._strategy())  # pylint: disable=no-member
        conjugator = data.draw(strategies.mappings(item.source_triangulation))  # Does this need to be mapping_classes?
        image = item.conjugate_by(conjugator.inverse())
        for method_name in dir(item):
            method = getattr(item, method_name)
            if hasattr(method, 'conjugacy_invariant'):
                image_method = getattr(image, method_name)
                item_property = method()
                image_property = image_method()
                self.assertEqual(item_property, image_property, msg='In %s: %s != %s' % (method_name, item_property, image_property))  # pylint: disable=no-member
