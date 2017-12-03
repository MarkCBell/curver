
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_mappingclassgroup import mcgs

import curver

@st.composite
def mapping_classes(draw, mcg=None):
    if mcg is None: mcg = draw(mcgs())
    word = draw(st.lists(elements=st.sampled_from(sorted(mcg)), max_size=10).map(lambda letters: '.'.join(letters)))
    return mcg(word)


class TestMappingClass(unittest.TestCase):
    @given(mapping_classes())
    @settings(max_examples=1, deadline=None)
    def test_pickle(self, h):
        self.assertEqual(h, pickle.loads(pickle.dumps(h)))
    
    @given(st.data())
    @settings(max_examples=10, deadline=None)
    def test_hash(self, data):
        mcg = data.draw(mcgs())
        g = data.draw(mapping_classes(mcg))
        h = data.draw(mapping_classes(mcg))
        self.assertTrue(hash(g) != hash(h) or g == h)
    
    @given(st.data())
    @settings(max_examples=1, deadline=None)
    def test_inverse(self, data):
        mcg = data.draw(mcgs())
        g = data.draw(mapping_classes(mcg))
        h = data.draw(mapping_classes(mcg))
        self.assertEqual(~(~g), g)
        self.assertEqual(~g * ~h, ~(h * g))
    
    @given(mapping_classes())
    @settings(max_examples=1, deadline=None)
    def test_order(self, h):
        self.assertLessEqual(h.order(), h.source_triangulation.max_order())
        self.assertEqual(h**(h.order()), h.source_triangulation.id_encoding())

