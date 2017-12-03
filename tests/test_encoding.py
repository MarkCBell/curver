
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_mappingclassgroup import mcgs
from test_triangulation import triangulations

import curver

@st.composite
def mapping_classes(draw, mcg=None):
    if mcg is None: mcg = draw(mcgs())
    word = draw(st.lists(elements=st.sampled_from(sorted(mcg)), max_size=10).map(lambda letters: '.'.join(letters)))
    return mcg(word)

@st.composite
def encodings(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    encoding = triangulation.id_encoding()
    num_flips = draw(st.integers(min_value=0, max_value=20))
    for _ in range(num_flips):
        T = encoding.target_triangulation
        edge = draw(st.sampled_from([edge for edge in T.edges if T.is_flippable(edge)]))
        flip = T.encode_flip(edge)
        encoding = flip * encoding
    
    return encoding


class TestEncoding(unittest.TestCase):
    @given(encodings())
    @settings(max_examples=1, deadline=None)
    def test_pickle(self, h):
        self.assertEqual(h, pickle.loads(pickle.dumps(h)))
    
    @given(st.data())
    @settings(max_examples=1, deadline=None)
    def test_inverse(self, data):
        g = data.draw(encodings())
        h = data.draw(encodings(g.target_triangulation))
        self.assertEqual(~(~g), g)
        self.assertEqual(~g * ~h, ~(h * g))
    
    @given(encodings())
    @settings(deadline=None)
    def test_intersection_matrix(self, h):
        matrix = h.intersection_matrix()
        matrix_transpose = [list(row) for row in zip(*matrix)]
        self.assertEqual(matrix_transpose, (~h).intersection_matrix())

class TestMappingClass(unittest.TestCase):
    @given(st.data())
    @settings(max_examples=10, deadline=None)
    def test_hash(self, data):
        mcg = data.draw(mcgs())
        g = data.draw(mapping_classes(mcg))
        h = data.draw(mapping_classes(mcg))
        self.assertTrue(hash(g) != hash(h) or g == h)
    
    @given(mapping_classes())
    @settings(max_examples=1, deadline=None)
    def test_order(self, h):
        self.assertLessEqual(h.order(), h.source_triangulation.max_order())
        self.assertEqual(h**(h.order()), h.source_triangulation.id_encoding())

