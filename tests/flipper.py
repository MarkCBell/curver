
import pytest
import unittest

from hypothesis import given, assume
import hypothesis.strategies as st

import curver

flipper = pytest.importorskip('flipper')  # Skip these tests if flipper is not available.

class TestFlipper(unittest.TestCase):
    @staticmethod
    def mapping_classes(data):
        name = data.draw(st.sampled_from(['S_1_1', 'S_1_2', 'S_2_1']))
        S_curver = curver.load(name)
        S_flipper = flipper.load(name)
        word = ''.join(data.draw(st.lists(elements=st.sampled_from(sorted(S_curver.mapping_classes)))))
        return S_curver(word), S_flipper(word)
    
    @given(st.data())
    def test_nielsen_thurston_type(self, data):
        g, h = self.mapping_classes(data)
        self.assertEqual(g.nielsen_thurston_type(), h.nielsen_thurston_type())
    
    @given(st.data())
    def test_dilatation(self, data):
        g, h = self.mapping_classes(data)
        assume(g.is_pseudo_anosov())
        self.assertEqual(g.dilatation(), h.dilatation())
    
    @given(st.data())
    def test_stratum(self, data):
        g, h = self.mapping_classes(data)
        assume(g.is_pseudo_anosov())
        
        flipper_stratum = sorted(c - 2 for c in h.stratum().values())
        self.assertEqual(g.stratum(), flipper_stratum)
