
import pytest
import unittest

from hypothesis import given
import hypothesis.strategies as st

import curver

flipper = pytest.importorskip('flipper')  # Skip these tests if flipper is not available.

class TestFlipper(unittest.TestCase):
    @given(st.data())
    def test_pseudo_anosov(self, data):
        name = data.draw(st.sampled_from(['S_1_1', 'S_1_2', 'S_2_1']))
        S = curver.load(name)
        SS = flipper.load(name)
        word = ''.join(data.draw(st.lists(elements=st.sampled_from(sorted(S.mapping_classes)))))
        g = S(word)
        h = SS(word)
        self.assertEqual(g.is_pseudo_anosov(), h.is_pseudo_anosov())
        
        if g.is_pseudo_anosov():
            self.assertEqual(g.dilatation(), h.dilatation())
            
            flipper_stratum = sorted(c - 2 for c in h.stratum().values())
            self.assertEqual(g.stratum(), flipper_stratum)

