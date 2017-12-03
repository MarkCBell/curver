
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

from test_triangulation import triangulations

import curver

@st.composite
def curves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    edge = draw(st.sampled_from(triangulation.edges))
    path = []
    seen = set()
    while edge not in seen:
        seen.add(edge)
        path.append(edge)
        edge = ~draw(st.sampled_from(triangulation.corner_lookup[edge.label].edges[1:]))
    start = path.index(edge)
    return triangulation.curve_from_cut_sequence(path[start:])

@st.composite
def multicurves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    curves = draw(st.lists(elements=curves(triangulation), min_size=1))
    return triangulation.sum(curves)


class TestCurve(unittest.TestCase):
    pass

