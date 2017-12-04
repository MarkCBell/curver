
from hypothesis import given, settings
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver

# TRIANGULATIONS = {(g, p): curver.load(g, p).triangulation.sig() for g in range(10) for p in range(1, 10) if 6*g + 3*p - 6 >= 3}
TRIANGULATIONS = {
    (0, 3): '3_p',
    (0, 4): '6_JDky1',
    (0, 5): '9_ZFO4OuIV',
    (0, 6): 'c_ZlgeM906o-354',
    (0, 7): 'f_ZBhisze7ZhPen+Kdl1',
    (0, 8): 'i_ZBCol3oodJti4fabgfrLdu1',
    (0, 9): 'l_ZBmFdZ5TmXJTHbErcDHjdC12vbv4',
    (1, 1): '3_o',
    (1, 2): '6_WKSv',
    (1, 3): '9_q3AtQIBp',
    (1, 4): 'c_qfj5RWuFLb802',
    (1, 5): 'f_qLjT1kC5ZTnGyy1-H',
    (1, 6): 'i_qLbeF0856iXdcySItPVI7O',
    (1, 7): 'l_qLH58G-B-C5woL+MsXnS96SB0Ur2',
    (1, 8): 'o_qLH3BJ5cZdEYqFJK8US2y-FobflUMuJXh',
    (1, 9): 'r_qLHzvY4m6fubS+XBc6Fe3q2oydqacuUFP8knNo4',
    (2, 1): '9_qBFcQEN4',
    (2, 2): 'c_q5YHALgCWTyw',
    (2, 3): 'f_q5sG87f308pRczx4k',
    (2, 4): 'i_q5s8zbvSopJhvzCSeTFSHs',
    (2, 5): 'l_q5sEtshJlZcP3VxnSg50en+0fYA1',
    (2, 6): 'o_q5sEjz8JquEV7VEKDk7V1jA2MXn78nNzc',
    (2, 7): 'r_q5sEXh18WxhR6p+c1t3FqBE5ADBr5HuhrOHDBf3',
    (2, 8): 'u_q5sEXJrSMGIad7plECHCF-ZDpcWLW83ZkNe4eiHHmFlD1',
    (2, 9): 'x_q5sEXd3XbSQm+-kXyrXCVYKYFcTz6Z5jHLHM+udBWiyxoM5Zbt1',
    (3, 1): 'f_0mGaGVParMxLYnCS3',
    (3, 2): 'i_0mGG3b5sZ48hMGl+XMxQf7',
    (3, 3): 'l_0mGGHziLfC85djn-maLAMR4st0K',
    (3, 4): 'o_0mGGHv87rZuOX6acFcu2Pxuu+DqmvfRb7',
    (3, 5): 'r_0mGGH-frv5VxbMIDZJg+PfsvsGkD+T4V2s3Up62',
    (3, 6): 'u_0mGGH-P8pozBmiNXbK0fZpo0ncTEiu04QdNjBn9ZJgn81',
    (3, 7): 'x_0mGGH-3ktQlFDgcyM9fkSxjeT1cb5E-PwJDw3I-lkIOMDUAcV41',
    (3, 8): 'A_0mGGH-3rGANCo9JE6TmZaKNI6kuZUobdsYwLbLayo8Ub3hD75-FeycbL1',
    (3, 9): 'D_0mGGH-3XUmszRC-f0Kyqkw4UM8-vbAhPgaq-V4Gb6SjpdRrqXpP6MLNexfnVuA4',
    (4, 1): 'l_0wtlMpgmgKA+k2ObfDEOlf+tdih1',
    (4, 2): 'o_0wtlMFe+inn8rzo8rBgSeeSJdPKrw1OFa',
    (4, 3): 'r_0wtlMFCguCalYSPETZQ9G4DvnpCtAr6eMZFLPR2',
    (4, 4): 'u_0wtlMF+WbZ81Ll+SKPg5sOpjSoG6gQWjbTLQZ38oPavs1',
    (4, 5): 'x_0wtlMF+eOhAEjOyEX96St+C4dYcwSbwozSel1VJ0du3wvIPKJk1',
    (4, 6): 'A_0wtlMF+uWL8zHWMw5In5jPZGKqe3qImPInE7fKCNV2uUrxwodSD9H7X42',
    (4, 7): 'D_0wtlMF+uW2PslyYLdtjcGMFYloqw8elhB75LSzd7FlORRH0vo8y8+7Bsxnr0lm5',
    (4, 8): 'G_0wtlMF+uKiqZQ4wW85xwYGeuWoKctihKeCaJXLhpKtRsAj00sS960b23XP0t67ospXIOl',
    (4, 9): 'J_0wtlMF+uKu-uOzqzZ+7kU2Lq4I3BEeaaVa57q2C0lhFS7iV1YBCH7IbI-L3HePt7iIGBe5WgA282',
    (5, 1): 'r_0wUXPZwUD9uA4K-fmxT93jrye5mW+6QUPZVGQU1',
    (5, 2): 'u_0wUXPZw7sxaMvf+HM6p0m50javFCeVyOYYsOtrBCvcH21',
    (5, 3): 'x_0wUXPZwTnJKfEIppmxHddHc+ER8yCqQ445v0d2Q101uZCcV6u01',
    (5, 4): 'A_0wUXPZwTTUvsjvktbsTgIgjJ7aqCDyJNKbtky0Ajvrz4SWEQ+5CjlC9F1',
    (5, 5): 'D_0wUXPZwTHMJ9w2zmjkrNnNtK7FwyVbF-mL4UHfta7eN54N2lZekGdxwvKMt1Om4',
    (5, 6): 'G_0wUXPZwTHYiyz+sdTFf3khZz195gt99PeEiX+MLv9mUa4O4yfOqo3b--xuvog9dP2fO9i',
    (5, 7): 'J_0wUXPZwTH0cWeOU-yEvWsxsiJOElP-Vml7Bi2JyY24itGFBv1Z0bif-hsK1c9zL-QETRIcKg9pP1',
    (5, 8): 'M_0wUXPZwTH0cRiVJPBRoeCnlyjorpiQI8coG+UQtINlFoP-I6sOedoZGzA8cPUK3lFgQIaV7vSIaoDajr1h',
    (5, 9): 'P_0wUXPZwTH0sVLX2dmkM054OZAYTbfU7OhNURAF5-lmsXhmepRzsMjqrme2wP+y2oFw16Ud-04jQ7uCBbFXMocr5F3',
    }

SIGNATURES = [TRIANGULATIONS[key] for key in sorted(TRIANGULATIONS)]

@st.composite
def triangulations(draw):
    sig = draw(st.sampled_from(SIGNATURES))
    return curver.triangulation_from_sig(sig)

MCGS = [curver.load(g, p) for g in range(0, 5) for p in range(1, 5) if 6*g + 3*p - 6 >= 3]

@st.composite
def mcgs(draw):
    return draw(st.sampled_from(MCGS))

@st.composite
def mapping_classes(draw, mcg=None):
    if mcg is None: mcg = draw(mcgs())
    word = draw(st.lists(elements=st.sampled_from(sorted(mcg)), max_size=10).map('.'.join))
    return mcg(word)

@st.composite
def encodings(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    rev_sequence = [triangulation.id_isometry()]
    num_flips = draw(st.integers(min_value=0, max_value=20))
    for _ in range(num_flips):
        T = rev_sequence[-1].target_triangulation
        edge = draw(st.sampled_from([edge for edge in T.edges if T.is_flippable(edge)]))
        flip = curver.kernel.EdgeFlip(T, T.flip_edge(edge), edge)
        rev_sequence.append(flip)
    
    return curver.kernel.Encoding(rev_sequence[::-1])

@st.composite
def homology_classes(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    algebraic = [draw(st.integers()) for _ in range(triangulation.zeta)]
    return curver.kernel.HomologyClass(triangulation, algebraic)



@st.composite
def arcs(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    
    # h = draw(encodings(triangulation))
    h = triangulation.id_encoding()
    edge = draw(st.sampled_from(triangulation.edges))
    arc = triangulation.edge_arc(edge)
    
    return (~h)(arc)

@st.composite
def multiarcs(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    
    # h = draw(encodings(triangulation))
    h = triangulation.id_encoding()
    geometric = draw(st.lists(elements=st.integers(max_value=0), min_size=h.target_triangulation.zeta, max_size=h.target_triangulation.zeta))
    multiarc = h.target_triangulation.lamination(geometric)
    
    return (~h)(multiarc)

# TODO: 1) Non-isolating curves are *very* rare in this strategy.
@st.composite
def curves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    
    # h = draw(encodings(triangulation))
    h = triangulation.id_encoding()
    
    edge = draw(st.sampled_from(h.target_triangulation.edges))
    path = []
    seen = set()
    while edge not in seen:
        seen.add(edge)
        path.append(edge)
        edge = ~draw(st.sampled_from(h.target_triangulation.corner_lookup[edge.label].edges[1:]))
    start = path.index(edge)
    multicurve = h.target_triangulation.lamination_from_cut_sequence(path[start:])
    curve = multicurve.peek_component()
    
    return (~h)(curve)

@st.composite
def multicurves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    pieces = draw(st.lists(elements=st.tuples(curves(triangulation), st.integers(min_value=1, max_value=100)).map(lambda pair: pair[0]*pair[1])))
    return triangulation.sum(pieces)

@st.composite
def laminations_old(draw, triangulation=None, min_weight=None, max_weight=None):
    if triangulation is None: triangulation = draw(triangulations())
    geometric = [None] * triangulation.zeta
    for index in triangulation.indices:
        _, a, b, _, c, d = triangulation.corner_lookup[index].indices + triangulation.corner_lookup[~index].indices
        a, b, c, d = [max(geometric[i], 0) if geometric[i] is not None else None for i in [a, b, c, d]]
        if all(w is not None for w in [a, b, c, d]):
            lower_ab, upper_ab = abs(a - b), a + b
            parity_ab = (a + b) % 2
            lower_cd, upper_cd = abs(c - d), c + d
            parity_cd = (c + d) % 2
            exclude_lower, exclude_upper = (max(lower_ab, lower_cd), min(upper_ab, upper_cd)) if parity_ab != parity_cd else (0, 0)
            e = draw(st.one_of(st.integers(max_value=exclude_lower), st.integers(min_value=exclude_upper)).map(
                lambda x: x + 1 if (x % 2 != parity_ab and lower_ab < x < upper_ab) or (x % 2 != parity_cd and lower_cd < x < upper_cd) else x
                ))
        elif all(w is not None for w in [a, b]):
            lower, upper = abs(a - b), a + b
            parity = (a + b) % 2
            e = draw(st.integers().map(lambda x: x + 1 if x % 2 != parity and lower < x < upper else x))
        elif all(w is not None for w in [c, d]):
            lower, upper = abs(c - d), c + d
            parity = (c + d) % 2
            e = draw(st.integers().map(lambda x: x + 1 if x % 2 != parity and lower < x < upper else x))
        else:
            e = draw(st.integers())
        geometric[index] = e
    return triangulation.lamination(geometric)

@st.composite
def laminations(draw, triangulation=None):
    return NotImplemented
    return draw(st.one_of(multicurve(triangulation), multiarc(triangulation)))

@st.composite
def permutations(draw, N=None):
    if N is None: N = draw(st.integers(min_value=1, max_value=10))
    return curver.kernel.Permutation(draw(st.permutations(range(N))))


class TestStrategiesHealth(unittest.TestCase):
    @given(triangulations())
    def test_triangulations(self, triangulation):
        self.assertIsInstance(triangulation, curver.kernel.Triangulation)

    @given(mcgs())
    def test_mcgs(self, mcg):
        self.assertIsInstance(mcg, curver.kernel.MappingClassGroup)

    @given(mapping_classes())
    def test_mapping_classes(self, h):
        self.assertIsInstance(h, curver.kernel.MappingClass)

    @given(encodings())
    def test_encodings(self, encoding):
        self.assertIsInstance(encoding, curver.kernel.Encoding)

    @given(homology_classes())
    def test_homology_classes(self, hc):
        self.assertIsInstance(hc, curver.kernel.HomologyClass)

    @given(arcs())
    def test_arcs(self, arc):
        self.assertIsInstance(arc, curver.kernel.Arc)

    @given(multiarcs())
    def test_multiarcs(self, multiarc):
        self.assertIsInstance(multiarc, curver.kernel.MultiArc)

    @given(curves())
    def test_curves(self, curve):
        self.assertIsInstance(curve, curver.kernel.Curve)

    @given(multicurves())
    def test_multicurves(self, multicurve):
        self.assertIsInstance(multicurve, curver.kernel.MultiCurve)

    @given(laminations_old())
    def test_laminations_old(self, lamination):
        self.assertIsInstance(lamination, curver.kernel.Lamination)

    @pytest.mark.skip('NotImplemented')
    @given(laminations())
    def test_laminations(self, lamination):
        self.assertIsInstance(lamination, curver.kernel.Lamination)

    @given(permutations())
    def test_permutations(self, perm):
        self.assertIsInstance(perm, curver.kernel.Permutation)

if __name__ == '__main__':
    unittest.main()

