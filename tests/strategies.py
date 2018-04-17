
from hypothesis import given
import hypothesis.strategies as st
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
    (2, 1): '9_qBFcQEN4',
    (2, 2): 'c_q5YHALgCWTyw',
    (2, 3): 'f_q5sG87f308pRczx4k',
    (3, 1): 'f_0mGaGVParMxLYnCS3',
    (3, 2): 'i_0mGG3b5sZ48hMGl+XMxQf7',
    }

SIGNATURES = [TRIANGULATIONS[key] for key in sorted(TRIANGULATIONS)]

@st.composite
def triangulations(draw):
    sig = draw(st.sampled_from(SIGNATURES))
    return curver.triangulation_from_sig(sig)

MCGS = [curver.load(g, p) for g, p in [(0, 3), (0, 4), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2)]]

@st.composite
def mcgs(draw):
    return draw(st.sampled_from(MCGS))

@st.composite
def mapping_classes(draw, mcg=None):
    if mcg is None: mcg = draw(mcgs())
    return mcg(draw(st.integers(min_value=0, max_value=10)))

@st.composite
def encodings(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    rev_sequence = [triangulation.id_isometry()]
    num_flips = draw(st.integers(min_value=0, max_value=20))
    for _ in range(num_flips):
        T = rev_sequence[-1].target_triangulation
        if draw(st.sampled_from([0, 0, 0, 0, 0, 1])) == 0:
            edge = draw(st.sampled_from([edge for edge in T.edges if T.is_flippable(edge)]))
            move = T.encode_flip(edge)[0]
        else:
            move = T.encode_relabel_edges([i if draw(st.booleans()) else ~i for i in draw(st.permutations(range(T.zeta)))])[0]
        rev_sequence.append(move)
    
    return curver.kernel.Encoding(rev_sequence[::-1])

@st.composite
def homology_classes(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    algebraic = [draw(st.integers()) for _ in range(triangulation.zeta)]
    return curver.kernel.HomologyClass(triangulation, algebraic)

@st.composite
def multiarcs(draw, triangulation=None, require_non_empty_boundary=False):
    if triangulation is None: triangulation = draw(triangulations())
    
    geometric = [0] * triangulation.zeta
    
    available_indices = set(triangulation.indices)
    num_arcs = draw(st.integers(min_value=1, max_value=len(available_indices)))
    for _ in range(num_arcs):
        index = draw(st.sampled_from(sorted(available_indices)))
        available_indices.remove(index)
        geometric[index] = draw(st.integers(max_value=-1))
    
    return triangulation.lamination(geometric)

@st.composite
def arcs(draw, triangulation=None):
    return draw(multiarcs(triangulation)).peek_component()

@st.composite
def multicurves(draw, triangulation=None):
    if triangulation is None: triangulation = draw(triangulations())
    
    # Build a special type of multiarc with non-empty boundary.
    geometric = [0] * triangulation.zeta
    
    indices = set()
    available_indices = set(triangulation.indices)
    
    tree = triangulation.dual_tree()
    tree = set([index for index in triangulation.indices if tree[index]])
    available_indices = available_indices - tree
    
    if draw(st.booleans()) and not any(g == 0 for (g, v) in triangulation.surface()):  # merge vertices:
        classes = curver.kernel.UnionFind(triangulation.vertices)
        for index in sorted(available_indices):
            arc = triangulation.edge_arc(index)
            if arc.connects_distinct_vertices():
                a, b = arc.vertices()
                if classes(a) != classes(b):
                    classes.union(a, b)
                    indices.add(index)
                    available_indices.remove(index)
    
    num_arcs = draw(st.integers(min_value=0 if indices else 1, max_value=len(available_indices)-1))
    for _ in range(num_arcs):
        index = draw(st.sampled_from(sorted(available_indices)))
        available_indices.remove(index)
        indices.add(index)
    
    # Set weights on these edges.
    for index in indices:
        geometric[index] = draw(st.integers(max_value=-1))
    
    multiarc = triangulation.lamination(geometric)
    boundary = multiarc.boundary()
    assert(not boundary.is_empty())
    return boundary

@st.composite
def curves(draw, triangulation=None):
    multicurve = draw(multicurves(triangulation))
    curve = multicurve.peek_component()
    return curve

@st.composite
def laminations(draw, triangulation=None):
    return draw(st.one_of(multicurves(triangulation), multiarcs(triangulation)))

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
    
    @pytest.mark.skip('Incomplete')
    @given(laminations())
    def test_laminations(self, lamination):
        self.assertIsInstance(lamination, curver.kernel.Lamination)
    
    @given(permutations())
    def test_permutations(self, perm):
        self.assertIsInstance(perm, curver.kernel.Permutation)


if __name__ == '__main__':
    unittest.main()

