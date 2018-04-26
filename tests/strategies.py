
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
def mapping_classes(draw, triangulation=None):
    return draw(encodings(triangulation, distribution=[2, 3]))

@st.composite
def mappings(draw, triangulation=None):
    return draw(encodings(triangulation, distribution=[0, 0, 0, 0, 1, 2, 3]))

@st.composite
def encodings(draw, triangulation=None, distribution=None):
    if triangulation is None: triangulation = draw(triangulations())
    if distribution is None: distribution = [0, 0, 0, 0, 1, 2, 3, 4]
    h = triangulation.id_encoding()
    move_types = draw(st.lists(elements=st.sampled_from(distribution), max_size=10))
    for move_type in move_types:
        T = h.target_triangulation
        if move_type == 0:  # EdgeFlip.
            edge = draw(st.sampled_from([edge for edge in T.edges if T.is_flippable(edge)]))
            move = T.encode_flip(edge)
        elif move_type == 1:  # Isometry.
            move = T.encode_relabel_edges([i if draw(st.booleans()) else ~i for i in draw(st.permutations(range(T.zeta)))])
        elif move_type == 2:  # Twist.
            curve = draw(curves(T))
            move = curve.encode_twist(power=draw(st.integers(min_value=-10, max_value=10).filter(lambda p: p)))
        elif move_type == 3:  # HalfTwist.
            arcs = [T.edge_arc(edge) for edge in T.edges if T.vertex_lookup[edge] != T.vertex_lookup[~edge]]
            if arcs:
                arc = draw(st.sampled_from(arcs))
                move = arc.encode_halftwist(power=draw(st.integers(min_value=-10, max_value=10).filter(lambda p: p)))
            else:
                move = T.id_encoding()
        else:  # move_type == 4:  # Crush.
            curve = draw(curves(T))
            move = curve.crush()
        
        h = move * h
    
    return h

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
    
    indices = set()
    available_indices = set(triangulation.indices)
    
    tree = triangulation.dual_tree()
    tree = set([index for index in triangulation.indices if tree[index]])
    available_indices = available_indices - tree
    
    classes = curver.kernel.UnionFind(triangulation.vertices)
    for component, (g, v) in triangulation.surface().items():
        available_component_indices = set([index for index in available_indices if index in component])
        if g != 0 and draw(st.booleans()):  # merge vertices:
            for index in sorted(available_component_indices):
                a, b = triangulation.vertex_lookup[index], triangulation.vertex_lookup[~index]
                if classes(a) != classes(b):
                    classes.union(a, b)
                    indices.add(index)
                    available_component_indices.remove(index)
        
        # Add in some of the remaining indices.
        indices.update(draw(st.sets(elements=st.sampled_from(sorted(available_component_indices)), min_size=0 if indices else 1, max_size=len(available_component_indices)-1)))
    
    # Set weights on these edges.
    geometric = [-1 if index in indices else 0 for index in triangulation.indices]
    
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
    
    @given(mappings())
    def test_mappings(self, h):
        self.assertIsInstance(h, curver.kernel.Mapping)
    
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

