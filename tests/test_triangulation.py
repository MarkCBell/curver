
from hypothesis import given
import hypothesis.strategies as st
import pickle
import pytest
import unittest

import curver

# TRIANGULATIONS = {(g, p): str(curver.load(g, p).triangulation) for g in range(10) for p in range(1, 10) if 6*g + 3*p - 6 >= 3}
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
    (1, 3): '9_q3AtQIBp', (4, 8): 'G_0wtlMF+uKiqZQ4wW85xwYGeuWoKctihKeCaJXLhpKtRsAj00sS960b23XP0t67ospXIOl',
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
    (4, 9): 'J_0wtlMF+uKu-uOzqzZ+7kU2Lq4I3BEeaaVa57q2C0lhFS7iV1YBCH7IbI-L3HePt7iIGBe5WgA282',
    (5, 1): 'r_0wUXPZwUD9uA4K-fmxT93jrye5mW+6QUPZVGQU1',
    (5, 2): 'u_0wUXPZw7sxaMvf+HM6p0m50javFCeVyOYYsOtrBCvcH21',
    }

SIGNATURES = [TRIANGULATIONS[key] for key in sorted(TRIANGULATIONS)]

@st.composite
def triangulations(draw):
    sig = draw(st.sampled_from(SIGNATURES)
    return curver.triangulation_from_sig(sig)


class TestTriangulation(unittest.TestCase):
    @given(triangulations())
    def test_pickle(self, triangulation):
        strn = pickle.dumps(triangulation)
        self.assertEqual(triangulation, pickle.loads(strn))
    
    @given(st.data())
    def test_hash(self, data):
        triangulation1 = data.draw(triangulations())
        triangulation2 = data.draw(triangulations())
        self.assertTrue(hash(triangulation1) != hash(triangulation2) or triangulation1 == triangulation2)
    
    @given(st.data())
    def test_flip(self, data):
        triangulation = data.draw(triangulations())
        edge = data.draw(st.sampled_from(triangulation.edges))
        self.assertEqual(triangulation.surface(), triangulation.flip_edge(edge).surface())
        
    @given(triangulations())
    def test_isometry(self, triangulation):
        identity = triangulation.id_isometry()
        self.assertIn(identity, triangulation.self_isometries())
        self.assertTrue(triangulation.is_isometric_to(triangulation))

    @given(triangulations())
    def test_sig(self, triangulation):
        self.assertEqual(triangulation, curver.triangulation_from_sig(triangulation.sig()))

