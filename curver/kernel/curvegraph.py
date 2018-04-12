
''' A module for representing the curve complex of a surface. '''

from itertools import combinations
from math import factorial
import networkx

import curver

class CurveGraph(object):
    ''' This represents the curve complex of a surface. '''
    def __init__(self, triangulation):
        self.triangulation = triangulation
        self.zeta = self.triangulation.zeta
        assert self.triangulation.is_connected()
        
        # Constants:
        # Universal:
        self.QUASICONVEXITY = 10  # Masur-Minsky.  # Ask Katie Vokes about this.
        self.BOUNDED_GEODESIC_IMAGE = 100  # BGIT.
        self.HYPERBOLICITY = 17  # Curve complex delta-hyperbolicity.
        # Uniform:
        self.D = (self.zeta**(2 * self.zeta * (2*self.HYPERBOLICITY + 2)**2))**2  # Bowditch bound on denominator [Bow08].
        # This explict bound comes from Theorems 6.2 & 6.3 of [Webb15] and Algorithm 3 of [BellWebb16].
        # TODO 4) Recheck this.
        self.XI = self.zeta // 2  # Actually should be 3g - 3 + n < 3g - 3 + 1.5n = zeta / 2
        self.M = 28 * factorial(self.XI) * self.HYPERBOLICITY * self.D * self.BOUNDED_GEODESIC_IMAGE
        E, n = self.triangulation.euler_characteristic, self.triangulation.num_vertices
        self.R = 18*E**2 - 30*E - 10*n  # Gadre and Tsai bound away from zero [GadreTsai].
        self.M2 = factorial(self.XI) * self.BOUNDED_GEODESIC_IMAGE * self.R
    
    def quasiconvex(self, a, b):
        ''' Return a polynomial-sized K--quasiconvex subset of the curve complex that contains a and b.
        
        From Algorithm 1 of [BellWebb16]_. '''
        
        # Uses Masur--Minsky theorem.
        
        assert isinstance(a, curver.kernel.Curve)
        assert isinstance(b, curver.kernel.Curve)
        assert a.triangulation == self.triangulation
        assert b.triangulation == self.triangulation
        
        conjugator_a = a.shorten()
        
        short_b = conjugator_a(b)
        conjugator_b = short_b.shorten()
        
        quasiconvex = set()
        for i in range(len(conjugator_b)+1):
            prefix = conjugator_b[i:]  # Get the first bit of conjugator_b.
            split_train_track = prefix(short_b)
            vertex_cycle = split_train_track.vertex_cycle()
            quasiconvex.add(conjugator_a.inverse()(prefix.inverse()(vertex_cycle)))
        
        return quasiconvex
    
    def tight_paths(self, a, b, length):
        ''' Return the set of all tight paths from a to b that are of the given length.
        
        From Algorithm 3 of [BellWebb16]_. '''
        
        assert isinstance(a, curver.kernel.MultiCurve)
        assert isinstance(b, curver.kernel.MultiCurve)
        assert a.triangulation == self.triangulation
        assert b.triangulation == self.triangulation
        assert length >= 0
        
        if length == 0:
            return set([(a,)]) if a == b else set()
        elif length == 1:
            return set([(a, b)]) if a.intersection(b) == 0 and a.no_common_component(b) else set()
        elif length == 2:
            m = a.boundary_union(b)  # m = \partial N(a \cup b).
            return set([(a, m, b)]) if not m.is_peripheral() and a.no_common_component(m) and b.no_common_component(m) else set()
        else:  # length >= 3.
            if not a.fills_with(b): return set()
            crush = a.crush()
            lift = crush.inverse()
            b_prime = crush(b)
            A_1 = set()
            for triangulation in b_prime.explore_ball(2*self.zeta*length + 2*self.zeta):
                for submultiarc in triangulation.sublaminations():
                    m_prime = submultiarc.boundary()
                    m = lift(m_prime)
                    A_1.add(m)
            
            P = set()
            for a_1 in A_1:
                for multipath in self.tight_paths(a_1, b, length-1):  # Recurse.
                    a_2 = multipath[0]
                    if a.boundary_union(a_2) == a_1:  # (a,) + multipath is tight:
                        P.add((a,) + multipath)
            
            return P
    
    def all_tight_geodesic_multicurves(self, a, b):
        ''' Return a set that contains all multicurves in any tight geodesic from a to b.
        
        From the first half of Algorithm 4 of [BellWebb16]_. '''
        
        assert isinstance(a, curver.kernel.Curve)
        assert isinstance(b, curver.kernel.Curve)
        assert a.triangulation == self.triangulation
        assert b.triangulation == self.triangulation
        
        guide = self.quasiconvex(a, b)  # U.
        L = 6*self.QUASICONVEXITY + 2  # See [Webb15].
        return set(multicurve for length in range(L+1) for c1 in guide for c2 in guide for path in self.tight_paths(c1, c2, length) for multicurve in path)
    
    def tight_geodesic(self, a, b):
        ''' Return a tight geodesic in the (multi)curve complex from a to b.
        
        From the second half of Algorithm 4 of [BellWebb16]_. '''
        
        assert isinstance(a, curver.kernel.Curve)
        assert isinstance(b, curver.kernel.Curve)
        assert a.triangulation == self.triangulation
        assert b.triangulation == self.triangulation
        
        # Build graph.
        vertices = list(self.all_tight_geodesic_multicurves(a, b))
        edges = [(u, v) for u, v in combinations(vertices, r=2) if u.intersection(v) == 0 and u.no_common_component(v)]
        G = networkx.Graph(edges)
        
        geodesic = networkx.shortest_path(G, a, b)  # Find a geodesic from self to other, however this might not be tight.
        
        for i in range(1, len(geodesic)-1):
            geodesic[i] = geodesic[i-1].boundary_union(geodesic[i+1])  # Tighten.
        
        return tuple(geodesic)
    
    def geodesic(self, a, b):
        ''' Return a geodesic in the curve complex from a to b.
        
        The geodesic will always come from a tight geodesic.
        From Algorithm 5 of [BellWebb16]_. '''
        
        assert isinstance(a, curver.kernel.Curve)
        assert isinstance(b, curver.kernel.Curve)
        assert a.triangulation == self.triangulation
        assert b.triangulation == self.triangulation
        
        return tuple(multicurve.peek_component() for multicurve in self.tight_geodesic(a, b))
    
    def distance(self, a, b):
        ''' Return the distance from a to b in the curve complex. '''
        
        # Could use self.tight_geodesic(a, b).
        return len(self.geodesic(a, b)) - 1

