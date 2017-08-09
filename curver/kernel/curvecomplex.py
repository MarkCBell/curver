
class CurveComplex(object):
	def __init__(self, triangulation):
		self.triangulation = triangulation
		# Put more constants here.
	
	def quasiconvex(self, a, b):
		''' Return a polynomial-sized K--quasiconvex subset of the curve complex that contains a and b. '''
		
		assert(isinstance(a, Curve))
		assert(isinstance(b, Curve))
		assert(a.triangulation == self.triangulation)
		assert(b.triangulation == self.triangulation)
		
		short, conjugator = a.shorten()
		
		train_track = conjugator(b).train_track()
		_, conjugator_tt = train_track.shorten()
		encodings = [conjugator_tt[i:] for i in range(len(conjugator_tt)+1)]
		return set([conjugator.inverse()(encoding.inverse()(encoding(train_track).vertex_cycle())) for encoding in encodings])
	
	def tight_paths(self, a, b, length):
		''' Return the set of all tight paths from a to b that are of the given length.
		
		From Algorithm 3 of [BellWebb16b]. '''
		assert(isinstance(a, MultiCurve))
		assert(isinstance(b, MultiCurve))
		assert(a.triangulation == self.triangulation)
		assert(b.triangulation == self.triangulation)
		assert(length >= 0)
		
		if length == 0:
			return set([(a,)]) if a == b else set()
		elif length == 1:
			return set([(a, b)]) if a.intersection(b) == 0 and a.no_common_component(b) else set()
		elif length == 2:
			m = a.boundary_union(b)  # m = \partial N(a \cup b).
			return set([(a, m, b)]) if not m.is_empty() and a.no_common_component(m) and b.no_common_component(m) else set()
		else:  # length >= 3.
			crush = a.crush()
			lift = crush.inverse()
			b_prime = crush(b)
			A_1 = set()
			for triangulation in b_prime.eplore_ball(2*self.zeta*length + 2*self.zeta):
				for submultiarc in triangulation.sublaminations():
					m_prime = submultiarc.boundary()
					m = lift(m_prime)
					A_1.add(m)
			
			P = set()
			for a_1 in A_1:
				for multipath in a_1.tight_paths(b, length-1):  # Recurse.
					a_2 = multipath[0]
					if a.boundary_union(a_2) == a_1:  # (a,) + multipath is tight:
						P.add((a,) + multipath)
			
			return P
	
	def all_tight_geodesic_multicurves(self, a, b):
		''' Return a set that contains all multicurves in any tight geodesic from a to b.
		
		From the first half of Algorithm 4 of [BellWebb16b]. '''
		
		assert(isinstance(a, Curve))
		assert(isinstance(b, Curve))
		assert(a.triangulation == self.triangulation)
		assert(b.triangulation == self.triangulation)
		
		guide = self.quasiconvex(a, b)  # U.
		L = 6*curver.kernel.constants.QUASICONVEXITY + 2  # See [Webb15].  !?!
		return set(multicurve for length in range(L+1) for c1 in guide for c2 in guide for path in self.tight_paths(c1, c2, length) for multicurve in path)
	
	def tight_geodesic(self, a, b):
		''' Return a tight geodesic in the (multi)curve complex from a to b.
		
		From the second half of Algorithm 4 of [BellWebb16b]. '''
		
		assert(isinstance(a, Curve))
		assert(isinstance(b, Curve))
		assert(a.triangulation == self.triangulation)
		assert(b.triangulation == self.triangulation)
		
		vertices = list(self.all_tight_geodesic_multicurves(a, b))
		edges = [(i, j) for i in range(len(vertices)) for j in range(i) if vertices[i].intersection(vertices[j]) == 0 and vertices[i].no_common_component(vertices[j])]
		
		G = networkx.Graph(edges)  # Build graph.
		indices = networkx.shortest_path(G, vertices.index(a), vertices.index(b))  # Find a geodesic from self to other.
		geodesic = [vertices[index] for index in indices]  # Get the geodesic, however this might not be tight.
		
		for i in range(1, len(geodesic)-1):
			geodesic[i] = geodesic[i-1].boundary_union(geodesic[i+1])  # Tighten.
		
		return tuple(geodesic)
	
	def geodesic(self, a, b):
		''' Return a geodesic in the curve complex from a to b.
		
		The geodesic will always come from a tight geodesic.
		From Algorithm 5 of [BellWebb16b]. '''
		
		assert(isinstance(a, Curve))
		assert(isinstance(b, Curve))
		assert(a.triangulation == self.triangulation)
		assert(b.triangulation == self.triangulation)
		
		return tuple(multicurve.peek_component() for multicurve in self.tight_geodesic(a, b))
	
	def distance(self, a, b):
		''' Return the distance from a to b in the curve complex. '''
		
		# Could use self.tight_geodesic(a, b).
		return len(self.geodesic(a, b)) - 1

