
''' A module for representing train tracks on triangulations.

Provides: TrainTrack. '''

import networkx

import curver
from curver.kernel.lamination import Shortenable  # Special import needed for subclassing.

class TrainTrack(Shortenable):
	''' A Lamination in which each triangle is tripod free. '''
	
	def is_short(self):
		return all(self.shorten_score(edge) == 0 for edge in self.triangulation.edges)
	
	def shorten_score(self, edge):
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		if not self.triangulation.is_flippable(edge):
			return 0
		
		a, b, c, d, e = self.triangulation.square(edge)
		da, db, dc, dd, de = [self.dual_weight(edgy) for edgy in self.triangulation.square(edge)]
		if self(e) <= 0:
			return 0
		if b == ~d and da > 0 and db == 0 and de == 0:
			return 0
		if a == ~c and da == 0 and db > 0 and de == 0:
			return 0
		if de > 0 or da < 0 or db < 0:
			return 1
		
		if (da == 0 and db == 0) or (da == 0 and de == 0) or (db == 0 and de == 0):
			if (da > 0 and dc > 0) or (db > 0 and dd > 0):
				return 2
			else:
				return 3
		
		return 4
	
	def mcomponents(self):
		''' Return a set of pairs (component, multiplicity). '''
		
		short, conjugator = self.shorten()
		
		components = []
		for edge in short.triangulation.positive_edges:  # Only need to check half of them.
			# Check for an Arc here.
			if short(edge) < 0:
				geometric = [-1 if index == edge.index else 0 for index in short.triangulation.indices]
				component, multiplicity = curver.kernel.Arc(short.triangulation, geometric), abs(short(edge))
				components.append((conjugator.inverse()(component), multiplicity))  # Map it back onto self.
			
			# Check for a curve here.
			if short.triangulation.is_flippable(edge):
				a, b, c, d, e = short.triangulation.square(edge)
				da, db, dc, dd, de = [short.dual_weight(edgy) for edgy in short.triangulation.square(edge)]
				if b == ~d and da > 0 and db == 0 and de == 0:
					geometric = [1 if index == b.index or index == e.index else 0 for index in short.triangulation.indices]
					component, multiplicity = curver.kernel.Curve(short.triangulation, geometric), short(e)
					components.append((conjugator.inverse()(component), multiplicity))  # Map it back onto self.
		
		return components
	
	def vertex_cycles(self):
		''' Return the set of vertex cycles of this train track. '''
		
		def connected_to(edge):
			''' Yield the edges you can reach by travelling out of the given edge. '''
			corner = self.triangulation.corner_lookup[edge.label]
			if self.dual_weight(corner[1]): yield ~corner[2]
			if self.dual_weight(corner[2]): yield ~corner[1]
		
		# Build graph.
		edges = [(edge, edgy) for edge in self.triangulation.edges for edgy in connected_to(edge)]
		G = networkx.DiGraph(edges)
		
		cycles = []
		for cycle in networkx.simple_cycles(G):
			geometric = [0] * self.zeta
			for edge in cycle:
				geometric[edge.index] += 1
			cycles.append(curver.kernel.Curve(self.triangulation, geometric))
		
		return cycles

