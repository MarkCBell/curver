
''' A module for representing more advanced ways of changing triangulations.

Provides: Crush, Lift. '''

import curver
from curver.kernel.moves import Move  # Special import needed for subclassing.

class Crush(Move):
	''' This represents the effect of crushing along a curve. '''
	def __init__(self, source_triangulation, target_triangulation, curve, matrix):
		super(Crush, self).__init__(source_triangulation, target_triangulation)
		
		assert(isinstance(curve, curver.kernel.Curve))
		assert(curve.is_short())
		assert(curve.triangulation == self.source_triangulation)
		
		self.curve = curve
		self.matrix = matrix
	
	def __str__(self):
		return 'Crush ' + str(self.curve)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve, self.matrix))
	
	def apply_lamination(self, lamination):
		geometric = list(lamination)
		if self.curve.weight() == 2:  # Non-isolating.
			# Get some edges.
			a = self.curve.parallel()
			_, b, e = self.source_triangulation.corner_lookup[a.label]
			_, c, d = self.source_triangulation.corner_lookup[~e.label]
			
			u = self.source_triangulation.vertex_lookup[a.label]  # = self.triangulation.vertex_lookup[~a.label].
			u_edges = curver.kernel.utilities.cyclic_slice(u, a, ~a)
			around_u = min(max(lamination.side_weight(edge), 0) for edge in u_edges)  # The number of components that go through a, around u and then back through ~a.
			out_u = sum(max(-lamination.side_weight(edge), 0) for edge in u_edges) + sum(max(-lamination(edge), 0) for edge in u_edges[1:])  # The number of arcs that come out of u (from a around to ~a).
			# Since around_u > 0 ==> out_u == 0 & out_u > 0 ==> around_u == 0, this is equivalent to around_u if around_u > 0 else -out_u
			geometric[b.index] = around_u - out_u
			
			v = self.source_triangulation.vertex_lookup[c.label]  # = self.triangulation.vertex_lookup[~c.label].
			v_edges = curver.kernel.utilities.cyclic_slice(v, c, ~c)
			around_v = min(max(lamination.side_weight(edge), 0) for edge in v_edges)
			out_v = sum(max(-lamination.side_weight(edge), 0) for edge in v_edges) + sum(max(-lamination(edge), 0) for edge in v_edges[1:])
			geometric[e.index] = around_v - out_v  # Same trick.
		else:  # self.curve is isolating.
			# Get some edges.
			a = self.curve.parallel()
			v = self.curve.triangulation.vertex_lookup[a.label]  # = self.triangulation.vertex_lookup[~a.label].
			_, b, e = self.source_triangulation.corner_lookup[a.label]
			
			v_edges = curver.kernel.utilities.cyclic_slice(v, a, ~a)  # The set of edges that come out of v from a round to ~a.
			around_v = min(max(lamination.side_weight(edge), 0) for edge in v_edges)
			out_v = sum(max(-lamination.side_weight(edge), 0) for edge in v_edges) + sum(max(-lamination(edge), 0) for edge in v_edges[1:])
			# around_v > 0 ==> out_v == 0; out_v > 0 ==> around_v == 0.
			
			# TODO: 1) WRONG!!
			twisting = min(max(lamination.side_weight(edge) - around_v, 0) for edge in v_edges[1:-1])
			drops = [max(lamination.side_weight(edge) - (0 if index in (0, len(v_edges)) else twisting) - around_v, 0) for index, edge in enumerate(v_edges)]
			left_tighten = [min(drops[:i+1]) for i in range(len(v_edges))]
			right_tighten = [min(drops[i:]) for i in range(len(v_edges))]
			cumulative_drops = [max(min(drops[:i+1]), min(drops[i:])) for i in range(len(v_edges))]
			
			# Unwind.
			for edge, cumulative_drop in zip(v_edges[1:], cumulative_drops):
				geometric[edge.index] -= twisting + cumulative_drop
			
			# Now have to reset b weight.
			geometric[b.index] = around_v - out_v  # Same trick as above.
		
		return self.target_triangulation.lamination(geometric)  # Have to promote.
	
	def apply_homology(self, homology_class):
		return NotImplemented  # I don't think we ever need this.
	
	def inverse(self):
		return Lift(self.target_triangulation, self.source_triangulation, self.curve, self.matrix)

class Lift(Move):
	''' This represents the inverse of crushing along a curve. '''
	def __init__(self, source_triangulation, target_triangulation, curve, matrix):
		super(Lift, self).__init__(source_triangulation, target_triangulation)
		
		assert(isinstance(curve, curver.kernel.Curve))
		assert(curve.triangulation == self.target_triangulation)
		
		self.curve = curve
		self.matrix = matrix
	
	def __str__(self):
		return 'Lift ' + str(self.curve)
	def __reduce__(self):
		return (self.__class__, (self.source_triangulation, self.target_triangulation, self.curve, self.matrix))
	
	def apply_lamination(self, lamination):
		# Really should check that the dual weights around a vertex are all non-negative.
		geometric = [sum(x * y for x, y in zip(row, lamination)) for row in self.matrix]  # Dot product.
		return self.target_triangulation.lamination(geometric)  # Have to promote.
	
	def apply_homology(self, homology_class):
		return NotImplemented  # I don't think we ever need this.
	
	def inverse(self):
		return Crush(self.target_triangulation, self.source_triangulation, self.curve, self.matrix)

