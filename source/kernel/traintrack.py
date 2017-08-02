
import curver
from curver.kernel.lamination import Lamination  # Special import needed for subclassing.

INFTY = float('inf')

class TrainTrack(Lamination):
	''' A Lamination in which each triangle is tripod free. '''
	
	def split(self, edge):
		# Splittable ==> flippable.
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		assert(self.dual_weight(edge.label) <= 0 or self.dual_weight(~edge.label) <= 0)
		
		return self.triangulation.encode_flip(edge)
	
	def score(self, edge):
		if isinstance(edge, curver.IntegerType): edge = self.triangulation.edge_lookup[edge]  # If given an integer instead.
		
		score = 0
		corner = self.triangulation.corner_lookup[edge.label]
		if self.dual_weight(corner[0]) > 0: score += 1
		if self.dual_weight(corner[1]) <= 0: score += 1
		if self.dual_weight(corner[2]) <= 0: score += 1
		if self(edge) <= 0: score += 1
		return score
	
	def mcomponents(self):
		
		train_track = self
		
		def short_curve(L, edge):
			a, b, c, d, e = L.triangulation.square(edge)
			geometric = [1 if i == b.index or i == e.index else 0 for i in range(self.zeta)]
			if b == ~d and \
				L.dual_weight(a) >= 0 and \
				L.dual_weight(b) == 0 and \
				L.dual_weight(c) >= 0 and \
				L.dual_weight(d) == 0 and \
				L.dual_weight(e) == 0 and \
				L.dual_weight(~e) == 0:
				multiplicity = L(e)
			else:
				multiplicity = 0
			return geometric, multiplicity
		
		components = []
		encoding = train_track.triangulation.id_encoding()
		
		# Remove all the obvious arcs. This reduces the number of places we have to look for new arcs later.
		for index in train_track.triangulation.indices:
			if train_track(index) < 0:
				geometric = [0 if i != index else -1 for i in range(train_track.zeta)]
				component, multiplicity = curver.kernel.Arc(train_track.triangulation, geometric), abs(train_track(index))
				components.append((encoding.inverse()(component), multiplicity))
				train_track = train_track - multiplicity * component
		
		# Remove all the obvious curves. This reduces the number of places we have to look for new curves later.
		for index in train_track.triangulation.indices:
			if train_track.triangulation.is_flippable(index):
				geometric, multiplicity = short_curve(train_track, index)
				if multiplicity > 0:
					component = curver.kernel.Curve(train_track.triangulation, geometric)
					components.append((encoding.inverse()(component), multiplicity))
					train_track = train_track - multiplicity * component
		
		extra = []
		while not train_track.is_empty():
			# This edge is always splittable.
			to_split = min(extra + train_track.triangulation.edges, key=train_track.score)
			
			move = train_track.split(to_split)
			encoding = move * encoding
			train_track = move(train_track)
			
			# The only place where an obvious arc could appear is along the new edge we have just introduced.
			if train_track(to_split) < 0:
				geometric = [0 if i != to_split else -1 for i in range(train_track.zeta)]
				component, multiplicity = curver.kernel.Arc(train_track.triangulation, geometric), abs(train_track(to_split))
				components.append((encoding.inverse()(component), multiplicity))
				train_track = train_track - multiplicity * component
			
			# The only places where a short curve could appear is across an edge adjacent to the one we just flipped.
			for index in [edge.index for edge in train_track.triangulation.square(to_split)]:
				if train_track.triangulation.is_flippable(index):
					geometric, multiplicity = short_curve(train_track, index)
					if multiplicity > 0:
						component = curver.kernel.Curve(train_track.triangulation, geometric)
						components.append((encoding.inverse()(component), multiplicity))
						train_track = train_track - multiplicity * component
			
			# Accelerate!!
			a, b, c, d, e = train_track.triangulation.square(to_split)
			extra = [~a, ~d]
		
		return components
	
	def vertex_cycle(self):
		return NotImplemented  # TODO: 2).
	
	def splitting_sequence(self):
		''' Return a sequence of Encodings taking this train track to its simplest form. '''

