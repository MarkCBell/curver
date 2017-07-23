''' A module for representing a triangulation of a punctured surface.

Provides five classes: Vertex, Edge, Triangle, Corner and Triangulation.
	A Vertex is a singleton.
	An Edge is an ordered pair of Vertices.
	A Triangle is an ordered triple of Edges.
	A Corner is a Triangle with a chosen side.
	A Triangulation is a collection of Triangles. '''

import curver

INFTY = float('inf')

def norm(value):
	''' A map taking an edges label to its index.
	
	That is, x and ~x should map to the same thing. '''
	
	return max(value, ~value)

class Edge(object):
	''' This represents an oriented edge, labelled with an integer.
	
	It is specified by its label and its inverse edge is labelled with ~its label.
	
	These are really just integers but with fancy printing and indexing set up on them. '''
	
	# Warning: This needs to be updated if the interals of this class ever change.
	__slots__ = ['label', 'index']
	
	def __init__(self, label):
		self.label = label
		self.index = norm(self.label)
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return ('' if self.sign() == +1 else '~') + str(self.index)
	def __reduce__(self):
		# Having __slots__ means we need to pickle manually.
		return (self.__class__, (self.label, None))
	def __hash__(self):
		return hash(self.label)
	def __eq__(self, other):
		if isinstance(other, Edge):
			return self.label == other.label
		else:
			return NotImplemented
	
	def __invert__(self):
		''' Return this edge but with reversed orientation. '''
		
		return Edge(~self.label)
	
	def sign(self):
		''' Return the sign (+/-1) of this edge. '''
		
		return +1 if self.label == self.index else -1

class Triangle(object):
	''' This represents a triangle.
	
	It is specified by a list of three edges, ordered anticlockwise.
	It builds its corners automatically. '''
	
	# Warning: This needs to be updated if the interals of this class ever change.
	__slots__ = ['edges', 'labels', 'indices', 'vertices', 'corners']
	
	def __init__(self, edges, rotate=None):
		assert(isinstance(edges, (list, tuple)))
		assert(all(isinstance(edge, Edge) for edge in edges))
		assert(len(edges) == 3)
		
		# Edges are ordered anti-clockwise. We will cyclically permute
		# these to a canonical ordering, the one where the edges are ordered
		# minimally by label.
		best_index = min(range(3), key=lambda i: edges[i].label) if rotate is None else rotate
		
		self.edges = edges[best_index:] + edges[:best_index]
		self.labels = [edge.label for edge in self]
		self.indices = [edge.index for edge in self]
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(tuple(self.edges))
	def __reduce__(self):
		# Having __slots__ means we need to pickle manually.
		return (self.__class__, (self.edges,))
	def __hash__(self):
		return hash(tuple(self.labels))
	def __eq__(self, other):
		if isinstance(other, Triangle):
			return self.labels == other.labels
		else:
			return NotImplemented
	def __len__(self):
		return 3  # This is needed for revered(triangle) to work.
	
	# Note that this is NOT the same convention as used in pieces.
	# There iterating and index accesses return vertices.
	def __iter__(self):
		return iter(self.edges)
	
	def __getitem__(self, index):
		return self.edges[index]

# Remark: In other places in the code you will often see L(triangulation). This is the space
# of laminations on triangulation with the coordinate system induced by the triangulation.

class Triangulation(object):
	''' This represents a triangulation of a punctured surface.
	
	It is specified by a list of Triangles. Its edges must be
	numbered 0, 1, ... and its vertices must be numbered 0, 1, ... '''
	def __init__(self, triangles):
		# We will sort the triangles into a canonical ordering, the one where the edges are ordered
		# minimally by label. This allows for fast comparisons.
		self.triangles = sorted(triangles, key=lambda t: t.labels)
		
		self.edges = [edge for triangle in self for edge in triangle.edges]
		self.positive_edges = [edge for edge in self.edges if edge.sign() == +1]
		self.labels = sorted([edge.label for edge in self.edges])
		self.indices = sorted([edge.index for edge in self.positive_edges])
		
		self.num_triangles = len(self.triangles)
		self.zeta = len(self.positive_edges)
		assert(self.zeta == self.num_triangles * 3 // 2)
		# Check that the edges have indices 0, ..., zeta-1.
		assert(set(self.labels) == set([i for i in range(self.zeta)] + [~i for i in range(self.zeta)]))
		
		self.edge_lookup = dict((edge.label, edge) for edge in self.edges)
		self.triangle_lookup = dict((edge.label, triangle) for triangle in self for edge in triangle.edges)
		self.corner_lookup = dict((edge.label, Triangle(triangle.edges, rotate=index)) for triangle in triangles for index, edge in enumerate(triangle))
		
		# Group the edges into vertex classes.
		# Here two edges are in the same class iff they have the same tail.
		unused = set(self.edges)
		self.vertex_classes = []
		while unused:
			new_vertex = [unused.pop()]
			while True:
				triangle = self.corner_lookup[new_vertex[-1].label]
				neighbour = ~triangle.edges[2]
				if neighbour in unused:
					new_vertex.append(neighbour)
					unused.remove(neighbour)
				else:
					break
			
			self.vertex_classes.append(tuple(new_vertex))
		
		self.vertex_lookup = dict((edge.label, vertex_class) for vertex_class in self.vertex_classes for edge in vertex_class)
		
		self.num_vertices = len(self.vertex_classes)
		
		self.euler_characteristic = self.num_vertices - self.zeta + self.num_triangles  # V - E + F.
		# NOTE: This assumes connected.
		self.genus = (2 - self.euler_characteristic) // 2
		
		# The maximum order of a periodic mapping class.
		# These bounds follow from the 4g + 4 bound on the closed surface [Primer reference]
		# and the Riemann removable singularity theorem which allows us to cap off the
		# punctures when the genus > 1 without affecting this bound.
		# NOTE: This assumes connected.
		if self.genus > 1:
			self.max_order = 4 * self.genus + 2
		elif self.genus == 1:
			self.max_order = max(self.num_vertices, 6)
		else:
			self.max_order = self.num_vertices
		
		# Two triangualtions are the same if and only if they have the same signature.
		self.signature = [e.label for t in self for e in t]
	
	@classmethod
	def from_tuple(cls, edge_labels):
		''' Return an Triangulation from a list of triples of edge labels.
		
		Let T be an ideal triangulaton of the punctured (oriented) surface S. Orient
		and edge e of T and assign an index i(e) in 0, ..., zeta-1. Now to each
		triangle t of T associate the triple j)t) := (j(e_1), j(e_2), j(e_3)) where:
			- e_1, e_2, e_3 are the edges of t, ordered acording to the orientation of t, and
			- j(e) = {  i(e) if the orientation of e agrees with that of t, and
					 { ~i(e) otherwise.
				Here ~x := -1 - x, the two's complement of x.
		
		We may describe T by the list [j(t) for t in T]. This function reconstructs
		T from such a list.
		
		edge_labels must be a list of triples of integers and each of
		0, ..., zeta-1, ~0, ..., ~(zeta-1) must occur exactly once. '''
		
		assert(isinstance(edge_labels, (list, tuple)))
		assert(all(isinstance(labels, (list, tuple)) for labels in edge_labels))
		assert(all(len(labels) == 3 for labels in edge_labels))
		assert(len(edge_labels) > 0)
		
		zeta = len(edge_labels) * 3 // 2
		
		# Check that each of 0, ..., zeta-1, ~0, ..., ~(zeta-1) occurs exactly once.
		flattened = set([label for labels in edge_labels for label in labels])
		for i in range(zeta):
			if i not in flattened:
				raise TypeError('Missing label %d' % i)
			if ~i not in flattened:
				raise TypeError('Missing label ~%d' % i)
		
		return cls([Triangle([Edge(label) for label in labels]) for labels in edge_labels])
	
	def __repr__(self):
		return str(self)
	def __str__(self):
		return str(list(self))
	def __iter__(self):
		return iter(self.triangles)
	def __getitem__(self, index):
		return self.triangles[index]
	def package(self):
		''' Return a small amount of info that create_triagulation can use to reconstruct this triangulation. '''
		return ([t.labels for t in self],)
	def __reduce__(self):
		# Triangulations are already pickleable but this results in a much smaller pickle.
		return (create_triangulation, (self.__class__,) + self.package())
	def __eq__(self, other):
		return self.signature == other.signature
	def __ne__(self, other):
		return not(self == other)
	
	def is_flippable(self, edge_label):
		''' Return if the given edge is flippable.
		
		An edge is flippable if and only if it lies in two distinct triangles. '''
		
		return self.triangle_lookup[edge_label] != self.triangle_lookup[~edge_label]
	
	def folded_boundary(self, edge_label):
		''' Return the edge bounding the once-punctured monogon containing edge_label.
		
		The given edge must not be flippable. '''
		
		assert(not self.is_flippable(edge_label))
		
		# As edge_label is not flippable the triangle containing it must be (edge_label, ~edge_label, boundary_edge).
		[boundary_edge] = [edge for edge in self.triangle_lookup[edge_label] if edge.label != edge_label and edge.label != ~edge_label]
		return boundary_edge
	
	def square_about_edge(self, edge_label):
		''' Return the four edges around the given edge.
		
		The chosen edge must be flippable. '''
		
		assert(self.is_flippable(edge_label))
		
		# Given the label e, return the edges a, b, c, d in order.
		#
		# #<----------#
		# |     a    ^^
		# |         / |
		# |  A     /  |
		# |       /   |
		# |b    e/   d|
		# |     /     |
		# |    /      |
		# |   /       |
		# |  /     B  |
		# | /         |
		# V/    c     |
		# #---------->#
		
		corner_A, corner_B = self.corner_lookup[edge_label], self.corner_lookup[~edge_label]
		return [corner_A.edges[1], corner_A.edges[2], corner_B.edges[1], corner_B.edges[2]]
	
	def flip_edge(self, edge_label):
		''' Return a new triangulation obtained by flipping the given edge.
		
		The chosen edge must be flippable. '''
		
		assert(self.is_flippable(edge_label))
		
		# Use the following for reference:
		# #<----------#     #-----------#
		# |     a    ^^     |\          |
		# |         / |     | \         |
		# |  A     /  |     |  \     A2 |
		# |       /   |     |   \       |
		# |b    e/   d| --> |    \e'    |
		# |     /     |     |     \     |
		# |    /      |     |      \    |
		# |   /       |     |       \   |
		# |  /     B  |     | B2     \  |
		# | /         |     |         \ |
		# V/    c     |     |          V|
		# #---------->#     #-----------#
		
		edge_map = dict((edge, Edge(edge.label)) for edge in self.edges)
		
		# Most triangles don't change.
		triangles = [Triangle([edge_map[edge] for edge in triangle]) for triangle in self if edge_label not in triangle.labels and ~edge_label not in triangle.labels]
		
		a, b, c, d = self.square_about_edge(edge_label)
		e = self.edge_lookup[edge_label]
		
		triangle_A2 = Triangle([edge_map[e], edge_map[d], edge_map[a]])
		triangle_B2 = Triangle([edge_map[~e], edge_map[b], edge_map[c]])
		
		return Triangulation(triangles + [triangle_A2, triangle_B2])
	
	def relabel_edges(self, label_map):
		''' Return a new triangulation obtained by relabelling the edges according to label_map. '''
		
		if isinstance(label_map, (list, tuple)):
			label_map = dict(enumerate(label_map))
		else:
			label_map = dict(label_map)
		
		# Build any missing labels.
		for i in self.indices:
			if i in label_map and ~i in label_map:
				pass
			elif i not in label_map and ~i in label_map:
				label_map[i] = ~label_map[~i]
			elif i in label_map and ~i not in label_map:
				label_map[~i] = ~label_map[i]
			else:
				raise curver.AssumptionError('Missing new label for %d.' % i)
		
		edge_map = dict((edge, Edge(label_map[edge.label])) for edge in self.edges)
		triangles = [Triangle([edge_map[edge] for edge in triangle]) for triangle in self]
		
		return Triangulation(triangles)
	
	def find_isometry(self, other, label_map):
		''' Return the isometry from this triangulation to other defined by label_map.
		
		label_map must be a dictionary mapping self.labels to other.labels. Labels may
		be omitted if they are determined by other given ones and these will be found
		automatically.
		
		Assumes (and checks) that such an isometry exists and is unique. '''
		
		assert(isinstance(label_map, dict))
		
		# Make a local copy as we may need to make a lot of changes.
		label_map = dict(label_map)
		
		source_orders = dict([(edge.label, len(vertex_class)) for vertex_class in self.vertex_classes for edge in vertex_class])
		target_orders = dict([(edge.label, len(vertex_class)) for vertex_class in other.vertex_classes for edge in vertex_class])
		# We do a depth first search extending the corner map across the triangulation.
		# This is a stack of labels that may still have consequences to check.
		to_process = [(edge_from_label, label_map[edge_from_label]) for edge_from_label in label_map]
		while to_process:
			from_label, to_label = to_process.pop()
			
			neighbours = [
				(~from_label, ~to_label),
				(self.corner_lookup[from_label].labels[1], other.corner_lookup[to_label].labels[1])
				]
			for new_from_label, new_to_label in neighbours:
				if new_from_label in label_map:
					# Check that this map is still consistent.
					if new_to_label != label_map[new_from_label]:
						raise curver.AssumptionError('This label_map does not extend to an isometry.')
				else:
					# Extend the map.
					if source_orders[new_from_label] != target_orders[new_to_label]:
						raise curver.AssumptionError('This label_map does not extend to an isometry.')
					label_map[new_from_label] = new_to_label
					to_process.append((new_from_label, new_to_label))
		
		if any(i not in label_map for i in self.labels):
			raise curver.AssumptionError('This label_map cannot be extended to an isometry.')
		
		return curver.kernel.Isometry(self, other, label_map)
	
	def isometries_to(self, other):
		''' Return a list of all isometries from this triangulation to other. '''
		
		assert(isinstance(other, Triangulation))
		
		if self.zeta != other.zeta:
			return []
		
		# !?! This needs to be modified to work on disconnected surfaces.
		
		# Isometries are determined by where a single triangle is sent.
		# We take a corner of smallest degree.
		source_cc = min(self.corner_classes, key=len)
		source_corner = source_cc[0]
		# And find all the places where it could be sent so there are as few as possible to check.
		target_corners = [corner for target_cc in other.corner_classes for corner in target_cc if len(target_cc) == len(source_cc)]
		
		isometries = []
		for target_corner in target_corners:
			try:
				isometries.append(self.find_isometry(other, {source_corner.label: target_corner.label}))
			except curver.AssumptionError:
				pass
		
		return isometries
	
	def self_isometries(self):
		''' Return a list of isometries taking this triangulation to itself. '''
		
		return self.isometries_to(self)
	
	def is_isometric_to(self, other):
		''' Return if there are any orientation preserving isometries from this triangulation to other. '''
		
		assert(isinstance(other, Triangulation))
		
		return len(self.isometries_to(other)) > 0
	
	# Curves we can build on this triangulation.
	def lamination(self, geometric, remove_peripheral=True):
		''' Return a new lamination on this surface assigning the specified weight to each edge. '''
		
		if remove_peripheral and False:  # !?!
			# Compute how much peripheral component there is on each corner class.
			# This will also check that the triangle inequalities are satisfied. When
			# they fail one of peripheral.values() is negative, which is non-zero and
			# so triggers the correction.
			def dual_weight(corner):
				''' Return double the weight of normal arc corresponding to the given corner. '''
				
				return geometric[corner.indices[1]] + geometric[corner.indices[2]] - geometric[corner.index]
			
			peripheral = dict((vertex, min(dual_weight(corner) for corner in self.corner_class_of_vertex(vertex))) for vertex in self.vertices)
			if any(peripheral.values()):  # Is there any to add / remove?
				geometric = [geometric[i] - sum(peripheral[v] for v in self.vertices_of_edge(i)) // 2 for i in range(self.zeta)]
		
		
		return curver.kernel.Lamination(self, geometric).promote()
	
	def empty_lamination(self):
		''' Return an empty curve on this surface. '''
		
		return self.lamination([0] * self.zeta, [0] * self.zeta)
	
	def key_curves(self):
		''' Return a list of curves which fill the underlying surface.
		
		As these fill, by Alexander's trick a mapping class is the identity
		if and only if it fixes all of them, including orientation. '''
		
		curves = []
		
		for index in self.indices:
			arc = curver.kernel.Arc(self, [1 if i == index else 0 for i in range(self.zeta)], [1 if i == index else 0 for i in range(self.zeta)])
			for component in arc.boundary().components():
				curves.append(component)
		
		return curves
	
	def id_isometry(self):
		''' Return the isometry representing the identity map. '''
		
		return curver.kernel.Isometry(self, self, dict((i, i) for i in self.labels))
	
	def id_encoding(self):
		''' Return an encoding of the identity map on this triangulation. '''
		
		return self.id_isometry().encode()
	
	def encode_flip(self, edge_label):
		''' Return an encoding of the effect of flipping the given edge.
		
		The given edge must be flippable. '''
		
		assert(self.is_flippable(edge_label))
		
		new_triangulation = self.flip_edge(edge_label)
		
		return curver.kernel.EdgeFlip(self, new_triangulation, edge_label).encode()
	
	def encode_spiral(self, edge_label, power):
		''' Return an encoding of the effect of spiraling about the given edge.
		
		The given edge must be spiralable. '''
		
		return curver.kernel.Spiral(self, self, edge_label, power).encode()
	
	def encode_relabel_edges(self, label_map):
		''' Return an encoding of the effect of flipping the given edge. '''
		
		if isinstance(label_map, (list, tuple)):
			label_map = dict(enumerate(label_map))
		else:
			label_map = dict(label_map)
		
		# Build any missing labels.
		for i in self.indices:
			if i in label_map and ~i in label_map:
				pass
			elif i not in label_map and ~i in label_map:
				label_map[i] = ~label_map[~i]
			elif i in label_map and ~i not in label_map:
				label_map[~i] = ~label_map[i]
			else:
				raise curver.AssumptionError('Missing new label for %d.' % i)
		
		new_triangulation = self.relabel_edges(label_map)
		
		return curver.kernel.Isometry(self, new_triangulation, label_map).encode()
	
	def encode(self, sequence):
		''' Return the encoding given by sequence.
		
		This consists of EdgeFlips, Isometries and LinearTransformations. Furthermore there are
		several conventions that allow these to be specified by a smaller amount of information.
		 - An integer x represents EdgeFlip(..., edge_label=x)
		 - A dictionary which has i or ~i as a key (for every i) represents a relabelling.
		 - A dictionary which is missing i and ~i (for some i) represents an isometry back to this triangulation.
		 - None represents the identity isometry.
		
		This sequence is read in reverse in order respect composition. For example:
			self.encode([1, {1: ~2}, 2, 3, ~4])
		is the mapping class which: flips edge ~4, then 3, then 2, then relabels
		back to the starting triangulation via the isometry which takes 1 to ~2 and
		then finally flips edge 1. '''
		
		assert(isinstance(sequence, (list, tuple)))
		
		h = None
		for item in reversed(sequence):
			if isinstance(item, curver.IntegerType):  # Flip.
				if h is None:
					h = self.encode_flip(item)
				else:
					h = h.target_triangulation.encode_flip(item) * h
			elif isinstance(item, dict):  # Isometry.
				if h is None:
					h = self.encode_relabel_edges(item)
				elif all(i in item or ~i in item for i in self.indices):
					h = h.target_triangulation.encode_relabel_edges(item) * h
				else:  # If some edges are missing then we assume that we must be mapping back to this triangulation.
					h = h.target_triangulation.find_isometry(self, item).encode() * h
			elif item is None:  # Identity isometry.
				if h is None:
					h = self.id_encoding()
				else:
					h = h.target_triangulation.id_encoding() * h
			elif isinstance(item, tuple) and len(item) == 2:  # Spiral
				if h is None:
					h = self.encode_spiral(item[0], item[1])
				else:
					h = h.target_triangulation.encode_spiral(item[0], item[1]) * h
			elif isinstance(item, curver.kernel.Encoding):  # Encoding.
				if h is None:
					h = item
				else:
					h = item * h
			elif isinstance(item, curver.kernel.Move):  # Move.
				if h is None:
					h = item.encode()
				else:
					h = item.encode() * h
			else:  # Other.
				if h is None:
					h = item.encode()
				else:
					h = item.encode() * h
		
		return h

def create_triangulation(cls, edge_labels):
	''' A helper function for pickling. '''
	
	return cls.from_tuple(edge_labels)

