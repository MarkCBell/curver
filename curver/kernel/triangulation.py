
''' A module for representing a triangulation of a punctured surface.

Provides: Edge, Triangle and Triangulation. '''

from math import factorial
from itertools import groupby

import curver
from curver.kernel.utilities import memoize  # Special import needed for decorating.

INFTY = float('inf')

def norm(number):
	''' A map taking an edges label to its index.
	
	That is, x and ~x should map to the same thing. '''
	
	return max(number, ~number)

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
		return (self.__class__, (self.label,))
	def __eq__(self, other):
		if isinstance(other, Edge):
			return self.label == other.label
		elif isinstance(other, curver.IntegerType):
			return self.label == other
		else:
			return NotImplemented
	def __ne__(self, other):
		return not (self == other)
	def __hash__(self):
		return hash(self.label)
	
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
	__slots__ = ['edges', 'labels', 'indices']
	
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
	def __eq__(self, other):
		if isinstance(other, Triangle):
			return self.labels == other.labels
		else:
			return NotImplemented
	def __ne__(self, other):
		return not (self == other)
	def __hash__(self):
		return hash(tuple(self.labels))
	def __len__(self):
		return 3  # This is needed for revered(triangle) to work.
	
	# Note that this is NOT the same convention as used in pieces.
	# There iterating and index accesses return vertices.
	def __iter__(self):
		return iter(self.edges)
	
	def __getitem__(self, index):
		return self.edges[index]
	def __contains__(self, other):
		return other in self.edges

# Remark: In other places in the code you will often see L(triangulation). This is the space
# of laminations on triangulation with the coordinate system induced by the triangulation.

class Triangulation(object):
	''' This represents a triangulation of a punctured surface.
	
	It is specified by a list of Triangles. Its edges must be numbered 0, 1, ... '''
	def __init__(self, triangles):
		# We will sort the triangles into a canonical ordering, the one where the edges are ordered
		# minimally by label. This allows for fast comparisons.
		self.triangles = sorted(triangles, key=lambda t: t.labels)
		
		self.edges = [edge for triangle in self for edge in triangle]
		self.positive_edges = [edge for edge in self.edges if edge.sign() == +1]
		self.labels = sorted([edge.label for edge in self.edges])
		self.indices = sorted([edge.index for edge in self.positive_edges])
		
		self.num_triangles = len(self.triangles)
		self.zeta = len(self.positive_edges)
		assert(self.zeta == self.num_triangles * 3 // 2)
		# Check that the edges have indices 0, ..., zeta-1.
		assert(set(self.labels) == set([i for i in range(self.zeta)] + [~i for i in range(self.zeta)]))
		
		self.edge_lookup = dict((edge.label, edge) for edge in self.edges)
		self.triangle_lookup = dict((edge.label, triangle) for triangle in self for edge in triangle)
		self.corner_lookup = dict((edge.label, Triangle(triangle.edges, rotate=index)) for triangle in triangles for index, edge in enumerate(triangle))
		
		# Group the edges into vertices and ordered anti-clockwise.
		# Here two edges are in the same class iff they have the same tail.
		unused = set(self.edges)
		self.vertices = []
		while unused:
			vertex = [unused.pop()]
			while True:
				neighbour = ~self.corner_lookup[vertex[-1].label][2]
				if neighbour in unused:
					vertex.append(neighbour)
					unused.remove(neighbour)
				else:
					break
			
			self.vertices.append(tuple(vertex))
		
		self.vertex_lookup = dict((edge.label, vertex) for vertex in self.vertices for edge in vertex)
		
		self.num_vertices = len(self.vertices)
		
		self.euler_characteristic = -self.zeta // 3  # = V - E + F since 3F = 2E and V = 0.
		
		# Two triangualtions are the same if and only if they have the same signature.
		self.signature = [e.label for t in self for e in t]
	
	@classmethod
	def from_tuple(cls, *edge_labels):
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
		
		if len(edge_labels) == 1: edge_labels = edge_labels[0]
		
		assert(isinstance(edge_labels, (list, tuple)))
		assert(all(isinstance(labels, (list, tuple)) for labels in edge_labels))
		assert(all(len(labels) == 3 for labels in edge_labels))
		assert(len(edge_labels) > 0)
		
		zeta = len(edge_labels) * 3 // 2
		
		# Check that each of 0, ..., zeta-1, ~0, ..., ~(zeta-1) occurs exactly once.
		flattened = set(label for labels in edge_labels for label in labels)
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
	def __hash__(self):
		return hash(tuple(self.signature))
	
	def max_order(self):
		''' Return the maximum order of a mapping class on this surface. '''
		
		# List of pairs of #vertices and #edges for each component.
		VE = [(len([vertex for vertex in self.vertices if list(vertex)[0] in component]), len(component) // 2) for component in self.components()]
		# List of pairs of genus and #vertices edges for each component.
		GV = [((2 - v + e // 3) // 2, v) for v, e in VE]
		
		def order(g, v):
			''' Return the maximum order of a periodic mapping class on S_{g, v}. '''
			# These bounds follow from the 4g + 4 bound on the closed surface [FarbMarg12]
			# and the Riemann removable singularity theorem which allows us to cap off the
			# punctures when the g > 1 without affecting this bound.
			if g > 1:
				return 4*g + 2
			elif g == 1:
				return max(v, 6)
			else:  # g == 0:
				return v
		
		# List of pairs of orders and multiplicity.
		OM = [(order(g, v), len(list(group))) for (g, v), group in groupby(sorted(GV))]
		
		product = 1
		for o, m in OM:
			product *= o * factorial(m)  # Actually need o*lcm(1, 2, ..., m), but this is easier (and not significantly larger?)
		
		return product
	
	def components(self):
		''' Return a list of sets of the edges in each component of self. '''
		
		classes = curver.kernel.UnionFind(self.edges)
		for edge in self.edges:
			classes.union(edge, ~edge)
		for triangle in self:
			classes.union(triangle)
		
		return list(classes)
	
	def dual_tree(self, avoid=None):
		''' Return a maximal tree in 1--skeleton of the dual of this triangulation.
		
		This are given as lists of Booleans signaling if each edge is in the tree.
		Note that when this surface is disconnected this tree is actually a forest.
		To make this unique / well-defined we return the numerically first one.
		
		If avoid is provided then none of these indices will be set in the dual tree. '''
		
		if avoid is None: avoid = set()
		
		# Kruskal's algorithm.
		dual_tree = [False] * self.zeta
		classes = curver.kernel.UnionFind(self.triangles)
		for index in self.indices:
			if index not in avoid:
				a, b = self.triangle_lookup[index], self.triangle_lookup[~index]
				if classes(a) != classes(b):
					classes.union(a, b)
					dual_tree[index] = True
		
		return dual_tree
	
	@memoize
	def homology_matrix(self):
		''' Return a matrix that kills the entries of the dual tree. '''
		
		dual_tree = self.dual_tree()
		
		M = []
		for index in self.indices:
			row = [0] * self.zeta
			if dual_tree[index]:
				edge = self.edge_lookup[index]
				while True:
					corner = self.corner_lookup[edge.label]
					edge = corner.edges[2]
					if not dual_tree[edge.index]:
						row[edge.index] -= edge.sign()
					else:
						edge = ~edge
					if edge.label == ~index: break
			else:
				row[index] = 1
			M.append(row)
		
		return list(zip(*M))  # Transpose the matrix. We need list so that this is not a generator in Python3.
	
	def is_flippable(self, edge):
		''' Return if the given edge is flippable.
		
		An edge is flippable if and only if it lies in two distinct triangles. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.edge_lookup[edge]  # If given an integer instead.
		
		return self.triangle_lookup[edge.label] != self.triangle_lookup[~edge.label]
	
	def folded_boundary(self, edge):
		''' Return the edge bounding the once-punctured monogon containing the given edge.
		
		The given edge must not be flippable. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.edge_lookup[edge]  # If given an integer instead.
		
		assert(not self.is_flippable(edge))
		
		# As edge_label is not flippable the triangle containing it must be (edge_label, ~edge_label, boundary_edge).
		[boundary_edge] = [edgy for edgy in self.triangle_lookup[edge.label] if edgy.index != edge.index]
		return boundary_edge
	
	def square(self, edge):
		''' Return the four edges around the given edge and the diagonal.
		
		The given edge must be flippable. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.edge_lookup[edge]  # If given an integer instead.
		
		assert(self.is_flippable(edge))
		
		# Given the e, return the edges a, b, c, d, e in order.
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
		
		corner_A, corner_B = self.corner_lookup[edge.label], self.corner_lookup[~edge.label]
		return [corner_A.edges[1], corner_A.edges[2], corner_B.edges[1], corner_B.edges[2], edge]
	
	# Build new triangulations:
	def flip_edge(self, edge):
		''' Return a new triangulation obtained by flipping the given edge.
		
		The chosen edge must be flippable. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.edge_lookup[edge]  # If given an integer instead.
		
		assert(self.is_flippable(edge))
		
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
		triangles = [Triangle([edge_map[edgy] for edgy in triangle]) for triangle in self if edge not in triangle and ~edge not in triangle]
		
		a, b, c, d, e = self.square(edge)
		
		if edge.sign() == +1:
			triangle_A2 = Triangle([edge_map[e], edge_map[d], edge_map[a]])
			triangle_B2 = Triangle([edge_map[~e], edge_map[b], edge_map[c]])
		else:  # edge.sign() == -1.
			triangle_A2 = Triangle([edge_map[~e], edge_map[d], edge_map[a]])
			triangle_B2 = Triangle([edge_map[e], edge_map[b], edge_map[c]])
		return Triangulation(triangles + [triangle_A2, triangle_B2])
	
	def all_encodings(self, num_flips):
		''' Yield all encodings that can be made using at most the given number of flips.
		
		Runs in exp(num_flips) time. '''
		
		if num_flips == 0:
			yield self.id_encoding()
		
		# TODO: 2) Make efficient by using the fact that disjoint flips commute.
		
		for edge in self.positive_edges:
			step = self.encode_flip(edge)
			for encoding in step.target_triangulation.all_encodings(num_flips-1):
				yield encoding * step
		
		return
	
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
		automatically. Additionally, if an entire component is omitted then we assume
		that the map is the identity on it.
		
		Assumes (and checks) that such an isometry exists and is unique. '''
		
		assert(isinstance(label_map, dict))
		
		# Make a local copy as we may need to make a lot of changes.
		label_map = dict(label_map)
		
		source_orders = dict([(edge.label, len(vertex)) for vertex in self.vertices for edge in vertex])
		target_orders = dict([(edge.label, len(vertex)) for vertex in other.vertices for edge in vertex])
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
		
		# If an entire component is omitted then assume the map is the identity on it.
		used_labels = set(label_map.keys() + label_map.values())
		for component in self.components():
			if not any(edge.label in used_labels for edge in component):
				for edge in component:
					label_map[edge.label] = edge.label
		
		if any(i not in label_map for i in self.labels):
			raise curver.AssumptionError('This label_map cannot be extended to an isometry.')
		
		return curver.kernel.Isometry(self, other, label_map)
	
	def isometries_to(self, other):
		''' Return a list of all isometries from this triangulation to other. '''
		
		assert(isinstance(other, Triangulation))
		
		if self.zeta != other.zeta:
			return []
		
		# TODO: 3) Modify this to work on disconnected surfaces.
		
		# Isometries are determined by where a single triangle is sent.
		# We take a corner of smallest degree.
		source_vertex = min(self.vertices, key=len)
		source_edge = source_vertex[0]
		# And find all the places where it could be sent so there are as few as possible to check.
		target_edges = [edge for target_vertex in other.vertices for edge in target_vertex if len(target_vertex) == len(source_vertex)]
		
		isometries = []
		for target_edge in target_edges:
			try:
				isometries.append(self.find_isometry(other, {source_edge.label: target_edge.label}))
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
	
	# Laminations we can build on this triangulation.
	def lamination(self, weights):
		''' Return a new lamination on this surface assigning the specified weight to each edge. '''
		
		assert(len(weights) == self.zeta)
		# Should check all dual weights.
		
		return curver.kernel.Lamination(self, weights).remove_peripheral().promote()
	
	def empty_lamination(self):
		''' Return an empty curve on this surface. '''
		
		return self.lamination([0] * self.zeta)
	
	def as_lamination(self):
		''' Return this triangulation as a lamination. '''
		
		return self.lamination([-1] * self.zeta)
	
	def edge_arcs(self):
		''' Return a list containing the Arc representing each Edge.
		
		As these fill, by Alexander's trick a mapping class is the identity if and only if it fixes all of them. '''
		
		return [curver.kernel.Arc(self, [0 if i != index else -1 for i in range(self.zeta)]) for index in self.indices]  # Could use self.lamination.
	
	def edge_homologies(self):
		''' Return a list containing the HomologyClass of each Edge. '''
		
		# Could skip those in self.dual_tree().
		return [curver.kernel.HomologyClass(self, [0 if i != index else 1 for i in range(self.zeta)]) for index in self.indices]
	
	def id_isometry(self):
		''' Return the isometry representing the identity map. '''
		
		return curver.kernel.Isometry(self, self, dict((i, i) for i in self.labels))
	
	def id_encoding(self):
		''' Return an encoding of the identity map on this triangulation. '''
		
		return self.id_isometry().encode()
	
	def encode_flip(self, edge):
		''' Return an encoding of the effect of flipping the given edge.
		
		The given edge must be flippable. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.edge_lookup[edge]  # If given an integer instead.
		
		assert(self.is_flippable(edge))
		
		new_triangulation = self.flip_edge(edge)
		
		return curver.kernel.EdgeFlip(self, new_triangulation, edge).encode()
	
	def encode_spiral(self, edge, power):
		''' Return an encoding of the effect of spiraling about the given edge.
		
		The given edge must be spiralable. '''
		
		if isinstance(edge, curver.IntegerType): edge = self.edge_lookup[edge]  # If given an integer instead.
		
		return curver.kernel.Spiral(self, self, edge, power).encode()
	
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
		 - A pair of integers (x, k) represents a Spiral about x to the power k.
		
		This sequence is read in reverse in order respect composition. For example:
			self.encode([1, {1: ~2}, 2, 3, ~4])
		is the mapping class which: flips edge ~4, then 3, then 2, then relabels
		back to the starting triangulation via the isometry which takes 1 to ~2 and
		then finally flips edge 1. '''
		
		assert(isinstance(sequence, (list, tuple)))
		assert(len(sequence) > 0)
		
		T = self
		h = None
		for item in reversed(sequence):
			if isinstance(item, curver.IntegerType):  # Flip.
				g = T.encode_flip(item)
			elif isinstance(item, dict):  # Isometry.
				if all(i in item or ~i in item for i in self.indices):
					g = T.encode_relabel_edges(item)
				else:  # If some edges are missing then we assume that we must be mapping back to this triangulation.
					g = T.find_isometry(self, item).encode()
			elif item is None:  # Identity isometry.
				g = T.id_encoding()
			elif isinstance(item, tuple) and len(item) == 2 and all(isinstance(x, curver.IntegerType) for x in item):  # Spiral. TODO 4) Update to Twist and HalfTwist!
				g = T.encode_spiral(item[0], item[1])
			elif isinstance(item, curver.kernel.Encoding):  # Encoding.
				g = item
			elif isinstance(item, curver.kernel.Move):  # Move.
				g = item.encode()
			else:  # Other.
				g = item.encode()
			
			if h is None:
				h = g
			else:
				h = g * h
			T = h.target_triangulation
		
		return h

def create_triangulation(cls, edge_labels):
	''' A helper function for pickling. '''
	
	return cls.from_tuple(edge_labels)

