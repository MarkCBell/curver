
''' A module for representing and manipulating finite subgroups of a mapping class group. '''

from collections import defaultdict, namedtuple
from copy import deepcopy
from fractions import Fraction
from itertools import groupby
try:
    from Queue import Queue
except ImportError:  # Python3.
    from queue import Queue

import curver
from curver.kernel.decorators import memoize, ensure

ConePoint = namedtuple('ConePoint', ['punctured', 'order', 'holonomy', 'preimages'])
Orbifold = namedtuple('Orbifold', ['euler_characteristic', 'preimages', 'cone_points'])

class FiniteSubgroup(object):
    ''' This represents a finite subgroup of a mapping class group. '''
    def __init__(self, mapping_classes, generators=None):
        self.mapping_classes = mapping_classes  # Dict: name |--> mapping class.
        self.generators = sorted(self.mapping_classes) if generators is None else generators
        self.triangulation = list(self.mapping_classes.values())[0].source_triangulation
        # asserts?
    
    def __str__(self):
        return '< ' + ', '.join(str(name) for name in self.generators) + ' >'
    def __repr__(self):
        return str(self)
    def __iter__(self):
        return iter(self.mapping_classes)
    def __len__(self):
        return len(self.mapping_classes)
    def __getitem__(self, item):
        return self.mapping_classes[item]
    def __eq__(self, other):
        if isinstance(other, FiniteSubgroup):
            return self.mapping_classes == other.mapping_classes and self.generators == other.generators
        else:
            return NotImplemented
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(tuple(self[letter] for letter in sorted(self)))
    
    @classmethod
    def from_generators(cls, generators):
        ''' Build the FiniteSubgroup from these generators (a dict mapping names to mapping classes).
        
        Currently this does not check that the subgroup generated is finite. '''
        
        mapping_classes = dict(generators)
        seen = set(generators.values())
        to_check = Queue()
        for word in sorted(generators):
            to_check.put(word)
        
        while not to_check.empty():
            word = to_check.get()
            current = mapping_classes[word]
            for letter in sorted(generators):
                generator = generators[letter]
                if letter + word not in mapping_classes:
                    neighbour = generator * current
                    if neighbour not in seen:
                        mapping_classes[letter + word] = neighbour
                        to_check.put(letter + word)
                        seen.add(neighbour)
                        # Once we know the correct bound we can add:
                        # if len(seen) > 84 * (g-1):
                        #   raise ValueError('Mapping classes do not generate a finite subgroup.')
        
        return cls(mapping_classes, list(generators))
    
    @memoize
    @ensure(lambda data: data.result.is_polygonalisation, lambda data: all(data.self[word](data.result) == data.result for word in data.self.generators))
    def invariant_polygonalisation(self):
        ''' Return a multiarc that is a polygonalisation and is invariant under self. '''
        
        H = [self[name] for name in self]
        # ||H|| := max(||h|| for h in H) which in this case is at most zeta ||self||.
        
        def orbit(a):
            ''' Yield the orbit of a under H. '''
            
            images = set()
            for h in H:
                image = h(a)
                if image not in images:
                    yield image
                    images.add(image)
        
        conjugator = self.triangulation.id_encoding()
        triangulation = conjugator.target_triangulation
        invariant_multiarc = triangulation.empty_lamination()
        while not invariant_multiarc.is_polygonalisation():  # Loops at most zeta times.
            # Start by constructing a list of lists of images of edges under elements in H. So
            #   H_images[i][j] = H[i](edge_j)
            # These make it easy to test whether the H--orbit of edge_j is embedded: this occurs iff
            #   H_images[i][j](j) <= 0 for every 0 <= i < |H|.
            # Note that since we will be needing this repeatedly, we will create the original here and
            # make a deepcopy at the start of each calculation later.
            original_H_images = [[h(arc) for arc in triangulation.edge_arcs()] for h in H]
            
            # Initially we have to check every edge, so we do this once here to avoid repeating it for every image in the next loop.
            for edge in triangulation.positive_edges:
                if invariant_multiarc(edge) == 0 and all(images[edge.index](edge) <= 0 for images in original_H_images):  # if not existing component and H--orbit is embedded.
                    arc = triangulation.edge_arc(edge)
                    invariant_multiarc = triangulation.disjoint_sum([invariant_multiarc] + list(orbit(arc)))  # Add it to the invariant arc.
                    break
            else:  # If that doesn't work then we can search the unicorn arcs.
                # Theorem: For any arc a there is an h in H such that there is a unicorn of a and h(a) whose H--orbit is embedded.
                # In fact if a is an arc that does not cut off a disk in S - invariant_multiarc then the obtained unicorn is also disjoint and not a component of invariant_multiarc.
                # Such an arc exists since invariant_multiarc is not a polygonalisation, and in fact one of the edges of triangulation must be one since invariant_multiarc is short.
                dual_tree = triangulation.dual_tree(avoid={edge for edge in triangulation.positive_edges if invariant_multiarc(edge) < 0})
                arc = triangulation.edge_arc([edge for edge in triangulation.positive_edges if edge.index not in dual_tree and invariant_multiarc(edge) == 0][0])
                
                done = False
                for image in orbit(arc):  # Loops at most |H| times.
                    # Check the unicorn arcs that can be made from arc and image.
                    H_images = deepcopy(original_H_images)  # Get a fresh copy to work with.
                    image_conjugator = image.shorten(drop=0)  # The Mosher sequence from image back to arc, this contains all the unicorn arcs.
                    # Theorem: Since arc is short, the set of arcs that appear in the Mosher flip sequence includes all unicorns made from arc and image.
                    for index, move in enumerate(reversed(image_conjugator)):  # Loops at most ||H|| times.
                        # Note that each of the following modifications runs in O(||H||).
                        # Currently, by induction, H_images[i][j] = (prefix * H[i] * ~prefix)(edge_j) where prefix = image_conjugator[-index:].
                        # We tackle this in four steps, starting with the easy side.
                        # 1) Update so that H_images[i][j] = (move * prefix * H[i] * ~prefix)(edge_j).
                        H_images = [[move(arcy) for arcy in images] for images in H_images]
                        # Now for the hard side. To make this easy, we use the fact that H_images[i].transpose() records the inverse map.
                        # 2) Update so that H_images[i][j] = (prefix * ~H[i] * ~prefix * ~move)(edge_j).
                        H_images = [[curver.kernel.Arc(move.source_triangulation, geometric) for geometric in zip(*images)] for images in H_images]
                        # 3) Update so that H_images[i][j] = (move * prefix * ~H[i] * ~prefix * ~move)(edge_j).
                        H_images = [[move(arcy) for arcy in images] for images in H_images]
                        # Finally flip back by using the same trick.
                        # 4) Update so that H_images[i][j] = (move * prefix * H[i] * ~prefix * ~move)(edge_j)
                        H_images = [[curver.kernel.Arc(move.target_triangulation, geometric) for geometric in zip(*images)] for images in H_images]
                        # Now H_images[i][j] = (next_prefix * H[i] * ~next_prefix)(edge_j) where next_prefix = image_conjugator[-index-1:].
                        
                        if not isinstance(move, curver.kernel.EdgeFlip):
                            continue
                        
                        edge = move.edge  # Only have to check the one new edge that has appeared.
                        # This arc is not already in invariant_multiarc so we only have to check ...
                        if any(images[edge.index](edge) > 0 for images in H_images):  # if H--orbit is not embedded.
                            continue
                        
                        arcy = move.target_triangulation.edge_arc(edge)
                        prefix = image_conjugator[-1-index:]
                        unicorn = prefix.inverse()(arcy)  # Pull it back.
                        invariant_multiarc = triangulation.disjoint_sum([invariant_multiarc] + list(orbit(unicorn)))  # Add it to the invariant arc.
                        # Break out and start again.
                        done = True
                        break
                    if done: break
                else:
                    raise RuntimeError('Unable to find invariant unicorn arc.')
            
            # Reshorten invariant_multiarc.
            next_conjugator = invariant_multiarc.shorten()
            conjugator = next_conjugator * conjugator
            invariant_multiarc = next_conjugator(invariant_multiarc)
            H = [next_conjugator * h * next_conjugator.inverse() for h in H]
            triangulation = conjugator.target_triangulation
        
        return conjugator.inverse()(invariant_multiarc)
    
    @memoize
    def quotient_orbifold_signature(self):
        ''' Return the signature of self.surface() / self.
        
        This records the covering map via the `order`, `holonomy` and `preimage` fields. '''
        
        polygonalisation = self.invariant_polygonalisation()
        
        conjugator = polygonalisation.shorten()
        short = conjugator(polygonalisation)
        
        H = dict((name, conjugator * self[name] * conjugator.inverse()) for name in self)
        
        # Some short names.
        triangulation = short.triangulation
        components = short.triangulation.components()
        surface = triangulation.surface()
        order = len(H)
        
        OrientedArc = namedtuple('OrientedArc', ['arc', 'hc', 'boundary'])
        # Build the oriented arcs.
        oriented = dict()
        for edge in triangulation.edges:
            if short(edge) < 0:
                arc = triangulation.edge_arc(edge)
                hc = triangulation.edge_homology(edge)
                if len(arc.vertices()) == 1:
                    [v] = arc.vertices()
                    v_edges = curver.kernel.utilities.cyclic_slice(v, edge, ~edge)
                    boundary = triangulation.curve_from_cut_sequence(v_edges[1:])
                else:  # two vertices:
                    boundary = arc.boundary()
                oriented[edge] = OrientedArc(arc, hc, boundary)
        oriented_arcs = list(oriented.values())
        
        # Build some useful maps.
        # A) For each pair of oriented arcs, the set of H (names) that map one to the other.
        pairs = dict((a, dict((b, set()) for b in oriented_arcs)) for a in oriented_arcs)
        for a in oriented_arcs:
            for name, h in H.items():
                b = OrientedArc(h(a.arc), h(a.hc), h(a.boundary))
                pairs[a][b].add(name)
        
        # B) The component of triangulation each oriented_arc lives in.
        component_lookup = dict((oriented[edge], component) for component in components for edge in component if edge in oriented)
        
        # C) The size of the orbit of each component of triangulation under the action of h.
        classes = curver.kernel.UnionFind(components)
        for a in oriented_arcs:
            for b in oriented_arcs:
                if pairs[a][b]:
                    classes.union(component_lookup[a], component_lookup[b])
        component_orbit_size = dict()
        for cls in classes:
            for component in cls:
                component_orbit_size[component] = len(components)
        
        # D) The (orbifold) Euler characteristic of each components quotient.
        euler_characteristic = dict((component, Fraction(surface[component].chi * component_orbit_size[component], order)) for component in components)
        
        # Compute the polygons cut out by short.
        # Remember to walk around their boundary in the correct direction.
        polygons = []
        used = set()
        for edge in triangulation.edges:
            if short(edge) < 0 and edge not in used:
                polygon_edges = [edge]
                while True:
                    # Correct direction requires taking the last oriented arc out of the vertex each time.
                    polygon_edges.append([edgy for edgy in curver.kernel.utilities.cyclic_slice(triangulation.vertex_lookup[~polygon_edges[-1]], ~polygon_edges[-1]) if short(edgy) < 0][-1])
                    if polygon_edges[-1] == polygon_edges[0]:  # if back where we started.
                        polygon_edges = polygon_edges[:-1]  # Remeber to discard the last edge as it duplicates the first.
                        break
                
                used = used.union(polygon_edges)  # Mark everything as used.
                polygons.append(polygon_edges)
        
        # There are three places to look for cone points:
        # 1) at vertices,
        # 2) at the midpoints of edges, and
        # 3) at the centres of polygons (note that these are not once-punctured polygons since we started with an invariant polygonalisation).
        candidates = [(True, [edge for edge in vertex if short(edge) < 0]) for vertex in triangulation.vertices] + \
            [(False, [edge, ~edge]) for edge in triangulation.positive_edges if short(edge) < 0] + \
            [(False, [edge for edge in polygon]) for polygon in polygons]
        
        cone_points = defaultdict(list)
        for punctured, edges in candidates:
            oriented_arcs = [oriented[edge] for edge in edges]
            start = oriented_arcs[0]
            places = [a for a in oriented_arcs[1:] + oriented_arcs[:1] if pairs[start][a]]  # ??
            
            c_order = len(places)
            holonomy = sorted(pairs[start][places[0]])
            preimages = order // c_order
            
            if punctured or len(places) > 1:  # Real cone point.
                cone_points[component_lookup[start]].append(ConePoint(punctured, c_order, holonomy, preimages))  # Put it in the correct component.
        
        # Remember to make all the data canonical by sorting.
        signature = sorted(Orbifold(euler_characteristic[component], component_orbit_size[component], sorted(cone_points[component])) for component in components)
        
        # Group.
        signature = [key for key, group in groupby(signature) for _ in range(len(list(group)) // key.preimages)]
        signature = [Orbifold(orbifold.euler_characteristic, orbifold.preimages, [key for key, group in groupby(orbifold.cone_points) for _ in range(len(list(group)) // key.preimages)]) for orbifold in signature]
        return signature

    def is_conjugate_to(self, other):
        ''' Return whether this subgroup is conjugate to other.
        
        That is, whether there is a mapping class that conjugates self[name] to other[name] for each name in self. '''
        
        assert isinstance(other, curver.kernel.FiniteSubgroup)
        
        if self.triangulation.surface() != other.triangulation.surface():  # Defined on different surfaces.
            return False
        
        if len(self) != len(other):  # Conjugacy invariant.
            return False
        
        return self.quotient_orbifold_signature() == other.quotient_orbifold_signature()  # Total conjugacy invariant.

