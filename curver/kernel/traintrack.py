
''' A module for representing train tracks on triangulations. '''

import networkx

import curver
from curver.kernel.lamination import Shortenable  # Special import needed for subclassing.
from curver.kernel.utilities import memoize  # Special import needed for decorating.

class TrainTrack(Shortenable):
    ''' A Lamination in which each triangle is tripod free. '''
    
    def shorten_strategy(self, edge):
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        if not self.triangulation.is_flippable(edge):
            return 0
        
        a, b, c, d, e = self.triangulation.square(edge)
        ad, bd, cd, dd, ed = [self.dual_weight(edgy) for edgy in self.triangulation.square(edge)]
        if self(e) <= 0:
            return 0
        if b == ~d and ad > 0 and bd == 0 and ed == 0:
            return 0
        if a == ~c and ad == 0 and bd > 0 and ed == 0:
            return 0
        
        num_bad = [ed > 0, ad < 0, bd < 0].count(True)
        num_good = [ed < 0, ad > 0, bd > 0].count(True)
        
        if num_bad > 0:
            return 0.25
        
        if num_good == 3:
            return 1
        elif num_good == 2:
            return 1
        else:  # num_good == 1:
            if (ad > 0 and cd > 0) or (bd > 0 and dd > 0):
                return 0.5
            else:
                return 0.75
    
    @memoize()
    def components(self):
        short, conjugator = self.shorten()
        
        components = dict()
        for edge in short.triangulation.positive_edges:  # Only need to check half of them.
            # Check for an Arc here.
            if short(edge) < 0:
                component, multiplicity = short.triangulation.edge_arc(edge), abs(short(edge))
                components[conjugator.inverse()(component)] = multiplicity  # Map it back onto self.
            elif short(edge) > 0:
                if short.triangulation.is_flippable(edge):  # Check for a non-peripheral curve here.
                    a, b, c, d, e = short.triangulation.square(edge)
                    da, db, dc, dd, de = [short.dual_weight(edgy) for edgy in short.triangulation.square(edge)]
                    if b == ~d and db == 0 and de == 0:
                        geometric = [1 if index == b.index or index == e.index else 0 for index in short.triangulation.indices]
                        component, multiplicity = curver.kernel.Curve(short.triangulation, geometric), short(e)
                        components[conjugator.inverse()(component)] = multiplicity  # Map it back onto self.
                else:  # Check for a peripheral curve here.
                    geometric = [1 if index == edge.index else 0 for index in short.triangulation.indices]
                    component, multiplicity = curver.kernel.Curve(short.triangulation, geometric), short(edge)
                    components[conjugator.inverse()(component)] = multiplicity  # Map it back onto self.
        
        return components
    
    def vertex_cycles(self):
        ''' Yield the vertex cycles of this train track.
        
        These are the curves carried by this train track that run over each brach at most twice.
        Be careful as there are often a *lot* of them. '''
        
        def connected_to(edge):
            ''' Yield the edges you can reach by travelling out of the given edge. '''
            corner = self.triangulation.corner_lookup[edge.label]
            if self.dual_weight(corner[1]): yield ~corner[2]
            if self.dual_weight(corner[2]): yield ~corner[1]
        
        # Build graph.
        edges = [(edge, edgy) for edge in self.triangulation.edges for edgy in connected_to(edge)]
        G = networkx.DiGraph(edges)
        
        for cycle in networkx.simple_cycles(G):
            curve = self.triangulation.lamination_from_cut_sequence(cycle)
            if isinstance(curve, curver.kernel.Curve):
                yield curve
    
    def vertex_cycle(self):
        ''' Return a vertex cycle of this train track. '''
        
        return next(iter(self.vertex_cycles()))

