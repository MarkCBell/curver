
''' A module for representing train tracks on triangulations. '''

import networkx

import curver
from curver.kernel.lamination import Lamination  # Special import needed for subclassing.
from curver.kernel.utilities import memoize  # Special import needed for decorating.

class TrainTrack(Lamination):
    ''' A Lamination in which each triangle is tripod free. '''
    
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

