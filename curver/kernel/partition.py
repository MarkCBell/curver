
''' A module for representing partition graphs. '''

from itertools import permutations, combinations_with_replacement
import networkx

class CurvePartitionGraph(object):
    ''' This represents a partition graph of a MultiCurve. '''
    def __init__(self, multicurve, graph):
        self.multicurve = multicurve
        self.graph = graph
    
    def __str__(self):
        G = self.graph
        N = G.nodes(data=True)
        
        return min(','.join('%d:%d' % (N[p]['genus'], N[p]['vertices']) for p in perm) + '-' + ','.join([':'.join(str(x) for x in sorted(d['weight'] for d in G.get_edge_data(i, j, default={}).values())) for i, j in combinations_with_replacement(perm, r=2)]) for perm in permutations(range(len(G))))
    
    def __repr__(self):
        return str(self)
    
    def __eq__(self, other):
        node_match = lambda v1, v2: v1['genus'] == v2['genus'] and v1['vertices'] == v2['vertices']
        edge_match = lambda e1, e2: sorted(e['weight'] for e in e1.values()) == sorted(e['weight'] for e in e2.values())
        return networkx.is_isomorphic(self.graph, other.graph, node_match=node_match, edge_match=edge_match)
    def __ne__(self, other):
        return not self == other
    
    def __hash__(self):
        return hash(tuple(sorted((d['genus'], d['vertices']) for _, d in self.graph.nodes(data=True))) + tuple(sorted([e['weight'] for _, _, e in self.graph.edges(data=True)])))

