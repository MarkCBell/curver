
import networkx
import curver

class CurvePartitionGraph(object):
    def __init__(self, multicurve, graph):
        self.multicurve = multicurve
        self.graph = graph
    
    def __str__(self):
        return '-'.join('%d,%d' % (g, v) for g, v in sorted((d['genus'], d['vertices']) for _, d in self.graph.nodes(data=True))) + \
                ' (' + '-'.join([str(e['weight']) for _, _, e in self.graph.edges(data=True)]) + ')'
    def __repr__(self):
        return str(self)
    
    def __eq__(self, other):
        node_match = lambda v1, v2: v1['genus'] == v2['genus'] and v1['vertices'] == v2['vertices']
        edge_match = lambda e1, e2: sorted(e['weight'] for e in e1.values()) == sorted(e['weight'] for e in e2.values())
        return networkx.is_isomorphic(self.graph, other.graph, node_match=node_match, edge_match=edge_match)
    def __ne__(self, other):
        return not (self == other)
    
    def __hash__(self):
        return hash(tuple(sorted((d['genus'], d['vertices']) for _, d in self.graph.nodes(data=True))) + tuple(sorted([e['weight'] for _, _, e in self.graph.edges(data=True)])))

