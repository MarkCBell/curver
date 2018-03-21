
''' A module for representing homology classes on triangulations. '''

import curver

class HomologyClass(object):
    ''' This represents a homology class of a triangulation (relative to its vertices). '''
    def __init__(self, triangulation, algebraic):
        assert isinstance(triangulation, curver.kernel.Triangulation)
        
        self.triangulation = triangulation
        self.zeta = self.triangulation.zeta
        self.algebraic = list(algebraic)
        assert len(self.algebraic) == self.zeta
    
    def __repr__(self):
        return str(self)
    def __str__(self):
        return str(self.algebraic)
    def __iter__(self):
        return iter(self.algebraic)
    def __call__(self, edge):
        ''' Return the geometric measure assigned to item. '''
        if isinstance(edge, curver.IntegerType): edge = curver.kernel.Edge(edge)  # If given an integer instead.
        
        return self.algebraic[edge.index] * edge.sign()
    def __eq__(self, other):
        return self.triangulation == other.triangulation and self.canonical().algebraic == other.canonical().algebraic
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        return hash(tuple(self.canonical().algebraic))
    def __add__(self, other):
        if isinstance(other, HomologyClass):
            if other.triangulation != self.triangulation:
                raise ValueError('Homology classes must be on the same triangulation to add them.')
            
            algebraic = [x + y for x, y in zip(self, other)]
            return HomologyClass(self.triangulation, algebraic)
        elif other == 0:  # So we can use sum.
            return self
        else:
            return NotImplemented
    def __radd__(self, other):
        return self + other
    def __neg__(self):
        algebraic = [-x for x in self]
        return HomologyClass(self.triangulation, algebraic)
    def __sub__(self, other):
        return self + (-other)
    def canonical(self):
        ''' Return the canonical form of this HomologyClass.
        
        This is the HomologyClass that is homologous to this one and has weight 0 on each edge of the standard dual tree of the underlying triangulation. '''
        
        return HomologyClass(self.triangulation, curver.kernel.matrix_vector(self.triangulation.homology_matrix(), self.algebraic))
    
    def is_canonical(self):
        ''' Return whether this homology class is already in canonical form.
        
        This is useful for constructing a basis of H_1(S). '''
        
        return self.canonical().algebraic == self.algebraic

