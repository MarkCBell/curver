
''' A module of data structures. '''

from collections import defaultdict, namedtuple
from itertools import chain, islice

import curver

class UnionFind(object):
    ''' A fast union--find data structure. Given items must be hashable. '''
    def __init__(self, items):
        self.parent = dict((item, item) for item in items)
        self.rank = dict((item, 0) for item in items)
    def __iter__(self):
        ''' Iterate through the groups of self. '''
        groups = defaultdict(list)
        for item in self.parent:
            groups[self(item)].append(item)
        return iter(groups.values())
    def __len__(self):
        return sum(1 if self.parent[item] == item else 0 for item in self.parent)
    def __repr__(self):
        return str(self)
    def __str__(self):
        return ', '.join('{' + ', '.join(str(item) for item in g) + '}' for g in self)
    def __call__(self, x):
        ''' Find the root of x. Two items are in the same group iff they have the same root. '''
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            x, self.parent[x] = self.parent[x], root
        return root
    def union2(self, x, y):
        ''' Combine the class containing x and the class containing y. '''
        rx, ry = self(x), self(y)
        if self.rank[x] > self.rank[y]:
            self.parent[ry] = rx
        elif self.rank[x] < self.rank[y]:
            self.parent[rx] = ry
        elif rx != ry:
            self.parent[ry] = rx
            self.rank[rx] += 1
    def union(self, *args):
        ''' Combine all of the classes containing the given items. '''
        if len(args) == 1: args = args[0]
        for item in args:
            self.union2(args[0], item)


Terminal = namedtuple('Terminal', ['value'])

class StraightLineProgram(object):
    ''' This represents a straight line program. '''
    def __init__(self, data):
        if isinstance(data, StraightLineProgram):  # Copy.
            data = data.graph
        elif isinstance(data, (list, tuple)):  # Wrap.
            if not data or any(not isinstance(children, (list, tuple)) for children in data) or any(not isinstance(child, (Terminal, curver.IntegerType)) for children in data for child in children):
                data = [[Terminal(child) for child in data]]
        else:
            raise ValueError('Unknown data.')
        
        self.graph = [tuple(children) for children in data]
        self.sinks = [item.value for lst in self.graph for item in lst if isinstance(item, Terminal)]  # !?!
        
        # If w is a descendant of v then w appears before v in self.indices.
        self.indices = []
        used = set()
        
        def dfs(v):
            ''' Perform a depth first traversal starting at the given index. '''
            
            for child in self(v):
                if not isinstance(child, Terminal) and child not in used:
                    dfs(child)
            used.add(v)
            self.indices.append(v)
        dfs(0)
        
        self.num_children = [None] * self.size()
        for index in self.indices:
            self.num_children[index] = sum(1 if isinstance(item, Terminal) else self.num_children[item] for item in self(index))
    
    def __str__(self):
        if len(self) <= 7:
            return str(list(self))
        else:
            return '[%s, %s, %s, ..., %s, %s, %s]' % tuple(chain(islice(self, 3), reversed(list(islice(reversed(self), 3)))))
        
    def __repr__(self):
        strn = []
        for index, item in enumerate(self.graph):
            strn.append('%d --> %s' % (index, item))
        
        return '\n'.join(strn)
    
    def __call__(self, item):
        return self.graph[item]
    def __len__(self):
        return self.num_children[0]
    def size(self):
        ''' Return the number of nodes in this graph. '''
        
        return len(self.graph)
    
    def __getitem__(self, value):
        if isinstance(value, slice):
            # TODO: 2) Implement this.
            return NotImplemented
        else:  # We are returning a single item.
            if value >= len(self) or value < -len(self):
                raise IndexError('index out of range')
            if value < 0: value = len(self) + value
            
            index = 0
            while True:
                for image in self(index):
                    if isinstance(image, Terminal):
                        if value == 0:
                            return image.value
                        else:
                            value -= 1
                    else:
                        if self.num_children[image] > value:
                            index = image
                            break
                        else:
                            value -= self.num_children[image]
        
        raise RuntimeError('Should not be able to reach here.')
    
    def __iter__(self):
        todo = [0]
        while todo:
            v = todo.pop()
            if isinstance(v, Terminal):
                yield v.value
            else:  # isinstance(v, curver.IntegerType):
                todo.extend(reversed(self(v)))
        return
    
    def __lshift__(self, index):
        return [[item if isinstance(item, Terminal) else item + index for item in lst] for lst in self.graph]
    def __rshift__(self, index):
        return [[item if isinstance(item, Terminal) else item - index for item in lst] for lst in self.graph]
    
    def __add__(self, other):
        return StraightLineProgram([[1, self.size()+1]] + (self << 1) + (other << self.size()+1))
    def __radd__(self, other):
        return StraightLineProgram([[1, other.size()+1]] + (other << 1) + (self << other.size()+1))
    
    def __mul__(self, other):
        if other == 0: return StraightLineProgram([])
        
        binary = [bool(int(x)) for x in bin(other)[2:]]
        binary_graph = [[i+2, i+2] for i in range(len(binary)-1)]
        binary_expansion = [i+1 for i, b in enumerate(binary) if b]
        
        return StraightLineProgram([binary_expansion] + binary_graph + (self << len(binary_graph)+1))
    def __rmul__(self, other):
        return self * other
    
    def __contains__(self, value):
        return Terminal(value) in self.sinks
    
    def reverse(self):
        ''' Return the StraightLineProgram that returns self[::-1]. '''
        
        return StraightLineProgram([lst[::-1] for lst in self.graph])
    def __reversed__(self):
        return iter(self.reverse())
    
    def map(self, function=lambda x: x):
        ''' Return the StraightLineProgram obtained by mapping the values of this one under the given function. '''
        
        return StraightLineProgram([[Terminal(function(child.value)) if isinstance(child, Terminal) else child for child in children] for children in self.graph])

