
''' A module of data structures. '''

from collections import defaultdict, namedtuple
from itertools import chain, groupby, islice
try:
    range = xrange
except ImportError:  # Python3.
    pass

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
        if isinstance(data, StraightLineProgram):
            return data
        
        self.graph = [tuple(children) for children in data]
        self.sinks = [item.value for lst in self.graph for item in lst if isinstance(item, Terminal)]  # !?!
        
        # If w is a descendant of v then w appears before v in self.indices.
        self.indices = []
        used = set()
        def dfs(i):
            for child in self(i):
                if not isinstance(child, Terminal) and not child in used:
                    dfs(child)
            used.add(i)
            self.indices.append(i)
        dfs(0)
        
        self.num_children = [None] * self.size()
        for index in self.indices:
            self.num_children[index] = sum(1 if isinstance(item, Terminal) else self.num_children[item] for item in self.graph[index])
    
    @classmethod
    def from_list(cls, data):
        return StraightLineProgram([[Terminal(item) for item in data]])
    
    def __str__(self):
        if len(self) <= 7:
            return str(list(self))
        else:
            return '[%s, %s, %s, ..., %s, %s, %s]' % tuple(chain(islice(self, 3), reversed(list(islice(reversed(self), 3)))))
            # return '[%s, %s, %s, ..., %s, %s, %s]' % (self[0], self[1], self[2], self[-3], self[-2], self[-1])
        
    def __repr__(self):
        return str(self)
        strn = []
        for index, item in enumerate(self.graph):
            strn.append('%d --> %s' % (index, item))
        
        return '\n'.join(strn)
    
    def __call__(self, item):
        return self.graph[item]
    def __len__(self):
        return self.num_children[0]
    def size(self):
        return len(self.graph)
    
    def __lshift__(self, index):
        return [[item if isinstance(item, Terminal) else item + index for item in lst] for lst in self.graph]
    def __rshift__(self, index):
        return [[item if isinstance(item, Terminal) else item - index for item in lst] for lst in self.graph]
    
    def __getitem__(self, value):
        if isinstance(value, slice):
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
    
    def __iter__(self):
        todo = [0]
        while todo:
            v = todo.pop()
            if isinstance(v, Terminal):
                yield v.value
            else:
                todo.extend(reversed(self(v)))
        return

    def __eq__(self, other):
        # TODO: Implement Plandowski's algorithm to run in poly(self.size() + other.size()) instead of poly(len(self) + len(other)).
        return len(self) == len(other) and all(x == y for x, y in zip(self, other))
    
    def __ne__(self, other):
        return not (self == other)
    
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
        return StraightLineProgram([lst[::-1] for lst in self.graph])
    def __reversed__(self):
        return iter(self.reverse())
    
    def count(self, value, root=None):
        counts = [None] * self.size()
        for index in self.indices:
            counts[index] = sum((1 if item.value == value else 0) if isinstance(item, Terminal) else counts[item] for item in self(index))
        return counts[0]

