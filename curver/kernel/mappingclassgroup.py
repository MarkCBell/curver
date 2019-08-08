
''' A module for representing triangulations along with laminations and mapping classes on them. '''

try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence
from random import choice
import re
from six.moves.queue import Queue

import curver

IS_NAME = re.compile(r'[a-z]\w*$')
IS_INT = re.compile(r'[+-]?\d+')
TOKENS = re.compile(r'[+-]?\d+|\^|\(|\)|[\w_\.]+')
LEADING_DOTS = re.compile(r'^\.*')

class MappingClassGroup(object):
    ''' This represents a triangulation along with a collection of named mapping classes on it.
    
    It can also be given a list of curves and arcs, in which case the twists and half-twists about
    these are also added to the list of known mapping classes.
    Most importantly this object can construct a mapping class from a string descriptor.
    See self.mapping_class for additional information. '''
    def __init__(self, pos_mapping_classes=None, curves=None, arcs=None):
        if pos_mapping_classes is None: pos_mapping_classes = dict()
        if curves is None: curves = dict()
        if arcs is None: arcs = dict()
        
        for name, arc in arcs.items():
            assert name not in pos_mapping_classes
            pos_mapping_classes[name] = arc.encode_halftwist()
        for name, curve in curves.items():
            assert name not in pos_mapping_classes
            pos_mapping_classes[name] = curve.encode_twist()
        
        assert pos_mapping_classes
        
        if not isinstance(pos_mapping_classes, dict):
            pos_mapping_classes = dict(curver.kernel.utilities.name_objects(pos_mapping_classes))
        
        self.triangulation = list(pos_mapping_classes.values())[0].source_triangulation
        self.zeta = self.triangulation.zeta
        
        assert all(isinstance(key, str) for key in pos_mapping_classes)
        assert all(isinstance(pos_mapping_class, curver.kernel.MappingClass) for pos_mapping_class in pos_mapping_classes.values())
        assert all(pos_mapping_class.source_triangulation == self.triangulation for pos_mapping_class in pos_mapping_classes.values())
        assert all(IS_NAME.match(name) for name in pos_mapping_classes)
        
        self.pos_mapping_classes = dict(pos_mapping_classes)
        self.neg_mapping_classes = dict((name.swapcase(), pos_mapping_class.inverse()) for name, pos_mapping_class in self.pos_mapping_classes.items())
        self.mapping_classes = dict(list(self.pos_mapping_classes.items()) + list(self.neg_mapping_classes.items()))
        
        self.arcs = arcs
        self.curves = curves
    
    def __repr__(self):
        return str(self)
    def __str__(self):
        pos_keys = sorted(self.pos_mapping_classes.keys(), key=curver.kernel.utilities.alphanum_key)
        return 'Mapping class group < %s > on %s' % (', '.join(pos_keys), self.triangulation)
    
    def __eq__(self, other):
        return self.triangulation == other.triangulation and self.mapping_classes == other.mapping_classes
    def __ne__(self, other):
        return not self == other
    
    def __getitem__(self, item):
        return self.mapping_classes[item]
    def __iter__(self):
        return iter(self.mapping_classes)
    
    def random_word(self, length, positive=True, negative=True, letters=None):
        ''' Return a random sequence of generators of the required length.
        
        The letters to choose from can be specified or, alternatively, the set
        of positive, negative or all (default) mapping classes can be used by using the
        flags postive and negative. '''
        
        if letters is None:
            pos_keys = sorted(self.pos_mapping_classes.keys())
            neg_keys = sorted(self.neg_mapping_classes.keys())
            
            if positive and negative:
                letters = pos_keys + neg_keys
            elif positive and not negative:
                letters = pos_keys
            elif not positive and negative:
                letters = neg_keys
            else:
                raise TypeError('At least one of positive and negative must be allowed')
        
        return [choice(letters) for _ in range(length)]
    
    def mapping_class(self, data, **kwargs):
        ''' Return a mapping class from data.
        
        Data can either be:
         * an iterable of mapping_class names,
         * an integer specifying the word length of a random mapping class, or
         * a string specifying the generators to be composed together.
        
        The string supports '^' powers, parentheses and optional '.' separators.
        Raises a ValueError if given a string that cannot be decomposed. '''
        
        if isinstance(data, curver.IntegerType):
            sequence = self.random_word(data, **kwargs)
        elif isinstance(data, str):
            SLP = curver.kernel.SLP  # Shorter alias.
            word = '(' + data + ')'  # This ensures that the last token is a ')' and so avoids a special case.
            word = re.sub(r'\s', '', word)  # Remove whitespace.
            
            MATCH_MCs = re.compile('|'.join(sorted(self.mapping_classes, key=len, reverse=True)))
            
            def decompose(word):
                ''' Break a word into a list of words that match MATCH_MCs. '''
                
                iterable = [x for subword in word.split('.') for x in MATCH_MCs.findall(subword)]
                if sum(len(x) for x in iterable) + word.count('.') < len(word):  # We were unable to decompose the entire word.
                    remaining = word
                    for index, item in enumerate(iterable):
                        remaining = LEADING_DOTS.sub('', remaining)  # Remove leading dots.
                        if remaining.startswith(item):
                            remaining = remaining[len(item):]
                        else:
                            raise ValueError('After extracting {}, the remaining {} of {} could not be decomposed'.format(iterable[:index], remaining, word))
                    remaining = LEADING_DOTS.sub('', remaining)  # Remove leading dots.
                    raise ValueError('After extracting {}, the remaining "{}" of "{}" could not be decomposed'.format(iterable, remaining, word))
                return iterable
            
            stack = [[]]
            tokens = TOKENS.findall(word)  # Break word into tokens.
            for index, token in enumerate(tokens):
                if IS_INT.match(token):
                    pass
                elif token == '^':
                    try:
                        power = int(tokens[index+1])
                    except (ValueError, IndexError):
                        raise ValueError('^ not followed by a power')
                    stack[-1][-1] = stack[-1][-1] * abs(power)
                    if power < 0: stack[-1][-1] = stack[-1][-1].reverse().map(lambda x: x.swapcase())
                elif token == '(':
                    stack[-1].append([])
                    stack.append(stack[-1][-1])
                elif token == ')':
                    stack.pop()
                    if not stack:
                        raise ValueError('Unbalanced parentheses')
                    stack[-1][-1] = SLP.sum(stack[-1][-1])
                else:
                    stack[-1].append(SLP(decompose(token)))
            if len(stack) > 1:
                raise ValueError('Unbalanced parentheses')
            
            sequence = stack[-1][-1]
        elif isinstance(data, Sequence):
            sequence = data
        else:
            raise TypeError('No method for generating a Sequence from this type')
        
        sequence = [item for letter in sequence for item in self.mapping_classes[letter]]
        return curver.kernel.MappingClass(sequence) if sequence else self.triangulation.id_encoding()
    
    def __call__(self, word, **kwargs):
        ''' A shortcut for self.mapping_class(...). '''
        return self.mapping_class(word, **kwargs)
    
    def lamination(self, geometric):
        ''' Return a new lamination on this surface assigning the specified weight to each edge. '''
        
        return self.triangulation(geometric)
    
    def cayley(self, generators, length):
        ''' Explore the Cayley graph for self with respect to the given generators.
        Yield the canonical names for all elements with at most the given length as they are encountered. '''
        
        identity = tuple()
        id_mapping_class = self('')
        convert = lambda X: (X[0], tuple(X[1].flatten()))  # Since numpy.ndarrays are not hashable we need a converter.
        elements = {convert((id_mapping_class.source_triangulation.as_lamination(), id_mapping_class.homology_matrix())): identity}
        good = set([identity])
        Q = Queue()
        Q.put(((id_mapping_class.source_triangulation.as_lamination(), id_mapping_class.homology_matrix()), identity))
        good = set([identity])
        yield identity
        while not Q.empty():
            image, word = Q.get()
            
            for generator in generators:
                next_word = (generator,) + word
                # Check all prefixes are good.
                if any(next_word[:i] not in good for i in range(1, len(next_word))):
                    continue
                
                lam, mat = image
                action = self.mapping_classes[generator]
                next_image = (action(lam), action.homology_matrix().dot(mat))
                key = convert(next_image)
                if key not in elements:
                    yield next_word
                    good.add(next_word)
                    elements[key] = next_word
                    if len(next_word) < length:
                        Q.put((next_image, next_word))

