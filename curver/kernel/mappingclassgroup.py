
''' A module for representing triangulations along with laminations and mapping classes on them. '''

from random import choice
import re

import curver

REGEX_IS_NAME = re.compile(r'[a-z]\w*$')

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
        assert all(REGEX_IS_NAME.match(name) for name in pos_mapping_classes)
        
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
        ''' Return a random word of the required length.
        
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
                raise TypeError('At least one of positive and negative must be allowed.')
        
        return '.'.join(choice(letters) for _ in range(length))
    
    def decompose_word(self, word):
        ''' Return a list of mapping_classes keys whose concatenation is word and the keys are chosen greedly.
        
        Raises a TypeError if the greedy decomposition fails. '''
        
        assert isinstance(word, str)
        
        # By sorting the available keys, longest first, we ensure that any time we
        # get a match it is as long as possible.
        available_letters = sorted(self.mapping_classes, key=len, reverse=True)
        decomposition = []
        for subword in word.split('.'):
            while subword:
                for letter in available_letters:
                    if subword.startswith(letter):
                        decomposition.append(letter)
                        subword = subword[len(letter):]
                        break
                else:
                    raise TypeError('After extracting %s, the remaining %s cannot be greedly decomposed as a concatination of self.mapping_classes.' % ('.'.join(decomposition), subword))
        
        return decomposition
    
    def mapping_class(self, word):
        ''' Return the mapping class corresponding to the given word or a random one of given length if given an integer.
        
        The given word is decomposed using self.decompose_word and the composition
        of the mapping classes involved is returned.
        
        Raises a TypeError if the word does not correspond to a mapping class. '''
        
        if not isinstance(word, str):
            word = self.random_word(word)
        
        # Remove any whitespace.
        word = word.replace(' ', '')
        
        # Check for balanced parentheses.
        counter = 0
        for letter in word:
            if letter == '(': counter += 1
            if letter == ')': counter -= 1
            if counter < 0: raise TypeError('Unbalanced parentheses.')
        if counter != 0: raise TypeError('Unbalanced parentheses.')
        
        # Expand out parenthesis powers.
        # This can fail with a TypeError.
        old_word = None
        while word != old_word:  # While a change was made.
            old_word = word
            for subword, power in re.findall(r'(\([\w_\.]*\))\^(-?\d+)', word):
                decompose = self.decompose_word(subword[1:-1])
                int_power = int(power)
                if int_power > 0:
                    replacement = '.'.join(decompose) * int_power
                else:
                    replacement = '.'.join(letter.swapcase() for letter in decompose[::-1]) * abs(int_power)
                word = word.replace(subword + '^' + power, replacement)
        
        # Remove any remaining parenthesis, these do not have a power and are treated as ^1
        word = word.replace('(', '').replace(')', '')
        
        # Expand out powers without parenthesis. Here we use a greedy algorithm and take the
        # longest mapping class occuring before the power. Note that we only do one pass and so
        # only all pure powers to be expanded once, that is 'aBBB^2^3' is not recognised.
        available_letters = sorted(self.mapping_classes, key=len, reverse=True)
        for letter in available_letters:
            for subword, power in re.findall(r'(%s)\^(-?\d+)' % letter, word):
                int_power = int(power)
                word = word.replace(subword + '^' + power, (letter if int_power > 0 else letter.swapcase()) * abs(int_power))
        
        # This can fail with a TypeError.
        sequence = [item for letter in self.decompose_word(word) for item in self.mapping_classes[letter]]
        return curver.kernel.MappingClass(sequence) if sequence else self.triangulation.id_encoding()
    
    def __call__(self, word):
        ''' A shortcut for self.mapping_class(...). '''
        return self.mapping_class(word)
    
    def lamination(self, geometric):
        ''' Return a new lamination on this surface assigning the specified weight to each edge. '''
        
        return self.triangulation(geometric)

