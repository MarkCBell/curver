
import unittest
import curver

# Things to test:
#  Relations in the mapping class group.

class TestS_1_1(unittest.TestCase):
    def setUp(self):
        self.S = curver.load('S_1_1')
        self.identity = self.S('')
    def assertRelation(self, word):
        self.assertEqual(self.S(word), self.identity)
    def test_mapping_class(self):
        h = self.S('ab')
        self.assertEqual(h.order(), 6)
        g = self.S('aba')
        self.assertEqual(g.order(), 4)
        involution = self.S('ababab')
        self.assertEqual(involution.order(), 2)
        identity = self.S('(ab)^6')
        self.assertEqual(identity.order(), 1)
    def test_relations(self):
        self.assertRelation('abaABA')
        self.assertRelation('abaBAB')

class TestS_1_2(unittest.TestCase):
    def setUp(self):
        self.S = curver.load('S_1_2')
        self.a = self.S.lamination([4, 3, 1, 2, 1, 3])
        self.b = self.S.lamination([195039444271276, 143576284551183, 259248710713939, 119795766434662, 23780518116521, 96015248318141])
        self.h = self.S('abc')
        self.g = self.S('abcaaabxx')
        self.identity = self.S('')
    def test_group(self):
        self.assertEqual(self.S('a^2b^3'), self.S('aabbb'))
    def test_components(self):
        self.assertEqual(self.a.triangulation.sum([multiplicity * component for component, multiplicity in self.a.components().items()]), self.a)
        self.assertTrue(isinstance(self.b, curver.kernel.Curve))
    def test_mapping_class(self):
        self.assertEqual(self.h, self.h)
        self.assertEqual(self.h.order(), 4)
        self.assertEqual((self.h**(self.h.order())), self.identity)
        self.assertEqual(self.S('xx'), self.S('(ab)^6'))
    def test_images(self):
        self.assertEqual((self.g**0)(self.a), self.S.lamination([4, 3, 1, 2, 1, 3]))
        self.assertEqual((self.g**1)(self.a), self.S.lamination([9, 10, 11, 8, 2, 6]))
        self.assertEqual((self.g**20)(self.a), self.b)
    def test_intersection(self):
        self.assertEqual(self.S('').intersection_matrix(), [[-1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]])

if __name__ == '__main__':
    unittest.main()

