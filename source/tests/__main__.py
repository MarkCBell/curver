
import unittest
import curver

class TestS_1_1(unittest.TestCase):
	def setUp(self):
		self.S = curver.load('S_1_1')
	def test_mapping_class(self):
		h = self.S('ab')
		self.assertEqual(h.order(), 6)
		g = self.S('aba')
		self.assertEqual(g.order(), 4)
		i = self.S('ababab')
		self.assertEqual(i.order(), 2)

class TestS_1_2(unittest.TestCase):
	def setUp(self):
		self.S = curver.load('S_1_2')
		self.a = self.S.lamination([17, 44, 12, 12, 41, 18])
		self.b = self.S.lamination([5719512847871531, 2642919316191732, 524642836301528, 2118276479890152, 5194870011570003, 6611446870756195])
		self.h = self.S('abc')
		self.g = self.S('abcaaabxx')
		self.identity = self.S('')
	def test_group(self):
		self.assertEqual(self.S('a^2b^3'), self.S('aabbb'))
	def test_mcomponents(self):
		self.assertEqual(sum(mult*comp for comp, mult in self.a.mcomponents()), self.a)
		self.assertTrue(isinstance(self.b, curver.kernel.MultiArc))
	def test_mapping_class(self):
		self.assertEqual(self.h, self.h)
		self.assertEqual(self.h.order(), 4)
		self.assertEqual((self.h**(self.h.order())), self.identity)
	def test_images(self):
		self.assertEqual((self.g**0)(self.a), self.S.lamination([17, 32, 0, 12, 29, 6]))
		self.assertEqual((self.g**1)(self.a), self.S.lamination([213, 143, 32, 59, 181, 216]))
		self.assertEqual((self.g**20)(self.a), self.b)
	def test_intersection(self):
		self.assertEqual(self.S('').intersection_matrix() == [[-1, 0, 0, 0, 0, 0], [0, -1, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, -1, 0], [0, 0, 0, 0, 0, -1]])

if __name__ == '__main__':
	unittest.main()

