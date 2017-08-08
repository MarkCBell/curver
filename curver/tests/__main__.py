
import os
import unittest

def main():
	loader = unittest.TestLoader()
	start_dir = os.path.dirname(__file__)
	suite = loader.discover(start_dir)
	
	runner = unittest.TextTestRunner()
	print('Running unit tests:')
	runner.run(suite)

if __name__ == '__main__':
	main()

