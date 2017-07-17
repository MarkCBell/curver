
import argparse
import curver.application

start = curver.application.start

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Curver GUI')
	parser.add_argument('load', nargs='?', help='path to load from when starting')
	args = parser.parse_args()
	
	start(load_from=args.load)

