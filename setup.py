
import os
try:
	from setuptools import setup
except ImportError:
	print('Unable to import setuptools, using distutils instead.')
	from distutils.core import setup

def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
	name='curver',
	version='0.1.0',
	description='For calculations in the curve complex',
	long_description=read('README'),
	author='Mark Bell',
	author_email='mcbell@illinois.edu',
	url='https://bitbucket.org/Mark_Bell/curver',
	# Remember to update these if the directory structure changes.
	packages=[
		'curver',
		'curver.kernel',
		'curver.application',
		'curver.tests',
		],
	package_data={
		'curver.application': ['icon/*', 'frames/*'],
		},
	install_requires=[
		'hypothesis',  # tests.hyp
		'networkx',  # Curve.tight_geodesic
		],
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Education',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Natural Language :: English',
		'Operating System :: OS Independent',
		'Programming Language :: Python',
		'Topic :: Scientific/Engineering :: Mathematics',
		],
	)

