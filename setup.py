
from __future__ import print_function

import os
try:
	from setuptools import setup
except ImportError:
	print('Unable to import setuptools, using distutils instead.')
	from distutils.core import setup

this_directory = os.path.dirname(__file__)
source_directory = os.path.join(this_directory, 'source')
exec(open(os.path.join(source_directory, 'version.py')).read())  # Load in the variable __version__.

setup(
	name='curver',
	version=__version__,
	description='For calculations in the curve complex',
	author='Mark Bell',
	author_email='mcbell@illinois.edu',
	url='https://bitbucket.org/Mark_Bell/curver',
	# Remember to update these if the directory structure changes.
	packages=[
		'curver',
		'curver.kernel',
		'curver.application',
		],
	package_dir={'curver': source_directory},
	package_data={
		'curver.application': ['icon/*', 'frames/*'],
		},
	install_requires=[
		'networkx',
		]
	)

