
try:
	from setuptools import setup
except ImportError:
	print('Unable to import setuptools, using distutils instead.')
	from distutils.core import setup

setup(
	name='curver',
	version='0.1.0',
	description='For calculations in the curve complex',
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
		'networkx',  # Curve.tight_geodesic
		]
	)

