#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''The setup script.'''

from setuptools import setup, find_packages

def read(file_name):
    with open(file_name) as open_file:
        return open_file.read()

requirements = [
    'networkx',
    'numpy',
]

test_requirements = [
    'pytest',
    'hypothesis',
    'pylint',
    'pytest',
    'pytest-cov'
]

setup(
    name='curver',
    version='0.0.1',
    description='For calculations in the curve complex',
    long_description=read('README.rst') + '\n\n' + read('HISTORY.rst'),
    author='Mark Bell',
    author_email='mcbell@illinois.edu',
    url='https://github.com/MarkCBell/curver',
    packages=find_packages(),
    package_data = {
        'curver.application': ['icon/*'],
        },
    include_package_data=True,
    install_requires=requirements,
    license='GNU General Public License v3',
    zip_safe=False,
    keywords='curver',
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
    test_suite='tests',
    tests_require=test_requirements,
)
