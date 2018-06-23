#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''The setup script.'''

from setuptools import setup, find_packages

requirements = [
    'networkx>=2.0',
    'numpy',
]

setup(
    name='curver',
    version='0.2.4',
    description='For calculations in the curve complex',
    long_description='See http://curver.readthedocs.io for the full README, LICENCE and documentation.',
    author='Mark Bell',
    author_email='mcbell@illinois.edu',
    url='https://github.com/MarkCBell/curver',
    packages=find_packages(),
    package_data={
        'curver.application': ['icon/*'],
        },
    include_package_data=True,
    install_requires=requirements,
    license='MIT License',
    zip_safe=False,
    keywords='curver',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics',
        ],
)
