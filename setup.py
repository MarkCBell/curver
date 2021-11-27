#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''The setup script.'''

from setuptools import setup, find_packages

requirements = [
    'decorator>=4.3.0',
    'networkx>=2.6',
    'numpy>=1.15.1',
    'realalg>=0.3.1',
]

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='curver',
    version='0.4.1',
    description='For calculations in the curve complex',
    long_description=readme(),
    author='Mark Bell',
    author_email='mcbell@illinois.edu',
    url='https://github.com/MarkCBell/curver',
    packages=find_packages(),
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
