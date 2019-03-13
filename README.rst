
Curver
======

.. image:: https://img.shields.io/pypi/v/curver.svg
    :target: https://pypi.org/project/curver/
    :alt: PyPI version

.. image:: https://img.shields.io/pypi/l/curver.svg
    :target: https://pypi.org/project/curver/
    :alt: PyPI license

.. image:: https://travis-ci.org/MarkCBell/curver.svg?branch=master
    :target: https://travis-ci.org/MarkCBell/curver
    :alt: Travis build status

.. image:: https://ci.appveyor.com/api/projects/status/kd8b36bkas7h9pp6/branch/master?svg=true
    :target: https://ci.appveyor.com/project/MarkCBell/curver/branch/master
    :alt: AppVeyor build status

.. image:: https://readthedocs.org/projects/curver/badge/?version=master
    :target: https://curver.readthedocs.io
    :alt: Documentation status

.. image:: https://img.shields.io/coveralls/github/MarkCBell/curver.svg?branch=master
    :target: https://coveralls.io/github/MarkCBell/curver?branch=master
    :alt: Coveralls status


Curver is a program for performing calculations in the curve complex.
It implements the Bell--Webb algorithm [BellWebb16]_ to determine the Nielsen--Thurston type of a mapping class.
This algorithm runs in polynomial time but the constants involved currently make this implementation impractical.

Curver officially supports Python 2.7 and 3.4 -- 3.7.
It also runs on PyPy and `Sage`_.

.. note:: The use of **Python 3** is *highly* preferred over Python 2.
    Consider upgrading your applications and infrastructure if you find yourself *still* using Python 2 in production today.
    If you are using Python 3, congratulations — you are indeed a person of excellent taste. — *Kenneth Reitz*

A taste of curver::

    >>> S = curver.load(0, 5)
    >>> S('s_0.s_1.s_0') == S('s_1.s_0.s_1')
    True
    >>> S('s_0.s_1.s_2.s_3').order(), S('s_0.s_1.s_3.s_2').order(), S('s_0.s_1.S_2.S_3').order()
    (5, 5, 5)
    >>> S('s_0.s_1.s_2.s_3').is_conjugate_to(S('s_0.s_1.s_3.s_2'))
    True
    >>> S('s_0.s_1.s_2.s_3').is_conjugate_to(S('s_0.s_1.S_2.S_3'))
    False

Features
--------

    - Solves the word problem for mapping class groups.
    - Performs Nielsen--Thurston classification of mapping classes.
    - `Solves the conjugacy problem for periodic mapping classes <https://periodic.herokuapp.com>`_.
    - Computes the asymptotic translation length of mapping classes on the curve complex.
    - Computes geodesics in the curve complex.
    - Computes quotient orbifolds and their quotient maps.
    - Computes the action of mapping classes on H_1.
    - Determines the topological type of multicurves.

.. _Sage: http://www.sagemath.org/
