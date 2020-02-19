
Curver
======

.. image:: https://img.shields.io/pypi/v/curver.svg
    :target: https://pypi.org/project/curver/
    :alt: PyPI version

.. image:: https://img.shields.io/pypi/l/curver.svg
    :target: https://pypi.org/project/curver/
    :alt: PyPI license

.. image:: https://api.travis-ci.com/MarkCBell/curver.svg?branch=master
    :target: https://travis-ci.com/MarkCBell/curver
    :alt: Travis build status

Curver is a program for performing calculations in the curve complex.
It implements the Bell--Webb algorithm to determine the Nielsen--Thurston type of a mapping class.
This algorithm runs in polynomial time but the constants involved currently make this implementation impractical.

Curver officially supports Python 3.6 -- 3.8.
Unoffically, it also runs on `PyPy`_ and `Sage`_ with some care.

Quickstart
----------

Curver is available on `PyPI`_, so it can be installed via::

    $ pip install curver --user --upgrade

Once installed, try it inside of Python::

    >>> import curver
    >>> S = curver.load(0, 5)
    >>> S('s_0.s_1.s_0') == S('s_1.s_0.s_1')
    True
    >>> f = S('s_0.s_1.s_2.s_3')
    >>> g = S('s_0.s_1.s_3.s_2')
    >>> h = S('s_0.s_1.S_2.S_3')
    >>> f.order(), g.order(), h.order()
    (5, 5, 5)
    >>> f.is_conjugate_to(g)
    True
    >>> f.is_conjugate_to(g)
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

External Links
--------------

* `PyPI`_
* `ReadTheDocs`_
* `GitHub`_
* `Travis`_
* `AppVeyor`_
* `Azure`_

.. _AppVeyor: https://ci.appveyor.com/project/MarkCBell/curver
.. _Azure: https://dev.azure.com/MarkCBell/curver
.. _GitHub: https://github.com/MarkCBell/curver
.. _PyPI: https://pypi.org/project/curver
.. _ReadTheDocs: http://curver.readthedocs.io
.. _Sage: http://www.sagemath.org
.. _Travis: https://travis-ci.com/MarkCBell/curver
.. _PyPy: https://pypy.org/

