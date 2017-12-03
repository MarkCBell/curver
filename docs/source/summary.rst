
.. image:: https://img.shields.io/pypi/v/curver.svg
    :target: https://pypi.python.org/pypi/curver
    :alt: PyPI version

.. image:: https://img.shields.io/pypi/l/curver.svg
    :target: https://pypi.python.org/pypi/curver
    :alt: PyPI license

.. image:: https://img.shields.io/travis/MarkCBell/curver.svg
    :target: https://travis-ci.org/MarkCBell/curver
    :alt: Curver build status

.. image:: https://img.shields.io/readthedocs/curver.svg
    :target: https://curver.readthedocs.io
    :alt: Documentation status

.. image:: https://img.shields.io/coveralls/github/MarkCBell/curver.svg
    :target: https://coveralls.io/github/MarkCBell/curver
    :alt: Coveralls status


Curver is a program for performing calculations in the curve complex.
It implements the Bell--Webb algorithm [BellWebb16]_ to determine the Nielsen--Thurston type of a mapping class.
This algorithm runs in polynomial time but the constants involved currently make this implementation impractical.

.. note:: The use of **Python 3** is *highly* preferred over Python 2.
    Consider upgrading your applications and infrastructure if you find yourself *still* using Python 2 in production today.
    If you are using Python 3, congratulations — you are indeed a person of excellent taste. — *Kenneth Reitz*

