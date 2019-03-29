
Installation
============

The first step to using any software package is getting it properly installed.

From PyPI
---------

Curver is available on the `Python Package Index`_.
The preferred method for installing the latest stable release is to use `pip`_ (included in Python 2.7.9+ and Python 3.4+ by default)::

    $ pip install curver

If you don't have Python installed, this `Python installation guide`_ can guide you through the process.
Consider using the ``--upgrade`` and ``--user`` flags to ensure that all required packages are upgraded and that curver is installed into a sensible place.

.. warning::
    In order to use the curver GUI on OS X, users must first update
    their copy of Tk/Tcl as described `here <https://www.python.org/download/mac/tcltk/>`_.
    Curver has been tested with `ActiveTcl 8.5.18 <https://www.activestate.com/activetcl/downloads>`_.

Since curver is under active development, you can install the latest development version via::

    $ pip install git+git://github.com/MarkCBell/curver.git@dev

Again, consider using the ``--upgrade`` and ``--user`` flags.

From source
-----------

Curver is free open source software and so it can be installed directly from its source code.

Obtaining the source
~~~~~~~~~~~~~~~~~~~~

The official git repository of Curver's source code is available on `GitHub <https://github.com/MarkCBell/curver>`_.
You can either clone the public repository::

    $ git clone git://github.com/MarkCBell/curver.git

Or, download the `tarball <https://github.com/MarkCBell/curver/tarball/master>`_::

    $ curl -OL https://github.com/MarkCBell/curver/tarball/master
    # optionally, zipball is also available (for Windows users).

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

Once you have a copy of the source, you can embed it in your own Python package, or install it into your site-packages easily::

    $ cd curver
    $ pip install . --user

.. _Python Package Index: https://pypi.org/project/curver/
.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

