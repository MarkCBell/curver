
Installation
============

The first step to using any software package is getting it properly installed.

From PyPI
---------

Curver is available on the `Python Package Index`_.
The preferred method for installing the latest stable release is to use `pip`_ (included in Python 2.7.9+ and Python 3.4+ by default)::

    $ pip install curver --user --upgrade

If you don't have Python installed, this `Python installation guide`_ can guide you through the process.

.. warning::
    In order to use the curver GUI on OS X, users must first update
    their copy of Tk/Tcl as described `here <https://www.python.org/download/mac/tcltk/>`_.
    Curver has been tested with `ActiveTcl 8.5.18 <http://www.activestate.com/activetcl/downloads>`_.

From source
-----------

Curver is under active development and there are several ways to obtain its source code.
As well as PyPI (and its mirrors), GitHub and BitBucket are the official distribution sources; alternatives are not supported.

Via git
~~~~~~~

A git repository of Curver's source code is available  on `GitHub <https://github.com/MarkCBell/curver>`_.
You can either clone the public repository::

    $ git clone git://github.com/MarkCBell/curver.git

Or, download the `tarball <https://github.com/MarkCBell/curver/tarball/master>`_::

    $ curl -OL https://github.com/MarkCBell/curver/tarball/master
    # optionally, zipball is also available (for Windows users).

Via mercurial
~~~~~~~~~~~~~

A mercurial repository of Curver's source code is available  on `BitBucket <https://bitbucket.org/Mark_Bell/curver>`_.
You can either clone the public repository::

You can clone the public repository::

    $ hg clone https://bitbucket.org/Mark_Bell/curver

Installing from source
~~~~~~~~~~~~~~~~~~~~~~

Once you have a copy of the source, you can embed it in your own Python
package, or install it into your site-packages easily::

    $ cd curver
    $ pip install . --user

.. _Python Package Index: https://pypi.python.org/pypi/curver
.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
