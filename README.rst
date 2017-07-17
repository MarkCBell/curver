
Curver is a program for performing calculations in the curve complex.
Eventually it will implement the Bell-Webb algorithm to determine the Nielsen-Thurston type of a mapping class.

Curver can be run as a Python 2 or Python 3 module.
For even more speed (~25% more) consider running curver with the -O optimise bytecode option.

Installation
============

`Curver <https://pypi.python.org/flipper>`_ is available on the `Python Package
Index <https://pypi.python.org>`_. The preferred method for installing the latest
stable release is to use `pip <http://pip.readthedocs.org/en/latest/installing.html>`_
(included in Python 2.7.9 and Python 3.4 by default)::

	> python -m pip install curver --user --upgrade

.. warning:: In order to use the curver GUI on OS X, users must first update
	their copy of Tk/Tcl as described `here <https://www.python.org/download/mac/tcltk/>`_.
	Curver has been tested with `ActiveTcl 8.5.18 <http://www.activestate.com/activetcl/downloads>`_.

Usage
=====

Once installed, start the curver GUI by using::

	> python -m curver.app


Citing
======

If you find curver useful in your research, please consider citing it. A suggested
reference is::

	Mark Bell. curver (Computer Software).
	pypi.python.org/pypi/curver, 2017. Version <<version number>>

the BibTeX entry::

	@Misc{curver,
		author = {Bell, Mark},
		title = {curver (Computer Software)},
		howpublished = {\url{pypi.python.org/pypi/curver}},
		year = {2017},
		note = {Version <<version number>>}
	}

or the BibItem::

	\bibitem{curver} Mark Bell: \emph{curver (Computer Software)},
		\url{pypi.python.org/pypi/curver}, (2017),
		Version <<version number>>.

Development
===========

The latest development version of curver is available from
`BitBucket <https://bitbucket.org/Mark_Bell/curver>`_.
Alternatively, you can clone the mercurial repository directly using
the command::

	> hg clone https://bitbucket.org/mark_bell/curver

And then install using the command::

	> python setup.py install --user

