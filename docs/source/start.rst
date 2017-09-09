
Getting Started
===============

Installation
~~~~~~~~~~~~

`Curver <https://pypi.python.org/curver>`_ is available on the `Python Package Index <https://pypi.python.org>`_.
The preferred method for installing the latest stable release is to use `pip <http://pip.readthedocs.org/en/latest/installing.html>`_ (included in Python 2.7.9 and Python 3.4 by default)::

	> python -m pip install curver --user --upgrade

.. warning::
	In order to use the curver GUI on OS X, users must first update
	their copy of Tk/Tcl as described `here <https://www.python.org/download/mac/tcltk/>`_.
	Curver has been tested with `ActiveTcl 8.5.18 <http://www.activestate.com/activetcl/downloads>`_.


A taste of curver
~~~~~~~~~~~~~~~~~

This sample curver program loads the mapping class group of the twice-punctured torus and computes the images of some arcs and curves under a mapping class::

	import curver
	S = curver.load('S_1_2')  # The twice-punctured torus.
	h = S('abC')  # The monodromy of the Whitehead link.
	a = S.lamination([0, -1, 0, 0, 0, 0])  # An arc on S.
	
	assert(isinstance(a, curver.kernel.Arc))
	
	print(h(a))  # Its image under h.
	
	c = a.boundary()  # The boundary of a regular neighbourhood.
	assert(h(c) == h(a).boundary())
	
	assert(not c.is_filling())  # A single curve cannot fill.
	assert(not c.fills_with(a))  # Even c \cup a does not fill.
	assert((h**4)(c).fills_with(a))  # But h^4(c) \cup a does.
	
	assert(h != h.inverse() and h != h**2)
	assert(h.inverse() == h**-1)
	assert(h.order() == 0)  # Since this has infinite order.
	
	
	twist = c.encode_twist()
	halftwist = a.encode_halftwist()
	assert(twist == halftwist**2)

It's often hard to visualise what is going on on these surfaces.
Fortunately curver can show us these curves (use Ctrl+W to quit)::

	curver.show(c, a, h(a))  # Start the GUI (see above warning).
	
	curver.show([(h**i)(a) for i in range(10)])

Future code
~~~~~~~~~~~

Unfortunately, these features are currently not implemented.
When they are, some of these will move to the taster section::

	import curver
	S = curver.load('S_1_2')
	h = S('abC')
	a = S.lamination([0, -1, 0, 0, 0, 0])
	c = a.boundary()
	
	print(a.fills_with(c))
	
	print(h.asymptotic_translation_length())  # Need crush.
	
	print(h.nielsen_thurston_type())  # Need asymptotic translation length.
	
	S_3_4 = curver.load('(3, 4)')  # Need to port Will Worden's code.

Development
~~~~~~~~~~~

The latest development version of curver is available from `BitBucket <https://bitbucket.org/Mark_Bell/curver>`_.
Alternatively, you can clone the `mercurial <https://www.mercurial-scm.org/>`_ repository directly using the command::

	> hg clone https://bitbucket.org/mark_bell/curver

And then install using the command::

	> python setup.py develop --user

In several places work is flagged TODO: #). The number determines the priority.
	1) == Major feature.
	2) == Performance (incl. polynomial time).
	3) == Minor feature.

