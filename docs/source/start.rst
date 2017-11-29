
Getting Started
===============

Installation
~~~~~~~~~~~~

`Curver <https://pypi.python.org/pypi/curver>`_ is available on the `Python Package Index <https://pypi.python.org>`_.
The preferred method for installing the latest stable release is to use `pip <http://pip.readthedocs.org/en/latest/installing.html>`_ (included in Python 2.7.9+ and Python 3.4+ by default):

.. code-block:: shell

    > pip install curver --user --upgrade

.. warning::
    In order to use the curver GUI on OS X, users must first update
    their copy of Tk/Tcl as described `here <https://www.python.org/download/mac/tcltk/>`_.
    Curver has been tested with `ActiveTcl 8.5.18 <http://www.activestate.com/activetcl/downloads>`_.


A taste of curver
~~~~~~~~~~~~~~~~~

This sample curver program loads the mapping class group of the twice-punctured torus and computes the images of some arcs and curves under a mapping class::

    >>> import curver
    >>> S = curver.load(1, 2)  # The twice-punctured torus.
    >>> print(S)
    Mapping class group < a_0, b_0, p_1, s_0, s_1 > on 6-WKSv
    
    >>> a = S.lamination([0, 0, 0, -1, 0, 0])  # An arc on S.
    >>> print(a)
    Arc [0, 0, 0, -1, 0, 0] on 6-WKSv
    
    >>> h = S('a_0.b_0.P_1')  # The monodromy of the Whitehead link.
    >>> print(h(a))  # The image of a under h.
    Arc [1, 1, 0, 0, 0, 0] on 6-WKSv
    
    >>> c = a.boundary()  # The boundary of a regular neighbourhood.
    >>> print(c)
    Curve [2, 2, 2, 0, 2, 2] on 6-WKSv
    >>> h(c) == h(a).boundary()
    True
    >>> c.is_filling()  # A single curve cannot fill S_{1,2}.
    False
    >>> c.fills_with(a)  # Even c \cup a does not fill.
    False
    >>> (h**4)(c).fills_with(a)  # But h^4(c) \cup a does.
    True
    
    >>> # Higher order functions on mapping classes.
    >>> h != h.inverse() and h != h**2
    True
    >>> h.inverse() == h**-1
    True
    >>> h.order()  # 0 == infinite order.
    0
    
    >>> twist = c.encode_twist()  # Build new mapping classes:
    >>> halftwist = a.encode_halftwist()
    >>> twist == halftwist**2
    True

It's often hard to visualise what is going on on these surfaces.
Fortunately curver can show us these curves (use Ctrl+W to quit)::

    >>> curver.show(c, a, h(a))  # Start the GUI (see above warning).
    
    >>> curver.show([(h**i)(a) for i in range(10)])

Curver can automatically create the mapping class group of any punctured surface.
It also includes all of the flipper / Twister example surfaces::

    >>> S = curver.load(3, 17)  # Genus 3 with 17 punctures.
    >>> S = curver.load('S_1_2')  # Flipper's twice-punctured torus.

Curver can also perform some hard calculations with mapping classes.
These run in polynomial time but, thanks to large constants, can still take a *very* long time::

    >>> print(h.nielsen_thurston_type())
    >>> print(h.asymptotic_translation_length())

Future code
~~~~~~~~~~~

Unfortunately, these features are currently not implemented.
When they are, some of these will move to the taster section::

    >>> print(a.fills_with(c))

