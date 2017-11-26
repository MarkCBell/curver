
Surfaces
========

Curver's :meth:`~curver.load.load` function can automatically build the Lickorish generators [FarbMarg12]_ for any surface.
This consists of four families of Dehn twists :math:`a_i, b_i, c_i, p_i` and a family of half-twists :math:`s_i`.
The :math:`p_i` twists are parallel to :math:`a_i` and are arranged as shown below:

.. image:: ./figures/surface.svg
   :height: 300
   :alt: MCG generators
   :target: _images/surface.svg
   :align: center

Of course this generating set is redundant::

    >>> S = curver.load(2, 5)
    >>> S('(a_0.b_0.c_0.b_1)^10') == S('(s_1.s_2.s_3.s_4)^5')
    True
    >>> S = curver.load(4, 4)
    >>> S('a_0.b_0.c_0.b_1.c_1.b_2.c_2.b_3')**18 == S('s_1.s_2.s_3')**4
    True
    >>> S = curver.load(5, 4)
    >>> h = S('b_3.c_2.b_2.a_2.c_1.b_1.b_2.c_1.c_2.b_2.b_3.c_2.a_2.b_2.c_1.b_1')  # From the proof of Theorem 4.14 [FarbMarg12].
    >>> h * S('a_1') * h**-1 == S('a_3')
    True
    >>> S = curver.load(3, 2)
    >>> S('a_0.b_0.c_0.b_1.c_1.b_2.p_1').order()
    8

As expected, when :math:`g = 0` only the half-twists are provided:

.. image:: ./figures/sphere.svg
   :height: 300
   :alt: MCG generators
   :target: _images/sphere.svg
   :align: center
