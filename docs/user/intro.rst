
Introduction
============

Curver was designed to test out the Bell--Webb algorithm [BellWebb16]_ to determine the Nielsen--Thurston type of a mapping class.
Its development has enabled several optimisations of the original algorithm to be found.
We hope that this will continue to the point at which this algorithm is actually practical to use.

Flipper
-------

The underlying framework of curver is based on `flipper <https://pypi.org/project/flipper/>`_.
This means that curver can do many of the things that flipper can, however, curver has several advantages including:

    - Curves and Arcs are now proper classes.
    - Homology is handled separately.
    - The ability to perform Dehn twists about all (multi)curves.
    - The ability to perform half-twists about all arcs.
    - Fast twists and half-twists.
    - Polynomial-time curve simplification.
    - Component extraction.
    - The ability to crush along (multi-)curves.
    - The ability to determine the topological type of a multicurve.
    - The ability to procedurally generate the mapping class group of any :ref:`surface <Surfaces>`.

At some point someone(’s graduate student) should port these features back to flipper.

If you are looking for a basic framework for manipulating curves, arcs and mapping classes then you should probably base your code on curver.

Curver License
--------------

    .. include:: ../../LICENSE
