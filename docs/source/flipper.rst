
Flipper
=======

The underlying framework of curver is based on `flipper <https://pypi.python.org/pypi/flipper>`_.
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
	- The ability to procedurally generate the mapping class group of any surface.

At some point someone(’s graduate student) should port these features back to flipper.

If you are looking for a basic framework for manipulating curves, arcs and mapping classes then you should probably base your code on curver.

