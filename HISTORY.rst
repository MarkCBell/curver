
History
=======

0.4.0 (2020-04-13)
------------------

* Officially removed Python 2 support.
* Disjoint sum supports a dictionay specifying multiplicities.
* Stripping leading and trailing identity isometries for performance.
* Rewritten to assume that laminations are the correct subclass.
* Switched to left_weight and right_weight.
* Removed Coveralls support.
* More figures for documention.
* Using six.
* MCG now use SLP to parse creation strings.
* Fixed bug where twists could occur even when drop=0.

0.3.5 (2019-07-22)
------------------

* Added Pachner 1 --> 3 move.
* Added MultiEdgeFlip Move.
* Crushing now respects the underlying triangulation.
* Using custom Half class.
* Now using create methods to build and link moves and their inverses simultaneously.
* Using IntegralLamination intermediate class.
* Testing curver.load.
* Dynamically building triangulations during testing.
* Memoize works even for methods that raise Exceptions.
* Standardised error message style.
* Fixed appveyor build process.
* Fixed py.test usage.

0.3.4 (2019-03-29)
------------------

* Better documentation including installation instructions.
* Travis now using stages.
* Docstring bugfix.

0.3.3 (2019-03-29)
------------------

* Extra decorators.
* Increased performance via increased memoization.
* Added link to herokuapp.
* Fixed typo in orbit calculations.

0.3.2 (2019-02-09)
------------------

* Extended quotient orbifolds to all finite subgroups.
* Extended the generating sets produced by curver.load().
* Can now identify which mapping classes are Dehn twists together with their twisting curve(s).
* Can now apply a mapping class together with a power.
* Improved EdgeFlip performance.
* Improved simplification performance.
* Simplified use of error classes.
* Fixed bug in scaling Laminations.

0.3.1 (2018-10-10)
------------------

* Now using a quadratic time algorithm to compute invariant multiarcs.
* Fixed bug in computing Mosher sequence.

0.3.0 (2018-10-02)
------------------

* Solves the conjugacy problem for periodic mapping classes.
* Computes quotient orbifold signature of periodic mapping classes.
* Finds invariant polygonalisations of periodic mapping classes.
* Switched to standard version number extraction.
* Switched documentation to use autosummary.
* Separated out Tox testing evnironments.
* Tox checks for print statements in the kernel.
* Mirgrated from NumPy matrices to arrays.

0.2.5 (2018-06-24)
------------------

* Typo hotfix.

0.2.4 (2018-06-24)
------------------

* Fixed API documentation build process and typos.
* Separated out Mappings from Encodings.
* Triangulation.surface now returns a dictionary.
* All moves can now be packaged.
* Standardised unittests.

0.2.3 (2018-04-18)
------------------

* Fixed twist homology action typo.

0.2.2 (2018-04-17)
------------------

* Minimise and shorten now only return encodings.
* Fixed shorten ordering bug.
* Fixed appveyor build process.

0.2.1 (2018-04-12)
------------------

* Generalised intersection code.
* Collated peripheral components code.
* Fixed documentation typos.
* More unittests.

0.2.0 (2018-04-11)
------------------

* New notion of short.
* Removed TrainTracks class, incorporated vertex cycles methods into MultiCurves.
* More performance.
* Simplified hypothesis strategy code.

0.1.2 (2018-02-19)
------------------

* Licence hotfix.

0.1.1 (2018-02-19)
------------------

* Fixed numpy dtype bug.
* Added first version of SLP data structure.
* Switched to MIT licence.
* Using tox-travis.
* Whitespace, formatting, Flake8.

0.1.0 (2017-12-11)
------------------

* First full release.

0.0.1 (2017-12-08)
------------------

* Test release on PyPI.
