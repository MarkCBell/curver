
''' A module for storing constants that should only ever be computed once.

We need to put these in a special place so that curver can compute these after everything
else has been initialised. '''

import curver

QUASICONVEXITY = 10

