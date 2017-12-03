
Development
~~~~~~~~~~~

The latest development version of curver is available from `GitHub <https://github.com/MarkCBell/curver>`_ via:

.. code-block:: shell

    > git clone https://github.com/MarkCBell/curver

Alternatively, you can clone the mercurial repository from `BitBucket <https://bitbucket.org/Mark_Bell/curver>`_ via:

.. code-block:: shell

    > hg clone https://bitbucket.org/mark_bell/curver

And then install using the command:

.. code-block:: shell

    > pip install --editable . --user

The packages required for development, including building this documentation, can be installed using the included ``requirements-dev.txt`` file:

.. code-block:: shell

    > pip install -r requirements-dev.txt

Curver includes unittests to aid with development.
These can be run using `tox <https://tox.readthedocs.io/>`_:

.. code-block:: shell

    > tox
    ...
    congratulations :)

This includes `pylint <https://www.pylint.org/>`_ which can be run on its own via ``tox -e lint``.
As part of its analysis pylint will highlight places in curves where work is flagged as ``TODO: #)``.
These numbers roughly correspond to:

1. Major feature.
2. Performance.
3. Minor feature.

Curver uses `Sphinx <http://www.sphinx-doc.org/>`_ to build its documentation.
This is automatically built by readthedocs when changes are pushed to the repository but can also be built locally by using:

.. code-block:: shell

    > cd ./docs
    > make html


