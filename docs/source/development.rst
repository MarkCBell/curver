
Development
~~~~~~~~~~~

The latest development version of curver is available from `BitBucket <https://bitbucket.org/Mark_Bell/curver>`_.
Alternatively, you can clone the `mercurial <https://www.mercurial-scm.org/>`_ repository directly using the command:

.. code-block:: shell

    > hg clone https://bitbucket.org/mark_bell/curver

And then install using the command:

.. code-block:: shell

    > python -m pip install --editable . --user

The packages required for development, including building this documentation, can be installed using the included ``requirements.txt`` file:

.. code-block:: shell

    > python -m pip install -r requirements.txt

Curver includes unittests to aid with development.
These can be run using `tox <https://tox.readthedocs.io/>`_:

.. code-block:: shell

    > tox
    ...
    congratulations :)

Curver uses `Sphinx <http://www.sphinx-doc.org/>`_ to build its documentation.
This is automatically built by readthedocs when changes are pushed to the repository but can also be built locally by using:

.. code-block:: shell

    > cd ./docs
    > make html
    > open ./build/html/index.html

In several places work is flagged TODO: #). The number determines the priority.
    1) == Major feature.
    2) == Performance.
    3) == Minor feature.

