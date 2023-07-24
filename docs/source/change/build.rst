#####################
Building and Testing
#####################

Here, you will find all the dependencies and the steps to build, test, and install! 

**************
Unit Testing
**************

During contribution, we would like to ensure that all our tests pass or add new
tests! To do this, we use ``pytest`` internally, and our testing codes rely
on it for mocking and parametrization.

Dependencies
-------------

You will need to install the following packages:

* pytest
* pytest-cov
* pytest-mock
* coverage

Running Tests
--------------

To run the tests, we recommend installing the package in ``develop`` mode first by running
the following codes:

.. code-block:: shell 

    python setup.py develop #Only need to run once
    pytest . --cov

Note that this approach allows you to run the tests as a package instead of from local files,
which is advantagous. The ``--cov`` option gives you a coverage report of our tests.

If you don't want to do the installation, you can run the following:

.. code-block:: shell 

    pytest -m pytest --cov PythonTemplate

---------


***************
Documentation
***************

Our documentation is automatically built by ReadTheDocs. However, if you inclined to build locally
or contribute to our docs, you're more than welcomed to do so!

Dependencies
-------------

You will need the following packages:

* sphinx
* sphinx-rtd-theme
* sphinx-git
* sphinxcontrib-autoprogram
* sphinx-autodoc-typehints

.. note::

    ``sphinx-git`` is available on ``PyPI`` only. If you use a ``conda`` environment, use ``pip`` to install.

Build Documentation
--------------------

You can simply run the following codes to build the documentation locally:

.. code-block:: shell

    cd docs
    make html

The resulting HTML will be in the ``build`` folder, and click on ``index.html``
and you're all set to view the documentation. We don't have unit tests to test the
documentation.