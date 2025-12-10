.. _installation:

=========================
Installation
=========================

This section describes two type of installation:

- `Installation for standard users`_: for processing and analyze the MASCDB archive.
- `Installation for contributors`_: who want to enrich the project (e.g., adding new descriptors, ...).

We recommend setting up a virtual environment before installing pymascdb.


.. _virtual_environment:

Virtual Environment Creation
============================

Although optional, using a virtual environment when installing pymascdb is recommended.

Virtual environments isolate dependencies, simplify package management, improve maintainability,
enhance security, and streamline your development workflow.

Below are two options for creating a virtual environment,
using `venv <https://docs.python.org/3/library/venv.html>`__ or
`conda <https://docs.conda.io/en/latest/>`__ (recommended).

**With conda:**

* Install `mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_, `miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ or `anaconda <https://docs.anaconda.com/anaconda/install/>`__ if you haven't already installed.
* Create a new conda environment (e.g., *pymascdb-py311*):

.. code-block:: bash

    conda create --name pymascdb-py311 python=3.11 --no-default-packages

* Activate the environment:

.. code-block:: bash

    conda activate pymascdb-py311

**With venv:**

* On Windows, create and activate a virtual environment:

.. code-block:: bash

    python -m venv pymascdb-pyXXX
    cd pymascdb-pyXXX/Scripts
    activate

* On macOS/Linux, create and activate a virtual environment:

.. code-block:: bash

    python3 -m venv pymascdb-pyXXX
    source pymascdb-pyXXX/bin/activate


.. _installation_standard:

Installation for standard users
==================================

The latest pymascdb stable version is available
on the `Python Packaging Index (PyPI) <https://pypi.org/project/pymascdb/>`__
and on the `conda-forge channel <https://anaconda.org/conda-forge/pymascdb>`__.

Therefore you can either install the package with pip or conda (recommended).
Please install the package in the virtual environment you created before !

**With conda:**

.. code-block:: bash

   conda install -c conda-forge pymascdb


.. note::
   In alternative to conda, if you are looking for a lightweight package manager you could use `micromamba <https://micromamba.readthedocs.io/en/latest/>`__.

**With pip:**

.. code-block:: bash

   pip install pymascdb


.. _installation_contributor:

Installation for contributors
================================

The latest pymascdb version is available on the GitHub repository `pymascdb <https://github.com/ltelab/pymascdb>`__.
You can install the package in editable mode, so that you can modify the code and see the changes immediately.
Here below we provide the steps to install the package in editable mode.

Clone the repository from GitHub
......................................

According to the :ref:`contributors guidelines <contributor_guidelines>`, you should first
`create a fork into your personal GitHub account <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`__.

Then create a local copy of the repository you forked with:

.. code-block:: bash

   git clone https://github.com/<your-account>/pymascdb.git
   cd pymascdb

Create the development environment
......................................

We recommend to create a dedicated conda environment for development purposes.
You can create a conda environment (i.e. with python 3.11) with:

.. code-block:: bash

	conda create --name pymascdb-dev-py311 python=3.11 --no-default-packages
	conda activate pymascdb-dev-py311

Install the package dependencies
............................................

.. code-block:: bash

	conda install --only-deps pymascdb


Install the package in editable mode
................................................

Install the pymascdb package in editable mode by executing the following command in the pymascdb repository's root:

.. code-block:: bash

	pip install -e ".[dev]"


Install code quality checks
..............................................

Install the pre-commit hook by executing the following command in the pymascdb repository's root:

.. code-block:: bash

   pre-commit install


Pre-commit hooks are automated scripts that run during each commit to detect basic code quality issues.
If a hook identifies an issue (signified by the pre-commit script exiting with a non-zero status), it halts the commit process and displays the error messages.

.. note::

	The versions of the software used in the pre-commit hooks is specified in the `.pre-commit-config.yaml <https://github.com/ltelab/pymascdb/blob/main/.pre-commit-config.yaml>`__ file. This file serves as a configuration guide, ensuring that the hooks are executed with the correct versions of each tool, thereby maintaining consistency and reliability in the code quality checks.


Further details about pre-commit hooks can be found in the Contributors Guidelines, specifically in the provided in the :ref:`Code quality control <code_quality_control>` section.


Optional dependencies
=======================

Specific functionality in pymascdb may require additional optional dependencies.
To unlock the full functionalities offered by pymascdb, it is recommended to install also the packages detailed here below.

IDEs
..............

For an improved development experience, consider installing the intuitive `Jupyter <https://jupyter.org/>`_ and
`Spyder <https://www.spyder-ide.org/>`_ Python Integrated Development Environments (IDEs):

.. code-block:: bash

   conda install -c conda-forge jupyter spyder


Speed Up Xarray Computations
...............................

To speed up arrays computations with xarray, install
`flox <https://flox.readthedocs.io/en/latest/>`_,
`numbagg <https://github.com/numbagg/numbagg>`_,
`bottleneck <https://bottleneck.readthedocs.io/en/latest/intro.html>`_ and
`opt-einsum <https://optimized-einsum.readthedocs.io/en/stable/>`_:

.. code-block:: bash

   conda install -c conda-forge flox numbagg bottleneck opt-einsum



Analyze MASCDB data on Jupyter Notebooks
==========================================

If you want to run pymascdb on a `Jupyter Notebook <https://jupyter.org/>`__,
you have to take care to set up the IPython kernel environment where pymascdb is installed.

For example, if your conda/virtual environment is named ``pymascdb-dev``, run:

.. code-block:: bash

   python -m ipykernel install --user --name=pymascdb-dev

When you will use the Jupyter Notebook, by clicking on ``Kernel`` and then ``Change Kernel``, you will be able to select the ``pymascdb-dev`` kernel.
