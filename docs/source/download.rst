.. _download:

Download
=======================================
*mascdb* can be downloaded `from github <https://github.com/ltelab/pymascdb/>`_
(a *pip* or *conda* package is expected in the future).

Follow the recommended installation procedure:

1. Clone the repository:

.. code-block:: bash

   git clone https://github.com/ltelab/pymascdb.git
   cd pymascdb

2. Create the appropriate *conda* environment (*mascdb*)

.. code-block:: bash

   conda env create -f environment.yml

3. Activate the environment and install the package using the local installer

.. code-block:: bash

   conda activate mascdb
   python setup.py install

In this way, *mascdb* can simply be loaded, from any directory, as:

.. code-block:: python

    import mascdb.api
    from mascdb.api import MASC_DB
