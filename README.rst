Bio2BEL HGNC |build| |coverage| |docs| |zenodo|
===============================================
This package creates:

- HGNC BEL namespace
- HGNC gene families BEL namespace,
- HGNC gene family membership BEL Script
- orthology BEL script

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
``bio2bel_hgnc`` can be installed easily from `PyPI <https://pypi.python.org/pypi/bio2bel_hgnc>`_ with the
following code in your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install bio2bel_hgnc

or from the latest code on `GitHub <https://github.com/bio2bel/hgnc>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/bio2bel/hgnc.git@master

Setup
-----
HGNC can be downloaded and populated from either the Python REPL or the automatically installed command line
utility.

Python REPL
~~~~~~~~~~~
.. code-block:: python

    >>> import bio2bel_hgnc
    >>> hgnc_manager = bio2bel_hgnc.Manager()
    >>> hgnc_manager.populate()

Command Line Utility
~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

    bio2bel_hgnc populate

Citation
--------
Gray KA, Yates B, Seal RL, Wright MW, Bruford EA. genenames.org: the HGNC resources in 2015. Nucleic Acids Res. 2015
Jan;43(Database issue):D1079-85. doi: 10.1093/nar/gku1071. http://www.ncbi.nlm.nih.gov/pubmed/25361968

Abstract
--------
The HUGO Gene Nomenclature Committee (HGNC) based at the European Bioinformatics Institute (EMBL-EBI) assigns unique
symbols and names to human genes. To date the HGNC have assigned over 39,000 gene names and, representing an increase
of over 5000 entries in the past two years. As well as increasing the size of our database, we have continued
redesigning our website http://www.genenames.org and have modified, updated and improved many aspects of the site
including a faster and more powerful search, a vastly improved HCOP tool and a REST service to increase the number of
ways users can retrieve our data. This article provides an overview of our current online data and resources, and
highlights the changes we have made in recent years.

Installation
------------
:code:`pip3 install git+https://github.com/bio2bel/hgnc.git`

Acknowledgements
----------------
- This package heavily relies on Andrej Konotopez's package `PyHGNC <https://github.com/lekono/pyhgnc>`_

.. |build| image:: https://travis-ci.org/bio2bel/hgnc.svg?branch=master
    :target: https://travis-ci.org/bio2bel/hgnc
    :alt: Build Status

.. |coverage| image:: https://codecov.io/gh/bio2bel/hgnc/coverage.svg?branch=master
    :target: https://codecov.io/gh/bio2bel/hgnc?branch=master
    :alt: Coverage Status

.. |documentation| image:: http://readthedocs.org/projects/bio2bel-hgnc/badge/?version=latest
    :target: http://bio2bel.readthedocs.io/projects/hgnc/en/latest/?badge=latest
    :alt: Documentation Status

.. |climate| image:: https://codeclimate.com/github/bio2bel/hgnc/badges/gpa.svg
    :target: https://codeclimate.com/github/bio2bel/hgnc
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/bio2bel_hgnc.svg
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/bio2bel_hgnc.svg
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/bio2bel_hgnc.svg
    :alt: MIT License

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1162644.svg
    :target: https://doi.org/10.5281/zenodo.1162644
