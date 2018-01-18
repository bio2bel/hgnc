Bio2BEL HGNC |build| |coverage| |docs|
======================================
This package creates:

- HGNC BEL namespace
- HGNC gene families BEL namespace,
- HGNC gene family membership BEL Script
- orthology BEL script

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

.. |docs| image:: http://readthedocs.org/projects/bio2bel-hgnc/badge/?version=latest
    :target: http://bio2bel.readthedocs.io/projects/hgnc/en/latest/?badge=latest
    :alt: Documentation Status
