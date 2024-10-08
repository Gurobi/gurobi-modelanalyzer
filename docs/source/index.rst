.. gurobi-modelanalyzer documentation master file, created by
   sphinx-quickstart on Thu Jun 15 08:54:02 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Gurobi Modelanalyzer
====================
gurobi-modelanalyzer is a Python package containing different modules
designed to provide additional insights and model characteristics in
addition to the functionality available in the Gurobi programming APIs.
The initial module consists of two functions designed to provide
explanations of ill conditioned basis matrices.  Additional modules
to provide detailed information about model data shall be provided in
future versions.

The included solcheck module is used to analyze feasibility of a
provided solution to a model instance.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   illcond
   solcheck
   contactus
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
