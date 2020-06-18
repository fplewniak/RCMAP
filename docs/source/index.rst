.. RCMAP documentation master file, created by
   sphinx-quickstart on Thu Jun 18 13:44:45 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RCMAP's documentation!
=================================
Residue Conservation in Multiple Alignment of Proteins (RCMAP) is a Python package to help manual annotation of protein
sequences by comparing them to a multiple alignment of reference sequences belonging to a functional family.

The RCMAP package provides the shell command ``evaluate_seq`` whose input is a multiple alignment file in FastA format,
containing reference sequences and one or more unknown sequences to annotate. It then displays for each unknown sequence
whether it is consistent at every user-specified position with the aminoacid conservation profile of the reference
sequences.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started.rst
   api/index.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
