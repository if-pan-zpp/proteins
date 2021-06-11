Analysis of legacy code
=======================

Sources
-------
The goal is to replicate (whilst possibly fixing bugs, if any are found) the program originally written in Fortran 77 language. More specifically, given were:

- the source code, in the form of a single file ``cg.f`` spanning 8.3k LOC;
- input parameters (``parameters.txt``, ``parametersMD05.txt``, ``parametersMDCG.txt``, ``parametersMJ96.txt``) in a custom format;
- a ``README.txt`` file, which in particular contains:

  - compilation instructions;
  - usage instructions;
  - descriptions of input files and (some) global variables;
  - descriptions of output files;

- three examples, in ascending order of required simulation time and complexity:

  - simulation of ubiquitin, given in PDB format;
  - simulation of ``9AAC``, given in custom format;
  - simulation of gluten;

- a research paper ``CPC14.pdf`` providing the high-level description of the underlying model.

Fortran 95 sources
^^^^^^^^^^^^^^^^^^
Later in development, we acquired an incomplete reimplementation of Fortran 77 program in Fortran 95. The extent to which this implementation is helpful is unclear.

Papers
^^^^^^
Aside from the code itself, we found some research papers utilizing the program to shed more light on certain aspects of the simulation that we found unclear:

- aforementioned ``CPC14.pdf``;
- `PDB reference <http://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf>`_;
- `Disordered peptide chains in an alpha-C-based coarse-grained model <https://arxiv.org/pdf/1807.07825.pdf>`_;
- `Dynamika molekularna białek strukturalnie nieuporządkowanych oraz ich agregatów w ramach modeli gruboziarnistych <http://www.ifpan.edu.pl/rn_ifpan/Mioduszewski-doktorat.pdf>`_.

Documentation
-------------
We found it more appropriate to document the code out-of-source. This is chiefly due to the significant "entanglement" of the reference source code. 

.. toctree::
   :maxdepth: 2
  
   overview
   main
   functions