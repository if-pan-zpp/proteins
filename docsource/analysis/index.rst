Code analysis
=============

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
~~~~~~~~~~~~~~~~~~

Later in development, we acquired an incomplete reimplementation of Fortran 77 program in Fortran 95. The extent to which this implementation is helpful is unclear.

Incremental design and code coverage
------------------------------------

We have elected to provide documentation and implementation of only inasmuch of the source code as to be capable of implementing consecutive examples. To facilitate it, we used an external program to determine which parts of the code are used for each example.

Fortran 77 source code documentation
------------------------------------

.. toctree::
   :maxdepth: 2

   main
   functions