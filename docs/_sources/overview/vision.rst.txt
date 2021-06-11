Vision (Nov 2020)
=================

Note: a more in-depth version of this document is available (in Polish).

Introduction
------------
This document describes the scope of the project of rewriting a protein simulation program based on a model created at IF PAN.

Definitions
^^^^^^^^^^^
- Model - a model of the dynamics of proteins created at IF PAN;
- Reference implementation - an existing implementation of the simulation based on the model, namely a single source file written in Fortran77;
- Rewrite - a program, whose implementation is the end goal of this project;
- *Pseudo-improper-dihedral* (PID) and quasi-adiabatic potentials - two of the available potentials for governing the dynamics of the unstructured parts of proteins.

Context
-------

Description of the product
^^^^^^^^^^^^^^^^^^^^^^^^^^
The program simulates structured, unstructured and partially structured proteins in a coarse-grained model. The reference implementation has been used in many research papers.

Functionality
^^^^^^^^^^^^^
The program simulates proteins in a coarse-grained model, i.e. entire residues are represented as single beads. For structured proteins, Go model for the potential is used, whereas for unstructured parts either PID or quasi-adiabatic potential is used.

- PID potential is more accurate, and more computationally expensive, wherefore it is chiefly used in the simulation of smaller systems;
- quasi-adiabatic potential is less accurate and less computationally expensive, recommended for the simulation of larger structures.

Additional features include performing the simulation on a 2D plane, introduction of repulsive or attractive walls, simulation of ATF stretching the proteins, inclusion of periodic boundary conditions. Also, restart files, contact maps, statistics and current geometry of the system can be periodically saved. The geometry is saved in a standard format, making it possible to view it by, for example, VMD.

The program is on an MIT license.

Drawbacks
^^^^^^^^^
Due to the obsolescence of Fortran77, long lifetime of the program and its development by many parties made its code hard to comprehend. The complexity of the code and a high degree of interdependence of its parts make modifying the program very difficult. Moreover, the degree of the parallelization of the program is insufficient -- at best, 4 CPU cores can be used to speed it up.

Quality of the source code
""""""""""""""""""""""""""
In its current form, the source code of the program is near-incomprehensible. The variable names, due to the need for their brevity caused by the hard limit on the number of columns in Fortran77, in most cases do not help understand their purpose. Some subroutines are extremely long, and most of the code is uncommented. Due to Fortran77, in many circumstances `goto` instruction must be used, which further obfuscates the code. Last, some fragments of the code have been commented out by the authors.

Alternatives
^^^^^^^^^^^^
There are many programs for simulating proteins, all of which approximate the true dynamics of the system in different ways, making direct comparisons between them difficult -- in short, they server different purposes.

Moreover, neither of the other programs includes PID and quasi-adiabatic potentials to be used. Although some of them allow for the introduction of custom forces, none of them provide a mechanism for creating contacts dynamically.

The goal of the project
-----------------------
The goal of the project is to create a rewrite of the reference implementation in C++ programming language, whilst preserving the correctness of the results. The rewrite should not be *significantly* slower than the reference implementation, and should include the possibility of parallelizing the computation with OpenMP. Optional goals include the implementation on MPI and with GPUs.

Benefits
^^^^^^^^
The rewrite is to be written using a modern programming language and modern development techniques, which should make its use and modification easier. Moreover, by unit testing the code, the results shall be more trustworthy, and the modification of the code less error-prone. Improved parallelization capabilities will allow for a reduction in simulation times, and for the simulation of even larger systems.

Requirements
------------

Key requirements
^^^^^^^^^^^^^^^^
In the order of decreasing priority:

1. Correctness of the results -- for all systems, the reference implementation and the rewrite should produce statistically equivalent trajectories, i.e. expected simulation statistics over many samples of the program should be equivalent;
2. Modularity and flexibility of the source code -- including full test coverage;
3. Readability of the code -- including full documentation, instruction for running a simulation, and descriptive error messages;
4. Improved parallelization -- the rewrite should in particular be able to use more than 4 CPUs;
5. Single-core performance -- the rewrite should not be slower by more than 100% when compared to the reference implementation, on 1 CPU.

Optional requirements
^^^^^^^^^^^^^^^^^^^^^
1. MPI integration;
2. GPU integration;
3. Improving performance -- by using faster equivalent algorithms.

"Human" aspect of the project
-----------------------------

Clients
^^^^^^^
+---------------------+-----------------------+---------------+
| Name                | Occupation            | Role          |
+=====================+=======================+===============+
| Łukasz Mioduszewski | PhD candidate, IF PAN | Contractor    |
+---------------------+-----------------------+---------------+
| Marek Cieplak       | Professor, IF PAN     | Contractor    |
+---------------------+-----------------------+---------------+
| Jacek Sroka         | Researcher, WMIM UW   | Proseminar TA |
+---------------------+-----------------------+---------------+

Users
^^^^^
Users of the product include researchers over at IF PAN and beyond, who may want to be involved in expanding the program in the future, as well as using it for research.

User environment
^^^^^^^^^^^^^^^^
The user should be able to run the program on a standard-issue PC. They should also be capable of modifying the program and releasing the modifications.

Limitations
-----------
Due to the use of floating-point arithmetic and the nondeterminism of the used algorithms, the results of the simulation may not necessarily be identical with the same input. Because of this, the correctness of the program may only be verified statistically.

Milestones
----------

Analysis of the legacy code
^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Deadline**: Dec 3, 2020

A complete description of the reference implementation is created. This description should include:

- Semantics of all procedures, functions and global variables;
- Delineation of the source code into independent modules;
- Profiling the code;
- List of simulation parameters.

Prototype I
^^^^^^^^^^^
**Deadline**: Jan 7, 2021

A program is created, capable of running the first example (ubiquitin, a completely structured protein). Improved parallelization should be apparent even at this stage. Required functionality includes:

- Kernel of the simulation;
- Potentials for structured parts;
- Parsing PDB files, including derivation of contact maps;
- Initial geometry of the model.

Prototype II
^^^^^^^^^^^^
**Deadline**: Feb 25, 2021

A program is extended as to run the second example simulation. Required functionality includes:

- Full output (restarts, statistics etc.);
- PID potential based on a given M-J matrix;
- ATF stretching implemented.

Complete program
^^^^^^^^^^^^^^^^
**Deadline**: May 27, 2021

Entire program, along with the complete documentation and with full test coverage, are created.