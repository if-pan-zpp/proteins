OpenMM implementation
=====================
Before implementing a custom architecture for the model, we have decided to investigate whether one could utilize existing tools, such as `OpenMM <http://openmm.org/>`_, for implementing it.

Advantages
----------
1. OpenMM includes MPI and GPU platforms by default, so one could achieve massive performance gains with it;
2. A vast array of existing force fields and other tools, such as different integrators, handling PBC, PDB file parsing, automatic Verlet/cell lists etc. are available by default;
3. There are many frontends for it, including Python and C++.

Challenges
----------
1. OpenMM is designed with all-atom simulations in mind -- in particular, doing coarse-grained simulations efficiently may prove to be non-trivial, especially when we wish to use other tools in the suite;
2. Walls, especially resizable wall and lattice-like ones, are not available by default -- there are, however, facilities for defining an appropriate force field;
3. Quasi-adiabatic and PID potentials are not included, and in particular there is no clear way of defining custom force fields based on *dynamically* created/destroyed contacts.
