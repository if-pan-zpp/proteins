Overview
========

High-level control flow
-----------------------
The structure of the simulation at the highest level is similar to other MD programs:

1. Initialize variables, fields, read protein files etc.
2. For every "sample" (in this program we loop over temperatures and trajectories):
   
  a. Reset geometry and other parameters;
  b. Do "equilibration" (simulation without walls and external forces);
  c. Simulate the system proper: compute forces and update the state with an integrator, in a loop.

(These steps, of course, are interrupted by reporting, saving restart files and current geometry etc.)

Geometry of the system
----------------------
By default, the simulation occurs in an infinite space. This we can change by introducing a "simulation box"; there are, for each pair of opposite faces of the box (X, Y, Z), following options:

- there is a "solid" (repulsive) wall of potential;
- there is a wall in the form of an immovable face-centered-cubic lattice of beads, which interacts quasi-adiabatically with the chains;
- there is a flat attractive wall, which is the same as above but without the beads;
- there is no wall, and in particular we should treat it as a periodic boundary.

The size of the box is determined either by a ``CRYST1`` record or derived automatically as to achieve a particular density (the box is a cube, then).

Moreover, we can "program" the walls to change size, in particular oscillate.

Types of proteins
-----------------
The program accepts two types of proteins:

- structured -- ones where the native structure in the form of xyz coordinates is already present -- in this case, positions of CA atoms, crystallographic data (including repetition of chains) and disulfide bonds are extracted; contact map is either derived or given by a provided file;
- (partially) unstructured -- where we are given chains of residues (but without the coordinates), possibly along with some predefined contacts.

The force fields that may appear in the simulation depend on the type of protein.

Force fields
------------