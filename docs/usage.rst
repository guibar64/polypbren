=====
Usage
=====

.. contents::

Command line
============

::
    polypbren [OPTIONS]

Options:

  -h, --help                 Print help message
  -v, --verbosity:N          Set the degree of verbosity. Levels are 0,1,2 (default: 0).
  -d, --distribution:FILE    Change the distribution file to FILE (default: distrib.in)
  -p, --parameters:FILE      Change the distribution file to FILE (default: polypbren.cfg)
  --version                  Print version


Input files
===========

The following input files are needed:

- A distribution file (default name is 'distrib.in', use option ``-d:MyDistribFile`` to set it), 
  The first line is the number of particle types, denoted *n*, 
  then it must be *n* lines of 3 columns: number of particles (integer), radius (nm),
  and surface charge density (e/nm^2). 
  
  The charge is of course unknown in
  the case of pH regulation, it can be considered as a first guess (the better it is the faster the convergence will be)
  .
- A parameter file (default is 'polypbren.cfg', use option ``-p:MyConfig`` to set it).
  It is composed of lines of the form "key = value"
  , array of values are of the form [val1, val2, ...], comments start with '#'.
  Examples are given in the `examples <https://github.com/guibar64/polypbren/tree/master/examples>`_ folder.
  See also the complete list of parameters below.


The output files are:

- distrib.out : gives the effective surface charges densities
  (same formatting as the  distribution files).
- distrib-bare.out : gives the bare surface charge densities when there is titration.
- peffs.dat : gives, in this order, the volume fraction, the effective
  inverse screening length (nm^-1), and the osmotic pressure (Pa).
- phi*p*.dat : gives the electrostatic reduced potential profile of
  particles of type *p*.

Input parameters
================


model
-----
Possible values:
  - *cellmod* for the Cell Model
  - *jellmod* for the Renormalized Jellium Model   

volume_fraction
---------------  
Volume fraction to reach iteratively for the Cell Model

extern_potential
----------------
Potential at the edge of the cell (or at infinity for the Renormalized Jellium Model)

cell_length
-----------
Length of the (initial) cell (- cell radius - particle radius).
For the RJM, while the value is infinite in theory, a finite value should be set.

ion_charges
-----------
List of charge of the ion species. (Up to 8 species are accepted)

ion_densities
-------------
Reservoir of the ion species. The number of items must match the one of IonCharge (Up to 8 species are accepted)

condition
---------
Condition at the surface of the colloid. Set to `pH`, for titrable surface, and `charge` otherwise.

ph_input
--------
Sets the input pH (real number)

temperature
-----------

Sets the temperature

site_density
------------
Sets the density of titrable sites ont the particle surface.

pka
---
Sets the pka of sites.

stern_length
------------
Sets the thickness of the Stern layer (in nm).

space_bin
---------
Sets the space resolution (in nm). Default is 0.05.

max_of_steps
------------
Maximum number of steps to reach volume fraction for the Renormalized Jellium Model or the edge potential for the Cell Model

max_of_iterations
-----------------
Max number of iterations to reach a given pH or surface charge.

max_of_internal_iterations
--------------------------
Max number of iterations of the internal solver (Poisson-Boltzmann
with fixed boundary conditions)

tol_init_value
-------------- 
Sets the tolerance (precision to reach) the titration/charge loop. 

tol_inner_loop
--------------
Sets the tolerance for the internal loop.
  
tol_cell_condition
------------------
Sets the tolerance for the potential at the edge (Cell Model).

tol_volume_fraction
-------------------
Sets the tolerance for the potential at the edge (Renormalized Jellium Model).

max_potential_zeff_jellium
--------------------------
For the RJM, sets the maximum potential (in absolute value) at which point a fit is performed to compute the renormalized charges.
Default is 0.5.
