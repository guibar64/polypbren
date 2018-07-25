========
Overview
========


Polypbren is a program that computes the renormalized parameters (renormalized charges and screening length) for dispersions of
charged spherical particles with arbitrary size distributions. The surface charge of the colloids can be fixed or regulated by pH.
The Poisson-Boltzmann Cell Model and Renormalized Jellium Model are implemented. See the `References`_ for more information on the scientific part.


.. contents::


Building
========

You will need the `Nim compiler <https://nim-lang.org/>`_ to compile this program. 
The installation instructions are `here <https://nim-lang.org/install.html>`_. 


The following command line

::
   nimble install https://github.com/guibar64/polypbren.git

will download the source files, compile an executable (``polypbren.exe`` on windows, ``polypbren`` elsewhere). 
The binary is installed by default in ``~/.nimble/bin``
(or elsewhere depending on your nimble setting and OS).

Alternatively the whole repository (including examples) can be downloaded,
and the main program can be manually compiled following the command lines,

::
   
    git clone https://github.com/guibar64/polypbren.git
    cd polypbren
    nim c -d:release src/polypbren 

The generated executable is in the  ``src`` directory.
The ``-d:release`` flag is for a release build, skip it if you want stack traces, checks etc... but it is then much slower.

Usage
=====

Assuming `polypbren` is in your path,

::
    polypbren -h

will print an help message.

See `Usage <usage.html>`_ for more information on input files and parameters.

Solver
======

The module `pbsolv <pbsolv.html>`_ implements a simple solver for the 1-D Poisson-Boltzmann equation.

References
==========

S. Alexander, P. M. Chaikin, P. Grant, G. J. Morales, P. Pincus, and D. Hone, “Charge
renormalization, osmotic pressure, and bulk modulus of colloidal crystals: Theory,” The
Journal of Chemical Physics 80, 5776–5781 (1984), `<http://scitation.aip.org/content/aip/journal/jcp/80/11/10.1063/1.446600>`_

A. Torres, G. Téllez, and R. van Roij, “The polydisperse cell model: Nonlinear screening
and charge renormalization in colloidal mixtures,” Journal of Chemical Physics, Volume
128, Issue 15, pp. 154906-154906-8 (2008), `<http://adsabs.harvard.edu/abs/2008JChPh.128o4906T>`_


E. Trizac and Y. Levin, “Renormalized jellium model for charge-stabilized colloidal sus-
pensions,” Phys. Rev. E 69, 031403 (2004),
`<http://link.aps.org/doi/10.1103/PhysRevE.69.031403>`_


J. M. Falcón-González and R. Castañeda-Priego, “Renormalized jellium mean-field ap-
proximation for binary mixtures of charged colloids,” Phys. Rev. E 83, 041401 (2011), 
`<http://link.aps.org/doi/10.1103/PhysRevE.83.041401>`_

License
=======

This project is licensed under the terms of the MIT license.
