=======
Example
=======

.. contents::

Let us take a binary colloidal dispersion of titrating silica particles with radii 8 nm and 12 nm. 
It is in equilibrium with a reservoir of monovalent salt at concentration 10 mM,
and pH 9.5. The number of particles are 200 and 100 respectively.
The system is solved within the Cell Model. The density of the dispersion is fixed by a common
potential at the edge of cells, set here at 0.5 (in kT/e units).
In a directory of the user's choice, a file `distrib.ac` is to be created with the content below:

.. code::

    2
    200 8  0.6
    100 12 0.6

Here an initial surface charge density of 0.6 e/nm^2 is defined.

The input parameters are defined in a file named `polypbren.cfg`.

.. code-block:: ini

    # polypbren.cfg
    model = cellmod             # Cell Model. Use 'jellium' for the Renormalized Jellium Model.
    condition = pH              # charge regulation condition. Use 'charge' for fixed charges 
                                # (set in the distribution file) 'distrib.ac'
    pH_input = 9.5              # the pH
    ion_charges = [1.0, -1.0]    # the charges of ion species in e.
    ion_densities = [10.0, 10.0] # the bulk concentration of each ion species in mM
    temperature = 300.0
    # The following parameters serve for the titration
    site_density = 5.55
    pka = 7.7
    stern_length = 0.107283

    extern_potential = 0.5     # the potential at the edge of the cells in reduced units.
                               # For the RJM, it is the potential at infinity.
    cell_length = 10.0         # the initial length (cell radius - particle radius) in nm. 
                               # Must be higher than the expected final values.


The program has to be run in the same directory with the following command(it may take some time...):

.. code::
    
    polypbren

This should give the output files with the following content:

-   distrib.ac.out:

    .. code::
    
        2
        200 8.0 6.3927895e-01
        100 12.0 6.3102973e-01

-   distrib.out:

    .. code::
    
        2
        200 8.0 2.1222226e-01
        100 12.0 1.9001660e-01

-   peffs.dat:

    .. code::
    
        1.9791311e-01 3.4822501e-01 2.6380472e+00


This example can be retrieved in the `examples/exph9p5bin <https://github.com/guibar64/polypbren/tree/master/examples/exph9p5bin>`_ directory.

The precision of the calculation, e.g. the spatial resolution and the maximum number of iterations,
can also be tuned, see `parameters <usage.html#input-parameters>`_.
