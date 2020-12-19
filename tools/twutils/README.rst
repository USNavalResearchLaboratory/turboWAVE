twutils - Python tools for turboWAVE PIC code
=============================================

Modules
-------

#. twutils.NLO: evaluate dispersion in explicit nonlinear optics model
#. twutils.plot: Prepare inputs for matplotlib/mayavi based on data in npy files
#. twutils.pre: Conversion between experimental and simulation parameters
#. twutils.QO.stationary: Stationary states of quantum mechanical wave equations
#. twutils.QO.ionization: Ionization potentials, rates, and units

Scripts
-------

#. twinstall: configures and installs core turboWAVE
#. twplot: command line plotter
#. twscan: serial-parallel or parallel-parallel parameter scans
#. wigner: command line plotter for Wigner distributions
#. twtest: script to run a test suite
#. os2tw: script to convert OSIRIS outputs to turboWAVE outputs

Known Issues
------------

The twinstall script depends on git.  The pypi package cannot account for this, while the conda package can.  So if one uses pip rather than conda, git must be installed beforehand.
