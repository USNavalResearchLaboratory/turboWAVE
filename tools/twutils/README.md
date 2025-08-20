twutils - Python tools for turboWAVE PIC code
=============================================

Modules
-------

1. twutils.NLO: evaluate dispersion in explicit nonlinear optics model
2. twutils.plot: Prepare inputs for matplotlib/mayavi based on data in npy files
3. twutils.pre: Conversion between experimental and simulation parameters
4. twutils.QO.stationary: Stationary states of quantum mechanical wave equations
5. twutils.QO.ionization: Ionization potentials, rates, and units

Scripts
-------

1. twinstall: configures and installs core turboWAVE
2. twplot: command line plotter
3. twscan: serial-parallel or parallel-parallel parameter scans
4. wigner: command line plotter for Wigner distributions
5. twtest: script to run a test suite
6. os2tw: script to convert OSIRIS outputs to turboWAVE outputs

Conda Build
------------

To build and test the conda package locally:
* Install `conda-build` in the base environment
* Edit `setup.py` and `meta.yaml` as appropriate, especially update the version
* From the `tools` directory type `conda build twutils` (may be slow)
* Create directory somewhere: `local-channel/noarch/`
* Copy `miniforge3/conda-bld/noarch/twutils-<x.x.x>-py_0.tar.bz2` to `local-channel/noarch/`
    - on Windows look for it in `$home\AppData\Local`
* `conda index .` from `local-channel`
* `conda create -n tw -c file://full/path/to/local-channel twutils`
