Input Files
===========

Everything about a turboWAVE simulation is determined by the input file.  The input file can have any name, though it is best to not use spaces.  The ``.tw`` extension can be made to trigger syntax highlights in certain editors.  If no input file name is specified turboWAVE assumes it is ``stdin``.

The turboWAVE input file is written in what amounts to a "little language", which is described further in :doc:`ref-input`.  Fortunately this language is very simple and does not impose a significant learning curve.

This tutorial discusses simulation units, and then points you to the example files and the reference pages.

Simulation Units
----------------

There is a normalization scheme used throughout most of turboWAVE (quantum modules use atomic or natural units).  Most of the input file parameters are given in the normalized units, although some conventional units can be used by means of simple macros.  Output files are written in normalized units without exception.

The normalization scheme can be thought of in many ways.  The fundamental observation is that any solution of the Vlasov equation can be scaled up or down to produce a family of solutions.  We may as well express all quantities in a way that does not commit to which particular member of the family we are talking about.  No particular scale is any better than another.

On the other hand, if an atomic process like ionization comes into the problem, then we must commit to a definite physical scale.  For this reason, turboWAVE input files allow you to specify a unit of density that fixes the physical scale.

If you prefer to think in terms of physical scales, then specify a unit density.  To normalize some quantity given in physical units, divide by the unit given in Table I.

.. csv-table:: Table I. Simulation Units to Gaussian or SI Units
	:header: "Quantity", "Unit", "cgs Symbol", "mks Symbol"

	"Number Density", "user specified", :math:`n_1`, :math:`n_1`
	"Velocity", "speed of light", "c", "c"
	"Charge", "electronic charge", "e", "e"
	"Mass", "electronic mass", "m", "m"
	"Time", "derived", :math:`\omega_p^{-1}=\sqrt{\frac{m}{4\pi n_1 e^2}}`, :math:`\omega_p^{-1}=\sqrt{\frac{\epsilon_0 m}{n_1 e^2}}`
	"Space", "derived", :math:`c/\omega_p`, :math:`c/\omega_p`
	"Energy", "derived", :math:`mc^2`, :math:`mc^2`
	"Momentum", "derived", :math:`mc`, :math:`mc`
	"Scalar Potential", "derived", :math:`mc^2/e`, :math:`mc^2/e`
	"Vector Potential", "derived", :math:`mc^2/e`, :math:`mc/e`
	"Electric Field", "derived", :math:`mc\omega_p/e`, :math:`mc\omega_p/e`
	"Magnetic Field", "derived", :math:`mc\omega_p/e`, :math:`m\omega_p/e`
	"Particle Number","derived", :math:`n_1(c/\omega_p)^3`, :math:`n_1(c/\omega_p)^3`

.. tip::

	If you have a particular radiation frequency that is of interest, setting the unit of density to the critical density for that radiation will make the radiation angular frequency unity in simulation units.

.. tip::

	You can make the unit of length correspond to a conventional unit by choosing the unit of density such that :math:`c/\omega_p` comes out to your chosen unit (e.g., 1 cm).  You can do something similar with the unit of time.  However, the unit of time and space cannot both be nice decimal numbers in conventional units, since the speed of light is not.

On the Unit of Particle Number
-------------------------------

One thing that can be confusing in the normalization scheme is that the unit of density is not connected to the unit of length in the expected way.  This comes about because we are simultaneously demanding that the speed of light be unity, and that the unit of time be connected with the density.  This can lead to erroneous normalization of certain quantities with a density buried in them.  As an example, take current density.  The naive construction might be to take the unit as :math:`(c/\omega_p)^{-3}ec` (wrong).  The correct construction is :math:`n_1ec`.

It can be useful to identify, in this connection, a unit of particle number, per Table I.  Take the unit of energy.  The normalization in table I refers to the energy of a real particle "within" a simulation macroparticle.  If we want the energy contained in the whole macroparticle, we must multiply by the unit of particle number.  In fact, any volume integrated energy that turboWAVE writes out must include this factor if you want to convert to conventional units.

Unit Conversion Macros
-----------------------

In order to allow the use of conventional units in the input file, there is a simple notation that causes the parser to automatically perform a conversion.  The macro is formed by prepending a decimal value with the :samp:`%` symbol, and postpending it with a string identifying the unit.  As an example, enter :samp:`%10ps` to indicate a value of 10 picoseconds.  During pre-processing of the input file, the macro will be replaced with the decimal value of the quantity in normalized units.  For the full list of available conversions see :ref:`unit-conv`.

TurboWAVE Example Files
-----------------------

You will learn the most by studying the examples.  They can be found in :samp:`{twroot}/core/examples` and are organized into several directories.  These are as follows.

#. :samp:`hydro`: Contains examples of hydrodynamic simulations that use the SPARC module.  If you plan to use turboWAVE for hydro simulations, you should especially understand the simple-shock examples, and the diffusion example.  The shock cases can be compared with an analytical theory.  There is a Mathematica notebook in the folder which can be used to explore the analytical solution.
#. :samp:`pic`: Contains examples of fully explicit PIC simulations.  This contains especially variants on laser driven wakefields.
#. :samp:`pgc`: Similar to the :samp:`pic` directory, except uses the ponderomotive guiding center approximation to model the laser fields.
#. :samp:`fluid`: Cold relativistic fluid approximation for laser wakefield and beatwave cases.
#. :samp:`nonlinear-optics`: Contains examples of the nonlinear optics model for laser radiation in crystals.
#. :samp:`quantum`: Contains examples of atomic level processes using the quantum optoelectronics modules.
#. :samp:`misc`: Some other examples.


TurboWAVE Reference
-------------------

Detailed exposition of the input file elements are in the reference pages.  Highlights include:

	* :doc:`ref-input` describes the input file "little language"
	* :doc:`ref-tools` describes some important tools that can be used in various modules
	* :doc:`ref-PIC` describes modules used for PIC simulations
	* :doc:`ref-fluid` describes modules used for fluid simulations
